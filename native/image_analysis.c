/*
 * Native image analysis, from pixels to focus position.
 *
 * The public workflow is deliberately staged:
 *   1. IA_AnalyzeFrame estimates the global background and finds stars.
 *   2. IA_ComputeStarROI chooses a sensor rectangle around selected stars.
 *   3. IA_MeasureLockedStars follows those same stars through a focus sweep.
 *   4. A stateful focus run fits intersecting branches and combines stars.
 *
 * Callers own both image memory and workspace memory.  This keeps allocation
 * and hardware policy in Forth while this file remains a deterministic
 * numerical engine.  Integer milli-units make the DLL ABI reproducible across
 * compilers; floating point is used only inside calculations.
 */

#include "image_analysis.h"

#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define IA_HISTOGRAM_BINS 65536u
#define IA_CANDIDATES_PER_STAR 8u
#define IA_RADIAL_BINS 1024u
#define IA_RADIAL_SCALE 16.0

typedef struct IACandidate {
    uint32_t x;
    uint32_t y;
    uint32_t peak;
    uint32_t area;
} IACandidate;

typedef struct IAWorkspaceView {
    uint32_t *histogram;
    IACandidate *candidates;
    IACandidate *selected;
    uint32_t candidate_capacity;
} IAWorkspaceView;

static int ia_debug_enabled;

static void ia_debug_log(const char *format, ...)
{
    FILE *file;
    time_t now;
    struct tm local;
    char timestamp[24];
    va_list arguments;

    if (!ia_debug_enabled) {
        return;
    }
#ifdef _WIN32
    if (fopen_s(&file, IA_DEBUG_LOG_PATH, "a") != 0) {
        return;
    }
#else
    file = fopen(IA_DEBUG_LOG_PATH, "a");
    if (file == NULL) {
        return;
    }
#endif
    now = time(NULL);
#ifdef _WIN32
    localtime_s(&local, &now);
#else
    localtime_r(&now, &local);
#endif
    strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", &local);
    fprintf(file, "%s ", timestamp);
    va_start(arguments, format);
    vfprintf(file, format, arguments);
    va_end(arguments);
    fputc('\n', file);
    fflush(file);
    fclose(file);
}

static uint32_t ia_clamp_u32(double value)
{
    if (value <= 0.0) {
        return 0;
    }
    if (value >= 4294967295.0) {
        return UINT32_MAX;
    }
    return (uint32_t)(value + 0.5);
}

static int32_t ia_clamp_i32(double value)
{
    if (value <= -2147483648.0) {
        return INT32_MIN;
    }
    if (value >= 2147483647.0) {
        return INT32_MAX;
    }
    return (int32_t)(value >= 0.0 ? value + 0.5 : value - 0.5);
}

static uint32_t ia_align_u32(uint32_t value, uint32_t alignment)
{
    return (value + alignment - 1u) & ~(alignment - 1u);
}

static uint32_t ia_workspace_bytes(uint32_t maximum_stars)
{
    uint64_t candidate_capacity;
    uint64_t bytes;

    if (maximum_stars == 0u) {
        return 0u;
    }
    candidate_capacity = (uint64_t)maximum_stars * IA_CANDIDATES_PER_STAR;
    bytes = (uint64_t)IA_HISTOGRAM_BINS * sizeof(uint32_t);
    bytes = ia_align_u32((uint32_t)bytes, (uint32_t)_Alignof(IACandidate));
    bytes += candidate_capacity * sizeof(IACandidate);
    bytes += (uint64_t)maximum_stars * sizeof(IACandidate);
    if (bytes > UINT32_MAX) {
        return 0u;
    }
    return (uint32_t)bytes;
}

static int ia_workspace_view(
    void *workspace,
    uint32_t workspace_bytes,
    uint32_t maximum_stars,
    IAWorkspaceView *view)
{
    uint8_t *cursor;
    uint32_t required;
    uint32_t histogram_bytes;

    required = ia_workspace_bytes(maximum_stars);
    if (required == 0u || workspace == NULL || workspace_bytes < required) {
        return IA_WORKSPACE_TOO_SMALL;
    }

    cursor = (uint8_t *)workspace;
    histogram_bytes = IA_HISTOGRAM_BINS * (uint32_t)sizeof(uint32_t);
    view->histogram = (uint32_t *)cursor;
    cursor += ia_align_u32(histogram_bytes, (uint32_t)_Alignof(IACandidate));
    view->candidate_capacity = maximum_stars * IA_CANDIDATES_PER_STAR;
    view->candidates = (IACandidate *)cursor;
    cursor += view->candidate_capacity * sizeof(IACandidate);
    view->selected = (IACandidate *)cursor;
    return IA_OK;
}

static uint32_t ia_histogram_quantile(
    const uint32_t *histogram,
    uint64_t sample_count,
    uint64_t numerator,
    uint64_t denominator)
{
    uint64_t target;
    uint64_t cumulative;
    uint32_t i;

    if (sample_count == 0u || denominator == 0u) {
        return 0u;
    }
    target = (sample_count * numerator) / denominator;
    cumulative = 0u;
    for (i = 0u; i < IA_HISTOGRAM_BINS; ++i) {
        cumulative += histogram[i];
        if (cumulative > target) {
            return i;
        }
    }
    return IA_HISTOGRAM_BINS - 1u;
}

static void ia_global_background(
    const uint16_t *pixels,
    uint32_t width,
    uint32_t height,
    uint32_t stride,
    uint32_t *histogram,
    uint32_t *median,
    double *noise)
{
    uint64_t count;
    uint32_t y;
    uint32_t x;

    memset(histogram, 0, IA_HISTOGRAM_BINS * sizeof(uint32_t));
    for (y = 0u; y < height; ++y) {
        const uint16_t *row = pixels + (size_t)y * stride;
        for (x = 0u; x < width; ++x) {
            ++histogram[row[x]];
        }
    }

    count = (uint64_t)width * height;
    *median = ia_histogram_quantile(histogram, count, 1u, 2u);
    memset(histogram, 0, IA_HISTOGRAM_BINS * sizeof(uint32_t));
    for (y = 0u; y < height; ++y) {
        const uint16_t *row = pixels + (size_t)y * stride;
        for (x = 0u; x < width; ++x) {
            uint32_t value = row[x];
            uint32_t deviation = value >= *median ? value - *median : *median - value;
            ++histogram[deviation];
        }
    }
    *noise = 1.4826 * ia_histogram_quantile(histogram, count, 1u, 2u);
}

static int ia_candidate_weaker(const IACandidate *left, const IACandidate *right)
{
    if (left->peak != right->peak) {
        return left->peak < right->peak;
    }
    if (left->area != right->area) {
        return left->area < right->area;
    }
    if (left->y != right->y) {
        return left->y > right->y;
    }
    return left->x > right->x;
}

static void ia_heap_sift_up(IACandidate *heap, uint32_t index)
{
    while (index > 0u) {
        uint32_t parent = (index - 1u) / 2u;
        if (!ia_candidate_weaker(&heap[index], &heap[parent])) {
            break;
        }
        {
            IACandidate temporary = heap[parent];
            heap[parent] = heap[index];
            heap[index] = temporary;
        }
        index = parent;
    }
}

static void ia_heap_sift_down(IACandidate *heap, uint32_t count, uint32_t index)
{
    for (;;) {
        uint32_t left = index * 2u + 1u;
        uint32_t right = left + 1u;
        uint32_t weakest = index;

        if (left < count && ia_candidate_weaker(&heap[left], &heap[weakest])) {
            weakest = left;
        }
        if (right < count && ia_candidate_weaker(&heap[right], &heap[weakest])) {
            weakest = right;
        }
        if (weakest == index) {
            break;
        }
        {
            IACandidate temporary = heap[index];
            heap[index] = heap[weakest];
            heap[weakest] = temporary;
        }
        index = weakest;
    }
}

static void ia_heap_add(
    IACandidate *heap,
    uint32_t capacity,
    uint32_t *count,
    IACandidate candidate)
{
    if (*count < capacity) {
        heap[*count] = candidate;
        ia_heap_sift_up(heap, *count);
        ++*count;
        return;
    }
    if (ia_candidate_weaker(&heap[0], &candidate)) {
        heap[0] = candidate;
        ia_heap_sift_down(heap, *count, 0u);
    }
}

static int ia_candidate_descending(const void *left_pointer, const void *right_pointer)
{
    const IACandidate *left = (const IACandidate *)left_pointer;
    const IACandidate *right = (const IACandidate *)right_pointer;

    if (left->peak != right->peak) {
        return left->peak < right->peak ? 1 : -1;
    }
    if (left->area != right->area) {
        return left->area < right->area ? 1 : -1;
    }
    if (left->y != right->y) {
        return left->y > right->y ? 1 : -1;
    }
    if (left->x != right->x) {
        return left->x > right->x ? 1 : -1;
    }
    return 0;
}

static uint32_t ia_detect_candidates(
    const uint16_t *pixels,
    uint32_t width,
    uint32_t height,
    uint32_t stride,
    uint32_t threshold,
    const IAConfig *config,
    IACandidate *candidates,
    uint32_t candidate_capacity)
{
    uint32_t count = 0u;
    uint32_t y;
    uint32_t x;

    if (width < 5u || height < 5u) {
        return 0u;
    }

    for (y = 2u; y + 2u < height; ++y) {
        const uint16_t *row = pixels + (size_t)y * stride;
        for (x = 2u; x + 2u < width; ++x) {
            uint32_t value = row[x];
            uint32_t area = 0u;
            int strict = 0;
            int local_maximum = 1;
            int32_t dy;
            int32_t dx;

            if (value <= threshold) {
                continue;
            }
            for (dy = -1; dy <= 1 && local_maximum; ++dy) {
                const uint16_t *near_row = pixels + (size_t)(y + dy) * stride;
                for (dx = -1; dx <= 1; ++dx) {
                    uint32_t neighbour;
                    if (dx == 0 && dy == 0) {
                        continue;
                    }
                    neighbour = near_row[x + dx];
                    if (neighbour > value) {
                        local_maximum = 0;
                        break;
                    }
                    if (neighbour < value) {
                        strict = 1;
                    }
                }
            }
            if (!local_maximum || !strict) {
                continue;
            }
            for (dy = -2; dy <= 2; ++dy) {
                const uint16_t *near_row = pixels + (size_t)(y + dy) * stride;
                for (dx = -2; dx <= 2; ++dx) {
                    if (near_row[x + dx] > threshold) {
                        ++area;
                    }
                }
            }
            if (area < config->minimum_area_pixels) {
                continue;
            }
            {
                IACandidate candidate;
                candidate.x = x;
                candidate.y = y;
                candidate.peak = value;
                candidate.area = area;
                ia_heap_add(candidates, candidate_capacity, &count, candidate);
            }
        }
    }

    qsort(candidates, count, sizeof(IACandidate), ia_candidate_descending);
    return count;
}

static uint32_t ia_select_separated(
    const IACandidate *candidates,
    uint32_t candidate_count,
    uint32_t minimum_separation,
    IACandidate *selected,
    uint32_t selected_capacity)
{
    uint64_t separation_squared = (uint64_t)minimum_separation * minimum_separation;
    uint32_t selected_count = 0u;
    uint32_t i;

    for (i = 0u; i < candidate_count && selected_count < selected_capacity; ++i) {
        uint32_t j;
        int separated = 1;
        for (j = 0u; j < selected_count; ++j) {
            int64_t dx = (int64_t)candidates[i].x - selected[j].x;
            int64_t dy = (int64_t)candidates[i].y - selected[j].y;
            uint64_t distance_squared = (uint64_t)(dx * dx + dy * dy);
            if (distance_squared < separation_squared) {
                separated = 0;
                break;
            }
        }
        if (separated) {
            selected[selected_count++] = candidates[i];
        }
    }
    return selected_count;
}

static void ia_local_background(
    const uint16_t *pixels,
    uint32_t width,
    uint32_t height,
    uint32_t stride,
    double center_x,
    double center_y,
    uint32_t inner_radius,
    uint32_t outer_radius,
    double *background,
    double *noise)
{
    int32_t left = (int32_t)floor(center_x) - (int32_t)outer_radius;
    int32_t right = (int32_t)ceil(center_x) + (int32_t)outer_radius;
    int32_t top = (int32_t)floor(center_y) - (int32_t)outer_radius;
    int32_t bottom = (int32_t)ceil(center_y) + (int32_t)outer_radius;
    double inner_squared = (double)inner_radius * inner_radius;
    double outer_squared = (double)outer_radius * outer_radius;
    double sum = 0.0;
    double sum_squared = 0.0;
    uint32_t count = 0u;
    int32_t y;
    int32_t x;

    if (left < 0) left = 0;
    if (top < 0) top = 0;
    if (right >= (int32_t)width) right = (int32_t)width - 1;
    if (bottom >= (int32_t)height) bottom = (int32_t)height - 1;

    for (y = top; y <= bottom; ++y) {
        const uint16_t *row = pixels + (size_t)y * stride;
        for (x = left; x <= right; ++x) {
            double dx = x - center_x;
            double dy = y - center_y;
            double radius_squared = dx * dx + dy * dy;
            if (radius_squared >= inner_squared && radius_squared <= outer_squared) {
                double value = row[x];
                sum += value;
                sum_squared += value * value;
                ++count;
            }
        }
    }

    if (count == 0u) {
        *background = 0.0;
        *noise = 0.0;
        return;
    }

    *background = sum / count;
    *noise = sqrt(fmax(0.0, sum_squared / count - *background * *background));
    if (*noise > 0.0) {
        double limit = *background + 3.0 * *noise;
        sum = 0.0;
        sum_squared = 0.0;
        count = 0u;
        for (y = top; y <= bottom; ++y) {
            const uint16_t *row = pixels + (size_t)y * stride;
            for (x = left; x <= right; ++x) {
                double dx = x - center_x;
                double dy = y - center_y;
                double radius_squared = dx * dx + dy * dy;
                if (radius_squared >= inner_squared && radius_squared <= outer_squared &&
                    row[x] <= limit) {
                    double value = row[x];
                    sum += value;
                    sum_squared += value * value;
                    ++count;
                }
            }
        }
        if (count > 0u) {
            *background = sum / count;
            *noise = sqrt(fmax(0.0, sum_squared / count - *background * *background));
        }
    }
}

static void ia_measure_candidate(
    const uint16_t *pixels,
    uint32_t width,
    uint32_t height,
    uint32_t stride,
    const IAConfig *config,
    const IACandidate *candidate,
    const IACandidate *all_candidates,
    uint32_t candidate_count,
    IAStar *star)
{
    double center_x = candidate->x;
    double center_y = candidate->y;
    double background = 0.0;
    double noise = 0.0;
    double flux = 0.0;
    double weighted_x = 0.0;
    double weighted_y = 0.0;
    double moment_xx = 0.0;
    double moment_xy = 0.0;
    double moment_yy = 0.0;
    double radial_flux[IA_RADIAL_BINS];
    double aperture_squared;
    uint32_t radius = config->aperture_radius_pixels;
    uint32_t positive_pixels = 0u;
    uint32_t saturated_pixels = 0u;
    uint32_t peak = 0u;
    int32_t left;
    int32_t right;
    int32_t top;
    int32_t bottom;
    int32_t y;
    int32_t x;
    uint32_t i;

    memset(star, 0, sizeof(*star));
    star->struct_size = sizeof(*star);

    ia_local_background(
        pixels, width, height, stride, center_x, center_y,
        config->annulus_inner_radius_pixels,
        config->annulus_outer_radius_pixels,
        &background, &noise);

    aperture_squared = (double)radius * radius;
    left = (int32_t)candidate->x - (int32_t)radius;
    right = (int32_t)candidate->x + (int32_t)radius;
    top = (int32_t)candidate->y - (int32_t)radius;
    bottom = (int32_t)candidate->y + (int32_t)radius;
    if (left < 0) left = 0;
    if (top < 0) top = 0;
    if (right >= (int32_t)width) right = (int32_t)width - 1;
    if (bottom >= (int32_t)height) bottom = (int32_t)height - 1;

    for (y = top; y <= bottom; ++y) {
        const uint16_t *row = pixels + (size_t)y * stride;
        for (x = left; x <= right; ++x) {
            double dx = x - center_x;
            double dy = y - center_y;
            if (dx * dx + dy * dy <= aperture_squared) {
                double signal = row[x] - background;
                if (signal > 0.0) {
                    flux += signal;
                    weighted_x += signal * x;
                    weighted_y += signal * y;
                }
            }
        }
    }
    if (flux > 0.0) {
        center_x = weighted_x / flux;
        center_y = weighted_y / flux;
    }

    ia_local_background(
        pixels, width, height, stride, center_x, center_y,
        config->annulus_inner_radius_pixels,
        config->annulus_outer_radius_pixels,
        &background, &noise);

    left = (int32_t)floor(center_x) - (int32_t)radius;
    right = (int32_t)ceil(center_x) + (int32_t)radius;
    top = (int32_t)floor(center_y) - (int32_t)radius;
    bottom = (int32_t)ceil(center_y) + (int32_t)radius;
    if (left < 0 || top < 0 || right >= (int32_t)width || bottom >= (int32_t)height) {
        star->flags |= IA_STAR_EDGE;
    }
    if (left < 0) left = 0;
    if (top < 0) top = 0;
    if (right >= (int32_t)width) right = (int32_t)width - 1;
    if (bottom >= (int32_t)height) bottom = (int32_t)height - 1;
    star->left = left;
    star->top = top;
    star->right = right;
    star->bottom = bottom;

    memset(radial_flux, 0, sizeof(radial_flux));
    flux = 0.0;
    for (y = top; y <= bottom; ++y) {
        const uint16_t *row = pixels + (size_t)y * stride;
        for (x = left; x <= right; ++x) {
            double dx = x - center_x;
            double dy = y - center_y;
            double radius_squared = dx * dx + dy * dy;
            if (radius_squared <= aperture_squared) {
                uint32_t value = row[x];
                double signal = value - background;
                if (value > peak) {
                    peak = value;
                }
                if (value >= config->saturation_adu) {
                    ++saturated_pixels;
                }
                if (signal > 0.0) {
                    uint32_t radial_bin;
                    double radial_distance = sqrt(radius_squared);
                    flux += signal;
                    moment_xx += signal * dx * dx;
                    moment_xy += signal * dx * dy;
                    moment_yy += signal * dy * dy;
                    ++positive_pixels;
                    radial_bin = (uint32_t)(radial_distance * IA_RADIAL_SCALE);
                    if (radial_bin >= IA_RADIAL_BINS) {
                        radial_bin = IA_RADIAL_BINS - 1u;
                    }
                    radial_flux[radial_bin] += signal;
                }
            }
        }
    }

    star->x_milli = ia_clamp_i32(center_x * IA_MILLI);
    star->y_milli = ia_clamp_i32(center_y * IA_MILLI);
    star->flux_adu = ia_clamp_u32(flux);
    star->peak_adu = peak;
    star->saturated_pixels = saturated_pixels;
    star->background_milli_adu = ia_clamp_i32(background * IA_MILLI);
    star->noise_milli_adu = ia_clamp_i32(noise * IA_MILLI);
    star->positive_pixels = positive_pixels;

    if (saturated_pixels > 0u) {
        star->flags |= IA_STAR_SATURATED;
    }
    if (positive_pixels < config->minimum_area_pixels || flux <= 0.0) {
        star->flags |= IA_STAR_TOO_SMALL;
    }

    if (flux > 0.0) {
        double half_flux = flux * 0.5;
        double cumulative = 0.0;
        double hfr = 0.0;
        double trace;
        double difference;
        double discriminant;
        double lambda_major;
        double lambda_minor;
        double sigma_major;
        double sigma_minor;
        double fwhm;
        double ellipticity;
        double snr;

        for (i = 0u; i < IA_RADIAL_BINS; ++i) {
            double previous = cumulative;
            cumulative += radial_flux[i];
            if (cumulative >= half_flux) {
                double fraction = radial_flux[i] > 0.0 ?
                    (half_flux - previous) / radial_flux[i] : 0.0;
                hfr = ((double)i + fraction) / IA_RADIAL_SCALE;
                break;
            }
        }

        moment_xx /= flux;
        moment_xy /= flux;
        moment_yy /= flux;
        trace = moment_xx + moment_yy;
        difference = moment_xx - moment_yy;
        discriminant = sqrt(fmax(0.0, difference * difference + 4.0 * moment_xy * moment_xy));
        lambda_major = fmax(0.0, (trace + discriminant) * 0.5);
        lambda_minor = fmax(0.0, (trace - discriminant) * 0.5);
        sigma_major = sqrt(lambda_major);
        sigma_minor = sqrt(lambda_minor);
        fwhm = 2.354820045 * sqrt(fmax(0.0, trace * 0.5));
        ellipticity = sigma_major > 0.0 ? 1.0 - sigma_minor / sigma_major : 0.0;
        snr = noise > 0.0 ? flux / (noise * sqrt((double)positive_pixels)) : flux;

        star->hfr_milli_pixels = ia_clamp_i32(hfr * IA_MILLI);
        star->hfd_milli_pixels = ia_clamp_i32(2.0 * hfr * IA_MILLI);
        star->fwhm_milli_pixels = ia_clamp_i32(fwhm * IA_MILLI);
        star->sigma_major_milli_pixels = ia_clamp_i32(sigma_major * IA_MILLI);
        star->sigma_minor_milli_pixels = ia_clamp_i32(sigma_minor * IA_MILLI);
        star->ellipticity_milli = ia_clamp_i32(ellipticity * IA_MILLI);
        star->theta_millidegrees = ia_clamp_i32(
            0.5 * atan2(2.0 * moment_xy, difference) * (180000.0 / 3.14159265358979323846));
        star->snr_milli = ia_clamp_i32(snr * IA_MILLI);
        if (star->snr_milli < (int32_t)config->minimum_snr_milli) {
            star->flags |= IA_STAR_LOW_SNR;
        }
    } else {
        star->flags |= IA_STAR_LOW_SNR;
    }

    for (i = 0u; i < candidate_count; ++i) {
        int64_t dx;
        int64_t dy;
        uint64_t distance_squared;
        if (all_candidates[i].x == candidate->x && all_candidates[i].y == candidate->y) {
            continue;
        }
        dx = (int64_t)all_candidates[i].x - candidate->x;
        dy = (int64_t)all_candidates[i].y - candidate->y;
        distance_squared = (uint64_t)(dx * dx + dy * dy);
        if (distance_squared <= (uint64_t)radius * radius &&
            distance_squared >=
                (uint64_t)config->minimum_separation_pixels *
                    config->minimum_separation_pixels &&
            (uint64_t)all_candidates[i].peak * IA_MILLI >=
                (uint64_t)candidate->peak *
                    config->blend_minimum_contrast_milli) {
            star->flags |= IA_STAR_BLENDED;
            break;
        }
    }
}

static int ia_int32_ascending(const void *left, const void *right)
{
    int32_t a = *(const int32_t *)left;
    int32_t b = *(const int32_t *)right;
    return (a > b) - (a < b);
}

static int32_t ia_median_i32(int32_t *values, uint32_t count)
{
    if (count == 0u) {
        return 0;
    }
    qsort(values, count, sizeof(int32_t), ia_int32_ascending);
    if ((count & 1u) != 0u) {
        return values[count / 2u];
    }
    return (int32_t)(((int64_t)values[count / 2u - 1u] + values[count / 2u]) / 2);
}

static int32_t ia_mad_i32(int32_t *values, uint32_t count, int32_t median)
{
    uint32_t i;
    for (i = 0u; i < count; ++i) {
        int64_t difference = (int64_t)values[i] - median;
        if (difference < 0) {
            difference = -difference;
        }
        values[i] = difference > INT32_MAX ? INT32_MAX : (int32_t)difference;
    }
    return ia_median_i32(values, count);
}

uint32_t IA_CALL IA_Version(void)
{
    ia_debug_log("IA_Version version=%u", IA_API_VERSION);
    return IA_API_VERSION;
}

void IA_CALL IA_SetDebug(uint32_t enabled)
{
    if (enabled != 0u) {
        ia_debug_enabled = 1;
        ia_debug_log("debug enabled");
    } else {
        ia_debug_log("debug disabled");
        ia_debug_enabled = 0;
    }
}

void IA_CALL IA_DefaultConfig(IAConfig *config)
{
    if (config == NULL) {
        return;
    }
    memset(config, 0, sizeof(*config));
    config->struct_size = sizeof(*config);
    config->detection_sigma_milli = 5000u;
    config->minimum_threshold_adu = 16u;
    config->saturation_adu = 65535u;
    config->minimum_area_pixels = 3u;
    config->minimum_separation_pixels = 8u;
    config->aperture_radius_pixels = 16u;
    config->annulus_inner_radius_pixels = 18u;
    config->annulus_outer_radius_pixels = 24u;
    config->minimum_snr_milli = 5000u;
    config->blend_minimum_contrast_milli = 500u;
}

uint32_t IA_CALL IA_WorkspaceBytes(uint32_t maximum_stars)
{
    return ia_workspace_bytes(maximum_stars);
}

int IA_CALL IA_SummarizeStars(
    const IAStar *stars,
    uint32_t star_count,
    uint32_t reject_flags,
    void *workspace,
    uint32_t workspace_bytes,
    IAFrameSummary *summary)
{
    int32_t *values;
    uint32_t required;
    uint32_t included = 0u;
    uint32_t i;

    if (stars == NULL || summary == NULL || summary->struct_size < sizeof(*summary)) {
        return IA_INVALID_ARGUMENT;
    }
    required = star_count * (uint32_t)sizeof(int32_t);
    if (star_count > 0u && (workspace == NULL || workspace_bytes < required)) {
        return IA_WORKSPACE_TOO_SMALL;
    }

    summary->measured_stars = star_count;
    summary->included_stars = 0u;
    summary->saturated_stars = 0u;
    summary->edge_stars = 0u;
    summary->blended_stars = 0u;
    summary->low_snr_stars = 0u;
    for (i = 0u; i < star_count; ++i) {
        if ((stars[i].flags & IA_STAR_SATURATED) != 0u) ++summary->saturated_stars;
        if ((stars[i].flags & IA_STAR_EDGE) != 0u) ++summary->edge_stars;
        if ((stars[i].flags & IA_STAR_BLENDED) != 0u) ++summary->blended_stars;
        if ((stars[i].flags & IA_STAR_LOW_SNR) != 0u) ++summary->low_snr_stars;
        if ((stars[i].flags & reject_flags) == 0u) ++included;
    }
    summary->included_stars = included;
    if (included == 0u) {
        summary->median_hfd_milli_pixels = 0;
        summary->mad_hfd_milli_pixels = 0;
        summary->median_fwhm_milli_pixels = 0;
        summary->mad_fwhm_milli_pixels = 0;
        summary->median_ellipticity_milli = 0;
        summary->median_snr_milli = 0;
        return IA_NO_STARS;
    }

    values = (int32_t *)workspace;
#define IA_SUMMARY_FIELD(field, median_target, mad_target) \
    do { \
        uint32_t j = 0u; \
        for (i = 0u; i < star_count; ++i) { \
            if ((stars[i].flags & reject_flags) == 0u) values[j++] = stars[i].field; \
        } \
        median_target = ia_median_i32(values, included); \
        if ((mad_target) != NULL) { \
            for (j = 0u, i = 0u; i < star_count; ++i) { \
                if ((stars[i].flags & reject_flags) == 0u) values[j++] = stars[i].field; \
            } \
            *(int32_t *)(mad_target) = ia_mad_i32(values, included, median_target); \
        } \
    } while (0)

    IA_SUMMARY_FIELD(hfd_milli_pixels,
        summary->median_hfd_milli_pixels, &summary->mad_hfd_milli_pixels);
    IA_SUMMARY_FIELD(fwhm_milli_pixels,
        summary->median_fwhm_milli_pixels, &summary->mad_fwhm_milli_pixels);
    IA_SUMMARY_FIELD(ellipticity_milli,
        summary->median_ellipticity_milli, NULL);
    IA_SUMMARY_FIELD(snr_milli,
        summary->median_snr_milli, NULL);
#undef IA_SUMMARY_FIELD

    return IA_OK;
}

int IA_CALL IA_AnalyzeFrame(
    const uint16_t *pixels,
    uint32_t width,
    uint32_t height,
    uint32_t stride_pixels,
    const IAConfig *config,
    void *workspace,
    uint32_t workspace_bytes,
    IAStar *stars,
    uint32_t star_capacity,
    IAFrameSummary *summary)
{
    IAWorkspaceView view;
    uint32_t median;
    double noise;
    double threshold_delta;
    uint32_t threshold;
    uint32_t candidate_count;
    uint32_t selected_count;
    uint32_t i;
    int status;

    ia_debug_log(
        "IA_AnalyzeFrame enter width=%u height=%u stride=%u capacity=%u workspace=%u",
        width, height, stride_pixels, star_capacity, workspace_bytes);
    if (pixels == NULL || config == NULL || stars == NULL || summary == NULL ||
        width == 0u || height == 0u || stride_pixels < width || star_capacity == 0u ||
        config->struct_size < sizeof(*config) || summary->struct_size < sizeof(*summary) ||
        config->aperture_radius_pixels == 0u ||
        config->annulus_inner_radius_pixels <= config->aperture_radius_pixels ||
        config->annulus_outer_radius_pixels <= config->annulus_inner_radius_pixels) {
        ia_debug_log("IA_AnalyzeFrame invalid argument");
        return IA_INVALID_ARGUMENT;
    }

    status = ia_workspace_view(workspace, workspace_bytes, star_capacity, &view);
    if (status != IA_OK) {
        ia_debug_log("IA_AnalyzeFrame workspace status=%d", status);
        return status;
    }

    ia_global_background(pixels, width, height, stride_pixels, view.histogram, &median, &noise);
    threshold_delta = noise * config->detection_sigma_milli / IA_MILLI;
    if (threshold_delta < config->minimum_threshold_adu) {
        threshold_delta = config->minimum_threshold_adu;
    }
    threshold = median + ia_clamp_u32(threshold_delta);
    if (threshold > 65535u) {
        threshold = 65535u;
    }

    candidate_count = ia_detect_candidates(
        pixels, width, height, stride_pixels, threshold, config,
        view.candidates, view.candidate_capacity);
    selected_count = ia_select_separated(
        view.candidates, candidate_count, config->minimum_separation_pixels,
        view.selected, star_capacity);

    memset(summary, 0, sizeof(*summary));
    summary->struct_size = sizeof(*summary);
    summary->detected_stars = candidate_count;
    summary->measured_stars = selected_count;
    summary->background_milli_adu = ia_clamp_i32((double)median * IA_MILLI);
    summary->noise_milli_adu = ia_clamp_i32(noise * IA_MILLI);
    summary->detection_threshold_adu = threshold;

    for (i = 0u; i < selected_count; ++i) {
        ia_measure_candidate(
            pixels, width, height, stride_pixels, config, &view.selected[i],
            view.candidates, candidate_count, &stars[i]);
    }

    status = IA_SummarizeStars(
        stars, selected_count, 0u, view.histogram,
        IA_HISTOGRAM_BINS * (uint32_t)sizeof(uint32_t), summary);
    summary->detected_stars = candidate_count;
    summary->background_milli_adu = ia_clamp_i32((double)median * IA_MILLI);
    summary->noise_milli_adu = ia_clamp_i32(noise * IA_MILLI);
    summary->detection_threshold_adu = threshold;
    ia_debug_log(
        "IA_AnalyzeFrame exit status=%d detected=%u measured=%u included=%u threshold=%u",
        selected_count == 0u ? IA_NO_STARS : status,
        candidate_count, selected_count, summary->included_stars, threshold);
    return selected_count == 0u ? IA_NO_STARS : status;
}

int IA_CALL IA_MeasureLockedStars(
    const uint16_t *pixels,
    uint32_t width,
    uint32_t height,
    uint32_t stride_pixels,
    const IAConfig *config,
    const IAStar *locked_stars,
    uint32_t star_count,
    uint32_t search_radius_pixels,
    IAStar *measurements)
{
    uint32_t i;

    ia_debug_log(
        "IA_MeasureLockedStars enter width=%u height=%u stride=%u stars=%u radius=%u",
        width, height, stride_pixels, star_count, search_radius_pixels);
    if (pixels == NULL || config == NULL || locked_stars == NULL ||
        measurements == NULL || width == 0u || height == 0u ||
        stride_pixels < width || config->struct_size < sizeof(*config)) {
        ia_debug_log("IA_MeasureLockedStars invalid argument");
        return IA_INVALID_ARGUMENT;
    }

    for (i = 0u; i < star_count; ++i) {
        int32_t expected_x = locked_stars[i].x_milli / IA_MILLI;
        int32_t expected_y = locked_stars[i].y_milli / IA_MILLI;
        int32_t left = expected_x - (int32_t)search_radius_pixels;
        int32_t right = expected_x + (int32_t)search_radius_pixels;
        int32_t top = expected_y - (int32_t)search_radius_pixels;
        int32_t bottom = expected_y + (int32_t)search_radius_pixels;
        IACandidate candidate;
        int32_t y;
        int32_t x;

        if (left < 0) left = 0;
        if (top < 0) top = 0;
        if (right >= (int32_t)width) right = (int32_t)width - 1;
        if (bottom >= (int32_t)height) bottom = (int32_t)height - 1;
        ia_debug_log(
            "IA_MeasureLockedStars star=%u expected=(%d,%d) window=(%d,%d)-(%d,%d)",
            i, expected_x, expected_y, left, top, right, bottom);
        candidate.x = (uint32_t)(expected_x < 0 ? 0 : expected_x);
        candidate.y = (uint32_t)(expected_y < 0 ? 0 : expected_y);
        candidate.peak = 0u;
        candidate.area = 0u;
        for (y = top; y <= bottom; ++y) {
            const uint16_t *row = pixels + (size_t)y * stride_pixels;
            for (x = left; x <= right; ++x) {
                if (row[x] > candidate.peak) {
                    candidate.peak = row[x];
                    candidate.x = (uint32_t)x;
                    candidate.y = (uint32_t)y;
                }
            }
        }
        ia_measure_candidate(
            pixels, width, height, stride_pixels, config, &candidate,
            &candidate, 1u, &measurements[i]);
        ia_debug_log(
            "IA_MeasureLockedStars star=%u measured=(%d,%d) peak=%u hfd=%d flags=%u",
            i, measurements[i].x_milli, measurements[i].y_milli,
            measurements[i].peak_adu, measurements[i].hfd_milli_pixels,
            measurements[i].flags);
    }
    ia_debug_log("IA_MeasureLockedStars exit status=%d",
        star_count == 0u ? IA_NO_STARS : IA_OK);
    return star_count == 0u ? IA_NO_STARS : IA_OK;
}

int IA_CALL IA_RecommendExposure(
    uint32_t current_milliseconds,
    uint32_t measured_peak_adu,
    uint32_t background_adu,
    uint32_t target_peak_adu,
    uint32_t minimum_milliseconds,
    uint32_t maximum_milliseconds,
    uint32_t maximum_change_milli,
    uint32_t *recommended_milliseconds)
{
    double signal;
    double target_signal;
    double ratio;
    double maximum_ratio;
    double recommended;

    if (recommended_milliseconds == NULL || current_milliseconds == 0u ||
        maximum_milliseconds < minimum_milliseconds ||
        measured_peak_adu <= background_adu || target_peak_adu <= background_adu ||
        maximum_change_milli < IA_MILLI) {
        return IA_INVALID_ARGUMENT;
    }

    signal = measured_peak_adu - background_adu;
    target_signal = target_peak_adu - background_adu;
    ratio = target_signal / signal;
    maximum_ratio = (double)maximum_change_milli / IA_MILLI;
    if (ratio > maximum_ratio) ratio = maximum_ratio;
    if (ratio < 1.0 / maximum_ratio) ratio = 1.0 / maximum_ratio;
    recommended = current_milliseconds * ratio;
    if (recommended < minimum_milliseconds) recommended = minimum_milliseconds;
    if (recommended > maximum_milliseconds) recommended = maximum_milliseconds;
    *recommended_milliseconds = ia_clamp_u32(recommended);
    return IA_OK;
}

static uint32_t ia_exposure_from_range(
    uint32_t current_milliseconds,
    double required_minimum,
    double allowed_maximum,
    uint32_t minimum_milliseconds,
    uint32_t maximum_milliseconds,
    double maximum_ratio)
{
    double recommended;

    if (required_minimum <= allowed_maximum) {
        recommended = current_milliseconds;
        if (recommended < required_minimum) recommended = required_minimum;
        if (recommended > allowed_maximum) recommended = allowed_maximum;
    } else {
        recommended = allowed_maximum;
    }
    if (recommended > current_milliseconds * maximum_ratio) {
        recommended = current_milliseconds * maximum_ratio;
    }
    if (recommended < current_milliseconds / maximum_ratio) {
        recommended = current_milliseconds / maximum_ratio;
    }
    if (recommended < minimum_milliseconds) recommended = minimum_milliseconds;
    if (recommended > maximum_milliseconds) recommended = maximum_milliseconds;
    return ia_clamp_u32(recommended);
}

int IA_CALL IA_ComputeStarROI(
    const IAStar *stars,
    uint32_t star_count,
    uint32_t image_width,
    uint32_t image_height,
    uint32_t margin_pixels,
    IARect *roi)
{
    int64_t left;
    int64_t top;
    int64_t right;
    int64_t bottom;
    uint32_t i;

    ia_debug_log(
        "IA_ComputeStarROI enter stars=%u image=%ux%u margin=%u",
        star_count, image_width, image_height, margin_pixels);
    if (stars == NULL || star_count == 0u || image_width == 0u || image_height == 0u ||
        roi == NULL) {
        ia_debug_log("IA_ComputeStarROI invalid argument");
        return IA_INVALID_ARGUMENT;
    }

    left = stars[0].x_milli / IA_MILLI;
    right = left;
    top = stars[0].y_milli / IA_MILLI;
    bottom = top;
    for (i = 1u; i < star_count; ++i) {
        int64_t x = stars[i].x_milli / IA_MILLI;
        int64_t y = stars[i].y_milli / IA_MILLI;
        if (x < left) left = x;
        if (x > right) right = x;
        if (y < top) top = y;
        if (y > bottom) bottom = y;
    }

    left -= margin_pixels;
    top -= margin_pixels;
    right += margin_pixels;
    bottom += margin_pixels;
    if (left < 0) left = 0;
    if (top < 0) top = 0;
    if (right >= image_width) right = image_width - 1u;
    if (bottom >= image_height) bottom = image_height - 1u;

    roi->x = (int32_t)left;
    roi->y = (int32_t)top;
    roi->width = (int32_t)(right - left + 1);
    roi->height = (int32_t)(bottom - top + 1);
    ia_debug_log("IA_ComputeStarROI exit roi=(%d,%d) %dx%d",
        roi->x, roi->y, roi->width, roi->height);
    return IA_OK;
}

#define IA_FOCUS_MAGIC 0x49414632u
#define IA_FOCUS_CANDIDATES 64u

typedef struct IAFocusPoint {
    int32_t focus_position;
    int32_t hfd_milli_pixels;
    uint32_t flags;
} IAFocusPoint;

typedef struct IALine {
    double intercept;
    double slope;
} IALine;

typedef struct IAVFit {
    int32_t focus_milli_steps;
    int32_t rms_residual_milli;
} IAVFit;

typedef struct IAFocusRun {
    uint32_t magic;
    IAFocusConfig config;
    IAConfig analysis_config;
    IAStar candidates[IA_FOCUS_CANDIDATES];
    IAStar locked[IA_FOCUS_MAX_STARS];
    IAStar measurements[IA_FOCUS_MAX_STARS];
    IAStar previous[IA_FOCUS_MAX_STARS];
    IAFocusPoint points[IA_FOCUS_MAX_STARS][IA_FOCUS_MAX_POINTS];
    IARect roi;
    IAFocusResult result;
    int32_t center_focus;
    int32_t first_focus;
    uint32_t star_count;
    uint32_t point_index;
    uint32_t retry_count;
    uint32_t previous_valid;
    uint32_t last_exposure_ms;
} IAFocusRun;

static int ia_focus_config_valid(const IAFocusConfig *config)
{
    return config != NULL &&
        config->struct_size >= sizeof(*config) &&
        config->star_limit > 0u &&
        config->star_limit <= IA_FOCUS_MAX_STARS &&
        config->point_count >= 2u * config->minimum_branch_points &&
        config->point_count <= IA_FOCUS_MAX_POINTS &&
        (config->point_count & 1u) != 0u &&
        config->focus_increment > 0u &&
        config->target_minimum_peak_adu <
            config->target_maximum_peak_adu &&
        config->minimum_exposure_ms > 0u &&
        config->minimum_exposure_ms <= config->maximum_exposure_ms &&
        config->initial_maximum_change_milli >= IA_MILLI &&
        config->step_maximum_change_milli >= IA_MILLI &&
        config->minimum_branch_points >= 2u;
}

static void ia_focus_request(
    IAFocusRun *run,
    IAFocusRequest *request,
    uint32_t action,
    uint32_t reason_flags,
    int32_t focus_position,
    uint32_t exposure_milliseconds)
{
    memset(request, 0, sizeof(*request));
    request->struct_size = sizeof(*request);
    request->action = action;
    request->reason_flags = reason_flags;
    request->focus_position = focus_position;
    request->exposure_milliseconds = exposure_milliseconds;
    request->roi = run->roi;
    ia_debug_log(
        "focus request action=%u reason=%u focus=%d exposure=%u roi=(%d,%d %dx%d)",
        action, reason_flags, focus_position, exposure_milliseconds,
        request->roi.x, request->roi.y,
        request->roi.width, request->roi.height);
}

static int ia_focus_exposure(
    const IAFocusRun *run,
    const IAStar *measurements,
    const IAStar *previous,
    uint32_t current_milliseconds,
    uint32_t maximum_change_milli,
    int improving,
    uint32_t *recommended_milliseconds)
{
    double required_minimum = 0.0;
    double allowed_maximum = HUGE_VAL;
    uint32_t usable = 0u;
    uint32_t i;

    for (i = 0u; i < run->star_count; ++i) {
        double background =
            measurements[i].background_milli_adu / (double)IA_MILLI;
        double current_hfd = measurements[i].hfd_milli_pixels;
        double predicted_hfd = current_hfd;
        double measured_signal;
        double target_minimum_signal;
        double target_maximum_signal;
        double exposure_scale;
        double star_minimum;
        double star_maximum;

        if (background < 0.0) background = 0.0;
        if (current_hfd <= 0.0 || measurements[i].peak_adu <= background) {
            continue;
        }
        if (previous != NULL && previous[i].hfd_milli_pixels > 0) {
            double change = fabs(
                current_hfd - previous[i].hfd_milli_pixels);
            predicted_hfd = improving ?
                fmax(current_hfd * 0.5, current_hfd - change) :
                fmin(current_hfd * 2.0, current_hfd + change);
        }
        measured_signal = measurements[i].peak_adu - background;
        target_minimum_signal =
            run->config.target_minimum_peak_adu - background;
        target_maximum_signal =
            run->config.target_maximum_peak_adu - background;
        if (target_minimum_signal <= 0.0 ||
            target_maximum_signal <= target_minimum_signal) {
            continue;
        }
        exposure_scale =
            current_milliseconds /
            (measured_signal * current_hfd * current_hfd);
        star_minimum =
            exposure_scale * target_minimum_signal *
            predicted_hfd * predicted_hfd;
        star_maximum =
            exposure_scale * target_maximum_signal *
            predicted_hfd * predicted_hfd;
        if (star_minimum > required_minimum) required_minimum = star_minimum;
        if (star_maximum < allowed_maximum) allowed_maximum = star_maximum;
        ++usable;
    }
    if (usable == 0u) {
        return IA_NO_STARS;
    }
    *recommended_milliseconds = ia_exposure_from_range(
        current_milliseconds, required_minimum, allowed_maximum,
        run->config.minimum_exposure_ms,
        run->config.maximum_exposure_ms,
        (double)maximum_change_milli / IA_MILLI);
    return IA_OK;
}

static int ia_lock_compact_stars(IAFocusRun *run, uint32_t measured_stars)
{
    uint8_t chosen[IA_FOCUS_CANDIDATES] = { 0 };
    uint32_t desired;
    uint32_t first = UINT32_MAX;
    int32_t left;
    int32_t right;
    int32_t top;
    int32_t bottom;
    uint32_t i;

    for (i = 0u; i < measured_stars; ++i) {
        if ((run->candidates[i].flags &
             run->config.reject_star_flags) == 0u) {
            first = i;
            break;
        }
    }
    if (first == UINT32_MAX) {
        return IA_NO_STARS;
    }
    desired = run->config.star_limit;
    run->locked[0] = run->candidates[first];
    chosen[first] = 1u;
    run->star_count = 1u;
    left = right = run->locked[0].x_milli;
    top = bottom = run->locked[0].y_milli;

    while (run->star_count < desired) {
        uint32_t best = UINT32_MAX;
        uint64_t best_area = UINT64_MAX;
        for (i = 0u; i < measured_stars; ++i) {
            int32_t candidate_left;
            int32_t candidate_right;
            int32_t candidate_top;
            int32_t candidate_bottom;
            uint64_t area;
            if (chosen[i] ||
                (run->candidates[i].flags &
                 run->config.reject_star_flags) != 0u) {
                continue;
            }
            candidate_left =
                run->candidates[i].x_milli < left ?
                run->candidates[i].x_milli : left;
            candidate_right =
                run->candidates[i].x_milli > right ?
                run->candidates[i].x_milli : right;
            candidate_top =
                run->candidates[i].y_milli < top ?
                run->candidates[i].y_milli : top;
            candidate_bottom =
                run->candidates[i].y_milli > bottom ?
                run->candidates[i].y_milli : bottom;
            area =
                (uint64_t)(candidate_right - candidate_left + 1) *
                (uint64_t)(candidate_bottom - candidate_top + 1);
            if (area < best_area) {
                best_area = area;
                best = i;
            }
        }
        if (best == UINT32_MAX) break;
        run->locked[run->star_count++] = run->candidates[best];
        chosen[best] = 1u;
        if (run->candidates[best].x_milli < left) {
            left = run->candidates[best].x_milli;
        }
        if (run->candidates[best].x_milli > right) {
            right = run->candidates[best].x_milli;
        }
        if (run->candidates[best].y_milli < top) {
            top = run->candidates[best].y_milli;
        }
        if (run->candidates[best].y_milli > bottom) {
            bottom = run->candidates[best].y_milli;
        }
    }
    return IA_OK;
}

static int ia_measure_focus_frame(
    IAFocusRun *run,
    const uint16_t *pixels,
    uint32_t width,
    uint32_t height,
    uint32_t stride_pixels,
    int32_t origin_x,
    int32_t origin_y)
{
    IAStar local_locked[IA_FOCUS_MAX_STARS];
    uint32_t i;
    int status;

    for (i = 0u; i < run->star_count; ++i) {
        local_locked[i] = run->locked[i];
        local_locked[i].x_milli -= origin_x * IA_MILLI;
        local_locked[i].y_milli -= origin_y * IA_MILLI;
    }
    status = IA_MeasureLockedStars(
        pixels, width, height, stride_pixels,
        &run->analysis_config, local_locked, run->star_count,
        run->config.search_radius_pixels, run->measurements);
    if (status != IA_OK) return status;
    for (i = 0u; i < run->star_count; ++i) {
        run->locked[i] = run->measurements[i];
        run->locked[i].x_milli += origin_x * IA_MILLI;
        run->locked[i].y_milli += origin_y * IA_MILLI;
    }
    return IA_OK;
}

static int ia_regression(
    const IAFocusPoint *points,
    uint32_t first,
    uint32_t count,
    IALine *line)
{
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_xx = 0.0;
    double sum_xy = 0.0;
    double denominator;
    uint32_t i;

    for (i = first; i < first + count; ++i) {
        double x = points[i].focus_position;
        double y = points[i].hfd_milli_pixels;
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
    }
    denominator = count * sum_xx - sum_x * sum_x;
    if (denominator == 0.0) return IA_INVALID_ARGUMENT;
    line->slope =
        (count * sum_xy - sum_x * sum_y) / denominator;
    line->intercept = (sum_y - line->slope * sum_x) / count;
    return IA_OK;
}

static int ia_fit_v(
    const IAFocusRun *run,
    uint32_t star_index,
    IAVFit *fit)
{
    IAFocusPoint valid[IA_FOCUS_MAX_POINTS];
    double best_error = HUGE_VAL;
    double best_focus = 0.0;
    uint32_t count = 0u;
    uint32_t split;
    uint32_t i;

    for (i = 0u; i < run->config.point_count; ++i) {
        const IAFocusPoint *point = &run->points[star_index][i];
        if ((point->flags & run->config.reject_star_flags) == 0u &&
            point->hfd_milli_pixels > 0) {
            valid[count++] = *point;
        }
    }
    if (count < 2u * run->config.minimum_branch_points) {
        return IA_NO_STARS;
    }

    for (split = run->config.minimum_branch_points;
         split + run->config.minimum_branch_points <= count;
         ++split) {
        IALine left;
        IALine right;
        double focus;
        double squared_error = 0.0;
        if (ia_regression(valid, 0u, split, &left) != IA_OK ||
            ia_regression(valid, split, count - split, &right) != IA_OK ||
            left.slope >= 0.0 || right.slope <= 0.0 ||
            left.slope == right.slope) {
            continue;
        }
        focus =
            (right.intercept - left.intercept) /
            (left.slope - right.slope);
        if (focus < valid[split - 1u].focus_position ||
            focus > valid[split].focus_position) {
            continue;
        }
        for (i = 0u; i < count; ++i) {
            const IALine *line = i < split ? &left : &right;
            double residual =
                valid[i].hfd_milli_pixels -
                (line->intercept +
                 line->slope * valid[i].focus_position);
            squared_error += residual * residual;
        }
        squared_error /= count;
        if (squared_error < best_error) {
            best_error = squared_error;
            best_focus = focus;
        }
    }
    if (best_error == HUGE_VAL) {
        return IA_NO_STARS;
    }
    fit->focus_milli_steps = ia_clamp_i32(best_focus * IA_MILLI);
    fit->rms_residual_milli = ia_clamp_i32(sqrt(best_error));
    return IA_OK;
}

static int ia_finish_focus_run(
    IAFocusRun *run,
    uint32_t recommended_exposure_ms)
{
    int32_t focus_values[IA_FOCUS_MAX_STARS];
    int32_t residual_values[IA_FOCUS_MAX_STARS];
    int32_t deviations[IA_FOCUS_MAX_STARS];
    uint32_t fitted = 0u;
    uint32_t i;

    for (i = 0u; i < run->star_count; ++i) {
        IAVFit fit;
        if (ia_fit_v(run, i, &fit) == IA_OK &&
            fit.rms_residual_milli <=
                (int32_t)run->config.maximum_residual_milli) {
            focus_values[fitted] = fit.focus_milli_steps;
            residual_values[fitted] = fit.rms_residual_milli;
            ++fitted;
        }
    }
    memset(&run->result, 0, sizeof(run->result));
    run->result.struct_size = sizeof(run->result);
    run->result.stars_locked = run->star_count;
    run->result.stars_fitted = fitted;
    run->result.recommended_exposure_milliseconds =
        recommended_exposure_ms;
    if (fitted == 0u) return IA_NO_STARS;

    run->result.focus_milli_steps =
        ia_median_i32(focus_values, fitted);
    for (i = 0u; i < fitted; ++i) {
        int64_t difference =
            (int64_t)focus_values[i] -
            run->result.focus_milli_steps;
        if (difference < 0) difference = -difference;
        deviations[i] =
            difference > INT32_MAX ? INT32_MAX : (int32_t)difference;
    }
    run->result.mad_focus_milli_steps =
        ia_median_i32(deviations, fitted);
    run->result.median_rms_residual_milli =
        ia_median_i32(residual_values, fitted);
    ia_debug_log(
        "focus complete locked=%u fitted=%u focus=%d mad=%d residual=%d verify_exposure=%u",
        run->result.stars_locked, run->result.stars_fitted,
        run->result.focus_milli_steps,
        run->result.mad_focus_milli_steps,
        run->result.median_rms_residual_milli,
        run->result.recommended_exposure_milliseconds);
    return IA_OK;
}

uint32_t IA_CALL IA_FocusRunBytes(void)
{
    return (uint32_t)sizeof(IAFocusRun) +
        ia_workspace_bytes(IA_FOCUS_CANDIDATES);
}

void IA_CALL IA_DefaultFocusConfig(IAFocusConfig *config)
{
    if (config == NULL) return;
    memset(config, 0, sizeof(*config));
    config->struct_size = sizeof(*config);
    config->star_limit = IA_FOCUS_MAX_STARS;
    config->point_count = 11u;
    config->focus_increment = 20u;
    config->target_minimum_peak_adu = 20000u;
    config->target_maximum_peak_adu = 50000u;
    config->minimum_exposure_ms = 100u;
    config->maximum_exposure_ms = 10000u;
    config->initial_maximum_change_milli = 4000u;
    config->step_maximum_change_milli = 1500u;
    config->maximum_retries = 2u;
    config->roi_margin_pixels = 64u;
    config->search_radius_pixels = 16u;
    config->reject_star_flags =
        IA_STAR_SATURATED | IA_STAR_EDGE | IA_STAR_BLENDED |
        IA_STAR_TOO_SMALL | IA_STAR_LOW_SNR;
    config->minimum_branch_points = 3u;
    config->maximum_residual_milli = 1000u;
}

int IA_CALL IA_BeginFocusRun(
    void *run_memory,
    uint32_t run_bytes,
    const IAFocusConfig *focus_config,
    const IAConfig *analysis_config,
    const uint16_t *pixels,
    uint32_t width,
    uint32_t height,
    uint32_t stride_pixels,
    int32_t center_focus_position,
    uint32_t exposure_milliseconds,
    IAFocusRequest *request)
{
    IAFocusRun *run = (IAFocusRun *)run_memory;
    IAFrameSummary summary;
    void *workspace;
    uint32_t recommended;
    int status;

    if (run == NULL || run_bytes < IA_FocusRunBytes() ||
        !ia_focus_config_valid(focus_config) ||
        analysis_config == NULL ||
        analysis_config->struct_size < sizeof(*analysis_config) ||
        pixels == NULL || width == 0u || height == 0u ||
        stride_pixels < width || exposure_milliseconds == 0u ||
        request == NULL || request->struct_size < sizeof(*request)) {
        return IA_INVALID_ARGUMENT;
    }
    memset(run, 0, sizeof(*run));
    run->magic = IA_FOCUS_MAGIC;
    run->config = *focus_config;
    run->analysis_config = *analysis_config;
    run->center_focus = center_focus_position;
    run->first_focus =
        center_focus_position -
        (int32_t)(focus_config->point_count / 2u *
                  focus_config->focus_increment);
    workspace = (uint8_t *)run_memory + sizeof(*run);
    summary.struct_size = sizeof(summary);
    status = IA_AnalyzeFrame(
        pixels, width, height, stride_pixels,
        &run->analysis_config, workspace,
        ia_workspace_bytes(IA_FOCUS_CANDIDATES),
        run->candidates, IA_FOCUS_CANDIDATES, &summary);
    if (status != IA_OK) return status;
    status = ia_lock_compact_stars(run, summary.measured_stars);
    if (status != IA_OK) return status;
    status = IA_ComputeStarROI(
        run->locked, run->star_count, width, height,
        run->config.roi_margin_pixels, &run->roi);
    if (status != IA_OK) return status;
    status = ia_focus_exposure(
        run, run->locked, NULL, exposure_milliseconds,
        run->config.initial_maximum_change_milli, 1, &recommended);
    if (status != IA_OK) return status;
    run->last_exposure_ms = recommended;
    ia_debug_log(
        "focus begin center=%d first=%d stars=%u points=%u",
        run->center_focus, run->first_focus,
        run->star_count, run->config.point_count);
    ia_focus_request(
        run, request, IA_FOCUS_CAPTURE, IA_FOCUS_REASON_NONE,
        run->first_focus, recommended);
    return IA_OK;
}

int IA_CALL IA_SubmitFocusFrame(
    void *run_memory,
    const uint16_t *pixels,
    uint32_t width,
    uint32_t height,
    uint32_t stride_pixels,
    int32_t origin_x,
    int32_t origin_y,
    int32_t focus_position,
    uint32_t exposure_milliseconds,
    IAFocusRequest *request)
{
    IAFocusRun *run = (IAFocusRun *)run_memory;
    const IAStar *previous;
    uint32_t reason = IA_FOCUS_REASON_NONE;
    uint32_t usable = 0u;
    uint32_t recommended = exposure_milliseconds;
    uint32_t i;
    int status;

    if (run == NULL || run->magic != IA_FOCUS_MAGIC ||
        pixels == NULL || width == 0u || height == 0u ||
        stride_pixels < width || exposure_milliseconds == 0u ||
        request == NULL || request->struct_size < sizeof(*request)) {
        return IA_INVALID_ARGUMENT;
    }
    status = ia_measure_focus_frame(
        run, pixels, width, height, stride_pixels, origin_x, origin_y);
    if (status != IA_OK) return status;
    for (i = 0u; i < run->star_count; ++i) {
        if ((run->measurements[i].flags &
             run->config.reject_star_flags) == 0u) {
            ++usable;
        }
        if ((run->measurements[i].flags & IA_STAR_SATURATED) != 0u) {
            reason |= IA_FOCUS_REASON_SATURATED;
        }
        if ((run->measurements[i].flags & IA_STAR_LOW_SNR) != 0u) {
            reason |= IA_FOCUS_REASON_TOO_FAINT;
        }
    }
    previous = run->previous_valid ? run->previous : NULL;
    if (usable == 0u) {
        reason |= IA_FOCUS_REASON_NO_STARS;
        if (run->retry_count < run->config.maximum_retries &&
            ia_focus_exposure(
                run, run->measurements, NULL, exposure_milliseconds,
                run->config.step_maximum_change_milli, 1,
                &recommended) == IA_OK) {
            ++run->retry_count;
            ia_focus_request(
                run, request, IA_FOCUS_RETAKE, reason,
                focus_position, recommended);
            return IA_OK;
        }
        ia_focus_request(
            run, request, IA_FOCUS_FAILED, reason,
            focus_position, exposure_milliseconds);
        return IA_OK;
    }

    for (i = 0u; i < run->star_count; ++i) {
        IAFocusPoint *point =
            &run->points[i][run->point_index];
        point->focus_position = focus_position;
        point->hfd_milli_pixels =
            run->measurements[i].hfd_milli_pixels;
        point->flags = run->measurements[i].flags;
    }
    run->retry_count = 0u;
    if (run->point_index + 1u >= run->config.point_count) {
        if (ia_focus_exposure(
            run, run->measurements, previous, exposure_milliseconds,
            run->config.step_maximum_change_milli, 1,
            &recommended) != IA_OK) {
            recommended = exposure_milliseconds;
        }
        if (ia_finish_focus_run(run, recommended) == IA_OK) {
            ia_focus_request(
                run, request, IA_FOCUS_COMPLETE, reason,
                (run->result.focus_milli_steps + 500) / IA_MILLI,
                run->result.recommended_exposure_milliseconds);
        } else {
            ia_focus_request(
                run, request, IA_FOCUS_FAILED,
                reason | IA_FOCUS_REASON_NO_CURVE,
                focus_position, exposure_milliseconds);
        }
        return IA_OK;
    }

    {
        int32_t next_focus =
            run->first_focus +
            (int32_t)((run->point_index + 1u) *
                      run->config.focus_increment);
        int improving =
            abs(next_focus - run->center_focus) <
            abs(focus_position - run->center_focus);
        if (ia_focus_exposure(
            run, run->measurements, previous, exposure_milliseconds,
            run->config.step_maximum_change_milli, improving,
            &recommended) != IA_OK) {
            recommended = exposure_milliseconds;
        }
        memcpy(
            run->previous, run->measurements,
            run->star_count * sizeof(IAStar));
        run->previous_valid = 1u;
        ++run->point_index;
        run->last_exposure_ms = recommended;
        ia_focus_request(
            run, request, IA_FOCUS_CAPTURE, reason,
            next_focus, recommended);
    }
    return IA_OK;
}

int IA_CALL IA_EndFocusRun(
    void *run_memory,
    IAFocusResult *result)
{
    IAFocusRun *run = (IAFocusRun *)run_memory;
    if (run == NULL || run->magic != IA_FOCUS_MAGIC ||
        result == NULL || result->struct_size < sizeof(*result) ||
        run->result.stars_fitted == 0u) {
        return IA_INVALID_ARGUMENT;
    }
    *result = run->result;
    return IA_OK;
}

int IA_CALL IA_DiagnoseFocusFrame(
    void *run_memory,
    const uint16_t *pixels,
    uint32_t width,
    uint32_t height,
    uint32_t stride_pixels,
    int32_t origin_x,
    int32_t origin_y,
    IAFocusDiagnostic *diagnostic)
{
    IAFocusRun *run = (IAFocusRun *)run_memory;
    int32_t values[IA_FOCUS_MAX_STARS];
    int32_t deviations[IA_FOCUS_MAX_STARS];
    uint32_t measured = 0u;
    uint32_t saturated = 0u;
    uint32_t i;
    int status;

    if (run == NULL || run->magic != IA_FOCUS_MAGIC ||
        diagnostic == NULL ||
        diagnostic->struct_size < sizeof(*diagnostic)) {
        return IA_INVALID_ARGUMENT;
    }
    status = ia_measure_focus_frame(
        run, pixels, width, height, stride_pixels, origin_x, origin_y);
    if (status != IA_OK) return status;
    for (i = 0u; i < run->star_count; ++i) {
        if ((run->measurements[i].flags & IA_STAR_SATURATED) != 0u) {
            ++saturated;
        }
        if ((run->measurements[i].flags &
             run->config.reject_star_flags) == 0u &&
            run->measurements[i].hfd_milli_pixels > 0) {
            values[measured++] =
                run->measurements[i].hfd_milli_pixels;
        }
    }
    memset(diagnostic, 0, sizeof(*diagnostic));
    diagnostic->struct_size = sizeof(*diagnostic);
    diagnostic->stars_measured = measured;
    diagnostic->saturated_stars = saturated;
    if (measured == 0u) return IA_NO_STARS;
    diagnostic->median_hfd_milli_pixels =
        ia_median_i32(values, measured);
    for (i = 0u; i < measured; ++i) {
        int64_t difference =
            (int64_t)values[i] -
            diagnostic->median_hfd_milli_pixels;
        if (difference < 0) difference = -difference;
        deviations[i] =
            difference > INT32_MAX ? INT32_MAX : (int32_t)difference;
    }
    diagnostic->mad_hfd_milli_pixels =
        ia_median_i32(deviations, measured);
    return IA_OK;
}
