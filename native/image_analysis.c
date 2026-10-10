#include "image_analysis.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

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
        if (distance_squared <= (uint64_t)radius * radius) {
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
    return IA_API_VERSION;
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

    if (pixels == NULL || config == NULL || stars == NULL || summary == NULL ||
        width == 0u || height == 0u || stride_pixels < width || star_capacity == 0u ||
        config->struct_size < sizeof(*config) || summary->struct_size < sizeof(*summary) ||
        config->aperture_radius_pixels == 0u ||
        config->annulus_inner_radius_pixels <= config->aperture_radius_pixels ||
        config->annulus_outer_radius_pixels <= config->annulus_inner_radius_pixels) {
        return IA_INVALID_ARGUMENT;
    }

    status = ia_workspace_view(workspace, workspace_bytes, star_capacity, &view);
    if (status != IA_OK) {
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

    if (pixels == NULL || config == NULL || locked_stars == NULL ||
        measurements == NULL || width == 0u || height == 0u ||
        stride_pixels < width || config->struct_size < sizeof(*config)) {
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
    }
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

    if (stars == NULL || star_count == 0u || image_width == 0u || image_height == 0u ||
        roi == NULL) {
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
    return IA_OK;
}

static double ia_focus_shape(const IAFocusModel *model, double delta_steps)
{
    if (model->kind == IA_FOCUS_PIECEWISE_V) {
        double slope_micro = delta_steps < 0.0 ?
            model->left_slope_micro_per_step :
            model->right_slope_micro_per_step;
        return model->baseline_milli +
            fabs(delta_steps) * slope_micro / IA_MILLI;
    }
    if (model->kind == IA_FOCUS_HYPERBOLA) {
        double slope_metric = delta_steps *
            model->hyperbola_slope_micro_per_step / IA_MILLI;
        double radius = model->hyperbola_radius_milli;
        return model->baseline_milli +
            sqrt(radius * radius + slope_metric * slope_metric);
    }
    return 0.0;
}

static int ia_valid_focus_model(const IAFocusModel *model)
{
    if (model == NULL || model->struct_size < sizeof(*model)) {
        return 0;
    }
    if (model->kind == IA_FOCUS_PIECEWISE_V) {
        return model->left_slope_micro_per_step > 0 &&
            model->right_slope_micro_per_step > 0;
    }
    if (model->kind == IA_FOCUS_HYPERBOLA) {
        return model->hyperbola_radius_milli > 0 &&
            model->hyperbola_slope_micro_per_step > 0;
    }
    return 0;
}

int IA_CALL IA_FitFocusModel(
    const IAFocusModel *model,
    const IAFocusSample *samples,
    uint32_t sample_count,
    uint32_t reject_sample_flags,
    uint32_t resolution_milli_steps,
    IAFocusFit *fit)
{
    int32_t minimum_position = INT32_MAX;
    int32_t maximum_position = INT32_MIN;
    int64_t first_center;
    int64_t last_center;
    int64_t center;
    int64_t best_center = 0;
    double best_offset = 0.0;
    double best_error = HUGE_VAL;
    uint32_t used = 0u;
    uint32_t i;

    if (samples == NULL || fit == NULL || fit->struct_size < sizeof(*fit) ||
        sample_count == 0u || resolution_milli_steps == 0u) {
        return IA_INVALID_ARGUMENT;
    }
    memset(fit, 0, sizeof(*fit));
    fit->struct_size = sizeof(*fit);
    if (!ia_valid_focus_model(model)) {
        fit->flags = IA_FOCUS_INVALID_MODEL;
        return IA_INVALID_ARGUMENT;
    }

    for (i = 0u; i < sample_count; ++i) {
        if ((samples[i].flags & reject_sample_flags) == 0u &&
            samples[i].weight_milli > 0u) {
            if (samples[i].focus_position < minimum_position) {
                minimum_position = samples[i].focus_position;
            }
            if (samples[i].focus_position > maximum_position) {
                maximum_position = samples[i].focus_position;
            }
            ++used;
        }
    }
    fit->points_used = used;
    if (used < 3u || minimum_position >= maximum_position) {
        fit->flags = IA_FOCUS_INSUFFICIENT_POINTS;
        return IA_NO_STARS;
    }

    first_center = (int64_t)minimum_position * IA_MILLI;
    last_center = (int64_t)maximum_position * IA_MILLI;
    for (center = first_center; center <= last_center; center += resolution_milli_steps) {
        double weighted_offset_sum = 0.0;
        double weight_sum = 0.0;
        double squared_error = 0.0;
        double offset;

        for (i = 0u; i < sample_count; ++i) {
            if ((samples[i].flags & reject_sample_flags) == 0u &&
                samples[i].weight_milli > 0u) {
                double weight = (double)samples[i].weight_milli / IA_MILLI;
                double delta = samples[i].focus_position - (double)center / IA_MILLI;
                double shape = ia_focus_shape(model, delta);
                weighted_offset_sum += weight * (samples[i].metric_milli - shape);
                weight_sum += weight;
            }
        }
        offset = weighted_offset_sum / weight_sum;
        for (i = 0u; i < sample_count; ++i) {
            if ((samples[i].flags & reject_sample_flags) == 0u &&
                samples[i].weight_milli > 0u) {
                double weight = (double)samples[i].weight_milli / IA_MILLI;
                double delta = samples[i].focus_position - (double)center / IA_MILLI;
                double residual = samples[i].metric_milli -
                    (ia_focus_shape(model, delta) + offset);
                squared_error += weight * residual * residual;
            }
        }
        squared_error /= weight_sum;
        if (squared_error < best_error) {
            best_error = squared_error;
            best_center = center;
            best_offset = offset;
        }
        if (last_center - center < resolution_milli_steps) {
            break;
        }
    }

    fit->focus_milli_steps = best_center > INT32_MAX ? INT32_MAX : (int32_t)best_center;
    fit->vertical_offset_milli = ia_clamp_i32(best_offset);
    fit->rms_residual_milli = ia_clamp_i32(sqrt(best_error));
    for (i = 0u; i < sample_count; ++i) {
        if ((samples[i].flags & reject_sample_flags) == 0u &&
            samples[i].weight_milli > 0u) {
            int64_t sample_focus = (int64_t)samples[i].focus_position * IA_MILLI;
            if (sample_focus < best_center) ++fit->points_left;
            if (sample_focus > best_center) ++fit->points_right;
        }
    }
    if (fit->points_left == 0u || fit->points_right == 0u) {
        fit->flags |= IA_FOCUS_ONE_SIDED;
    }
    if (best_center == first_center ||
        best_center + resolution_milli_steps > last_center) {
        fit->flags |= IA_FOCUS_AT_BOUNDARY;
    }
    return IA_OK;
}

int IA_CALL IA_CombineFocusFits(
    const IAFocusFit *fits,
    uint32_t fit_count,
    uint32_t reject_fit_flags,
    void *workspace,
    uint32_t workspace_bytes,
    IAFocusSummary *summary)
{
    int32_t *values;
    uint32_t included = 0u;
    uint32_t i;

    if (fits == NULL || summary == NULL || summary->struct_size < sizeof(*summary)) {
        return IA_INVALID_ARGUMENT;
    }
    if (fit_count > 0u &&
        (workspace == NULL || workspace_bytes < fit_count * sizeof(int32_t))) {
        return IA_WORKSPACE_TOO_SMALL;
    }

    memset(summary, 0, sizeof(*summary));
    summary->struct_size = sizeof(*summary);
    values = (int32_t *)workspace;
    for (i = 0u; i < fit_count; ++i) {
        if ((fits[i].flags & reject_fit_flags) == 0u) {
            values[included++] = fits[i].focus_milli_steps;
        }
    }
    summary->fits_included = included;
    if (included == 0u) {
        return IA_NO_STARS;
    }
    summary->median_focus_milli_steps = ia_median_i32(values, included);
    for (i = 0u, included = 0u; i < fit_count; ++i) {
        if ((fits[i].flags & reject_fit_flags) == 0u) {
            values[included++] = fits[i].focus_milli_steps;
        }
    }
    summary->mad_focus_milli_steps = ia_mad_i32(
        values, included, summary->median_focus_milli_steps);
    for (i = 0u, included = 0u; i < fit_count; ++i) {
        if ((fits[i].flags & reject_fit_flags) == 0u) {
            values[included++] = fits[i].rms_residual_milli;
        }
    }
    summary->median_rms_residual_milli = ia_median_i32(values, included);
    return IA_OK;
}
