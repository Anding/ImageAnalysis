#include "image_analysis.h"
#include "analysis_parameters.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

enum { RADIAL_BINS = 2 * APERTURE_RADIUS_PIXELS * RADIAL_BINS_PER_PIXEL + 1 };

typedef struct Candidate {
    double flux;
    StarRecord star;
} Candidate;

/* Brightest first; positional tie-breaks are ascending, not discovery order. */
static int order(double af, double ax, double ay, double bf, double bx, double by)
{
    if (af != bf) return af > bf ? -1 : 1;
    if (ax != bx) return ax < bx ? -1 : 1;
    if (ay != by) return ay < by ? -1 : 1;
    return 0;
}

static int star_order(const void *left, const void *right)
{
    const StarRecord *a = left, *b = right;
    return order(a->flux, a->x, a->y, b->flux, b->x, b->y);
}

static int candidate_order(const Candidate *a, const Candidate *b)
{
    return order(a->flux, a->star.x, a->star.y, b->flux, b->star.x, b->star.y);
}

/* Risks describe this exposure only; they never decide detection or curve validity. */
static void quality(StarRecord *star, double sky, double mad)
{
    double spread = star->qxx + star->qyy;
    double elongation = hypot(star->qxx - star->qyy, 2 * star->qxy);
    double noise = fmax(MIN_NOISE_ADU, 1.4826 * mad);
    double significance = fmax(0, star->max_pixel - sky) / noise;
    double safe_peak = SAFE_PEAK_FRACTION * 65535;
    double faint_scale = significance > 0 ?
        pow(TARGET_PEAK_SIGMA / significance, 2) : 1;
    double bright_scale = star->max_pixel > 0 ?
        TARGET_PEAK_FRACTION * 65535 / star->max_pixel : 1;

    star->roundness = spread + elongation > 0 ?
        sqrt(fmax(0, (spread - elongation) / (spread + elongation))) : 0;
    star->over_exposed_risk = fmin(1, fmax(0,
        (star->max_pixel - safe_peak) / (65535 - safe_peak)));
    star->under_exposed_risk = fmin(1, fmax(0, 1 - significance / TARGET_PEAK_SIGMA));
    /* A calculated baseline balances sky significance with ADC headroom.
     * If requirements conflict, headroom wins; subsequent trials supply evidence.
     */
    star->exposure_scale = fmin(bright_scale, fmax(1, faint_scale));
}

/* Keep only the best candidates in a small, ordered, bounded array. */
static void retain(Candidate *candidates, int *count, int limit, Candidate candidate)
{
    int i = *count;
    if (i == limit) {
        Candidate last = candidates[i - 1];
        if (candidate_order(&candidate, &last) >= 0)
            return;
        --i;
    } else {
        ++*count;
    }
    while (i > 0) {
        Candidate previous = candidates[i - 1];
        if (candidate_order(&candidate, &previous) >= 0)
            break;
        candidates[i] = previous;
        --i;
    }
    candidates[i] = candidate;
}

static int inside(double x, double y, uint32_t width, uint32_t height)
{
    return x >= APERTURE_RADIUS_PIXELS && y >= APERTURE_RADIUS_PIXELS &&
        x + APERTURE_RADIUS_PIXELS <= width - 1.0 &&
        y + APERTURE_RADIUS_PIXELS <= height - 1.0;
}

/* Refine the detection centroid with all positive signal in its aperture. */
static int measure(const uint16_t *pixels, uint32_t width, uint32_t height,
    double sky, double cx, double cy, StarRecord *star)
{
    double flux = 0, sx = 0, sy = 0;
    double radial[RADIAL_BINS] = {0};
    double xx = 0, yy = 0, xy = 0, cumulative = 0;
    int left, right, top, bottom, x, y;
    unsigned bin;

    if (!inside(cx, cy, width, height))
        return 0;
    left = (int)ceil(cx - APERTURE_RADIUS_PIXELS);
    right = (int)floor(cx + APERTURE_RADIUS_PIXELS);
    top = (int)ceil(cy - APERTURE_RADIUS_PIXELS);
    bottom = (int)floor(cy + APERTURE_RADIUS_PIXELS);
    for (y = top; y <= bottom; ++y) {
        for (x = left; x <= right; ++x) {
            double dx = x - cx, dy = y - cy;
            double signal = pixels[(size_t)y * width + x] - sky;
            if (dx * dx + dy * dy <= APERTURE_RADIUS_PIXELS * APERTURE_RADIUS_PIXELS &&
                signal > 0) {
                flux += signal;
                sx += signal * x;
                sy += signal * y;
            }
        }
    }
    if (flux == 0)
        return 0;
    cx = sx / flux;
    cy = sy / flux;
    if (!inside(cx, cy, width, height))
        return 0;

    left = (int)ceil(cx - APERTURE_RADIUS_PIXELS);
    right = (int)floor(cx + APERTURE_RADIUS_PIXELS);
    top = (int)ceil(cy - APERTURE_RADIUS_PIXELS);
    bottom = (int)floor(cy + APERTURE_RADIUS_PIXELS);
    memset(star, 0, sizeof(*star));
    flux = sx = sy = 0;
    for (y = top; y <= bottom; ++y) {
        for (x = left; x <= right; ++x) {
            double dx = x - cx, dy = y - cy;
            double r2 = dx * dx + dy * dy;
            uint16_t value;
            double signal;
            if (r2 > APERTURE_RADIUS_PIXELS * APERTURE_RADIUS_PIXELS)
                continue;
            value = pixels[(size_t)y * width + x];
            if (value > star->max_pixel)
                star->max_pixel = value;
            if (value == 65535)
                ++star->saturated_pixel_count;
            signal = value - sky;
            if (signal <= 0)
                continue;
            flux += signal;
            sx += signal * dx;
            sy += signal * dy;
            xx += signal * dx * dx;
            yy += signal * dy * dy;
            xy += signal * dx * dy;
            ++star->positive_pixel_count;
        }
    }
    if (flux == 0)
        return 0;

    /* Central moments use the centroid of the final aperture's actual pixels. */
    star->x = cx + sx / flux;
    star->y = cy + sy / flux;
    if (!inside(star->x, star->y, width, height)) {
        memset(star, 0, sizeof(*star));
        return 0;
    }
    star->flux = flux;
    star->qxx = fmax(0, xx / flux - (sx / flux) * (sx / flux));
    star->qyy = fmax(0, yy / flux - (sy / flux) * (sy / flux));
    star->qxy = xy / flux - sx * sy / (flux * flux);

    /* Same aperture, radial distances about the returned subpixel centroid. */
    for (y = top; y <= bottom; ++y) {
        for (x = left; x <= right; ++x) {
            double dx = x - cx, dy = y - cy;
            double signal = pixels[(size_t)y * width + x] - sky;
            if (dx * dx + dy * dy <= APERTURE_RADIUS_PIXELS * APERTURE_RADIUS_PIXELS &&
                signal > 0) {
                dx = x - star->x;
                dy = y - star->y;
                bin = (unsigned)(sqrt(dx * dx + dy * dy) * RADIAL_BINS_PER_PIXEL);
                if (bin >= RADIAL_BINS)
                    bin = RADIAL_BINS - 1;
                radial[bin] += signal;
            }
        }
    }
    for (bin = 0; bin < RADIAL_BINS; ++bin) {
        double previous = cumulative;
        cumulative += radial[bin];
        if (cumulative >= flux * 0.5) {
            star->hfd = 2 * (bin + (flux * 0.5 - previous) / radial[bin])
                / RADIAL_BINS_PER_PIXEL;
            break;
        }
    }
    /* Reject single-pixel peaks and sources outside the useful HFD range. */
    return star->hfd >= MIN_HFD_PIXELS && star->hfd <= MAX_HFD_PIXELS;
}

/* Validate before shortlisting so bright extended objects cannot consume slots. */
static int discover(const uint16_t *pixels, uint32_t width, uint32_t height,
    double sky, double mad, double minimum_signal, Candidate *candidates)
{
    size_t count = (size_t)width * height, seed, head, tail, capacity = 1024;
    size_t *queue;
    unsigned char *visited;
    int selected = 0;

    visited = calloc(count, 1);
    queue = malloc(capacity * sizeof(*queue));
    if (!visited || !queue) {
        free(visited);
        free(queue);
        return -2;
    }

    /* Discovery flux includes only the above-threshold component pixels. */
    for (seed = 0; seed < count; ++seed) {
        double flux = 0, sx = 0, sy = 0;
        if (visited[seed] || pixels[seed] <= sky + minimum_signal)
            continue;
        head = 0;
        tail = 1;
        queue[0] = seed;
        visited[seed] = 1;
        while (head < tail) {
            size_t index = queue[head++];
            int x = (int)(index % width), y = (int)(index / width);
            int nx, ny;
            double signal = pixels[index] - sky;
            flux += signal;
            sx += signal * x;
            sy += signal * y;
            for (ny = y > 0 ? y - 1 : 0;
                 ny <= y + 1 && ny < (int)height; ++ny) {
                for (nx = x > 0 ? x - 1 : 0;
                     nx <= x + 1 && nx < (int)width; ++nx) {
                    size_t neighbour = (size_t)ny * width + nx;
                    if (visited[neighbour] || pixels[neighbour] <= sky + minimum_signal)
                        continue;
                    if (tail == capacity) {
                        size_t next = capacity <= SIZE_MAX / 2 / sizeof(*queue)
                            ? capacity * 2 : 0;
                        size_t *grown = next ? realloc(queue, next * sizeof(*queue)) : NULL;
                        if (!grown) {
                            free(queue);
                            free(visited);
                            return -2;
                        }
                        queue = grown;
                        capacity = next;
                    }
                    visited[neighbour] = 1;
                    queue[tail++] = neighbour;
                }
            }
        }
        if (tail >= MIN_COMPONENT_PIXELS) {
            Candidate candidate;
            candidate.flux = flux;
            if (measure(pixels, width, height, sky, sx / flux, sy / flux,
                        &candidate.star)) {
                quality(&candidate.star, sky, mad);
                retain(candidates, &selected, MAX_STARS, candidate);
            }
        }
    }
    free(queue);
    free(visited);
    return selected;
}

int CALL AnalyzeStars(const uint16_t *pixels, uint32_t width,
    uint32_t height, uint32_t median, uint32_t mad,
    StarRecord *stars)
{
    Candidate candidates[MAX_STARS];
    int selected, i;

    if (!stars)
        return -1;
    memset(stars, 0, sizeof(*stars));
    if (!pixels || !width || !height ||
        width > INT32_MAX || height > INT32_MAX ||
        (size_t)width > SIZE_MAX / height / sizeof(*pixels))
        return -1;
    selected = discover(pixels, width, height, median, mad,
        fmax(MIN_SIGNAL_INCREMENT_ADU, MAD_MULTIPLIER * (double)mad), candidates);
    if (selected < 0)
        return selected;
    for (i = 0; i < selected; ++i)
        stars[i] = candidates[i].star;
    qsort(stars, selected, sizeof(*stars), star_order);
    memset(stars + selected, 0, sizeof(*stars));
    return selected;
}
