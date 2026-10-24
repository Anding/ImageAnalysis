#include "image_analysis.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

enum { APERTURE = 16, RADIAL_SCALE = 16,
    RADIAL_BINS = 2 * APERTURE * RADIAL_SCALE + 1 };

typedef struct Candidate {
    double x, y, flux;
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

/* Keep only the best candidates in a small, ordered, bounded array. */
static void retain(Candidate *candidates, int *count, int limit, Candidate candidate)
{
    int i = *count;
    if (i == limit) {
        Candidate last = candidates[i - 1];
        if (order(candidate.flux, candidate.x, candidate.y,
                  last.flux, last.x, last.y) >= 0)
            return;
        --i;
    } else {
        ++*count;
    }
    while (i > 0) {
        Candidate previous = candidates[i - 1];
        if (order(candidate.flux, candidate.x, candidate.y,
                  previous.flux, previous.x, previous.y) >= 0)
            break;
        candidates[i] = previous;
        --i;
    }
    candidates[i] = candidate;
}

static int inside(double x, double y, uint32_t width, uint32_t height)
{
    return x >= APERTURE && y >= APERTURE &&
        x + APERTURE <= width - 1.0 && y + APERTURE <= height - 1.0;
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
    left = (int)ceil(cx - APERTURE);
    right = (int)floor(cx + APERTURE);
    top = (int)ceil(cy - APERTURE);
    bottom = (int)floor(cy + APERTURE);
    for (y = top; y <= bottom; ++y) {
        for (x = left; x <= right; ++x) {
            double dx = x - cx, dy = y - cy;
            double signal = pixels[(size_t)y * width + x] - sky;
            if (dx * dx + dy * dy <= APERTURE * APERTURE && signal > 0) {
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

    left = (int)ceil(cx - APERTURE);
    right = (int)floor(cx + APERTURE);
    top = (int)ceil(cy - APERTURE);
    bottom = (int)floor(cy + APERTURE);
    memset(star, 0, sizeof(*star));
    flux = sx = sy = 0;
    for (y = top; y <= bottom; ++y) {
        for (x = left; x <= right; ++x) {
            double dx = x - cx, dy = y - cy;
            double r2 = dx * dx + dy * dy;
            uint16_t value;
            double signal;
            if (r2 > APERTURE * APERTURE)
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
            if (dx * dx + dy * dy <= APERTURE * APERTURE && signal > 0) {
                dx = x - star->x;
                dy = y - star->y;
                bin = (unsigned)(sqrt(dx * dx + dy * dy) * RADIAL_SCALE);
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
                / RADIAL_SCALE;
            break;
        }
    }
    return 1;
}

/* Discover throughout the bitmap; aperture metrics are a separate phase. */
static int discover(const uint16_t *pixels, uint32_t width, uint32_t height,
    double sky, int limit, Candidate *candidates)
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
        if (visited[seed] || pixels[seed] <= sky + 16)
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
                    if (visited[neighbour] || pixels[neighbour] <= sky + 16)
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
        if (tail >= 3) {
            Candidate candidate = {sx / flux, sy / flux, flux};
            if (inside(candidate.x, candidate.y, width, height))
                retain(candidates, &selected, limit, candidate);
        }
    }
    free(queue);
    free(visited);
    return selected;
}

int CALL AnalyzeStars(const uint16_t *pixels, uint32_t width,
    uint32_t height, uint32_t background, uint32_t max_stars,
    StarRecord *stars)
{
    Candidate candidates[MAX_STARS];
    int selected, i, found = 0;

    if (!stars)
        return -1;
    memset(stars, 0, sizeof(*stars));
    if (!pixels || !width || !height ||
        width > INT32_MAX || height > INT32_MAX || max_stars > MAX_STARS ||
        (size_t)width > SIZE_MAX / height)
        return -1;
    if (max_stars == 0)
        return 0;
    selected = discover(pixels, width, height, background, (int)max_stars, candidates);
    if (selected < 0)
        return selected;
    for (i = 0; i < selected; ++i) {
        if (measure(pixels, width, height, background, candidates[i].x, candidates[i].y,
                    stars + found))
            ++found;
    }
    qsort(stars, found, sizeof(*stars), star_order);
    memset(stars + found, 0, sizeof(*stars));
    return found;
}
