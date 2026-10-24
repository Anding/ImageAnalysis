#include "image_analysis.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

enum { APERTURE = 16, RADIAL_SCALE = 16,
    RADIAL_BINS = 2 * APERTURE * RADIAL_SCALE + 1 };

static double background(const uint32_t *histogram, size_t count)
{
    uint64_t cumulative = 0;
    uint32_t bin;
    for (bin = 0; bin < 65536; ++bin) {
        cumulative += histogram[bin];
        if (cumulative > count / 2)
            return bin;
    }
    return 65535;
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

int CALL AnalyzeStars(const uint16_t *pixels, uint32_t width,
    uint32_t height, const uint32_t *histogram, StarRecord *stars)
{
    size_t count, seed, head, tail, capacity = 1024;
    size_t *queue;
    unsigned char *visited;
    double sky;
    int found = 0;

    if (!stars)
        return -1;
    memset(stars, 0, sizeof(*stars));
    if (!pixels || !histogram || !width || !height ||
        width > INT32_MAX || height > INT32_MAX ||
        (size_t)width > SIZE_MAX / height)
        return -1;
    count = (size_t)width * height;
    sky = background(histogram, count);
    visited = calloc(count, 1);
    queue = malloc(capacity * sizeof(*queue));
    if (!visited || !queue) {
        free(visited);
        free(queue);
        return -2;
    }

    /* One eight-connected component per candidate; no deblending or ranking. */
    for (seed = 0; seed < count && found < MAX_STARS; ++seed) {
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
                            memset(stars, 0, sizeof(*stars));
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
        if (tail >= 3 &&
            measure(pixels, width, height, sky, sx / flux, sy / flux, stars + found))
            ++found;
    }
    memset(stars + found, 0, sizeof(*stars));
    free(queue);
    free(visited);
    return found;
}
