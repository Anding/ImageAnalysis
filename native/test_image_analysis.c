/* Build with image_analysis.c; no equipment, files or third-party libraries. */
#include "image_analysis.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(test) do { if (!(test)) { \
    fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #test); exit(1); \
} } while (0)
#define NEAR(value, expected, tolerance) CHECK(fabs((value) - (expected)) < (tolerance))

static StarRecord stars[MAX_STARS + 2];

static void fill(uint16_t *pixels, unsigned count, uint16_t value)
{
    unsigned i;
    for (i = 0; i < count; ++i)
        pixels[i] = value;
}

static int analyze_background(uint16_t *pixels, unsigned width, unsigned height,
    unsigned background, unsigned limit)
{
    int result;
    memset(stars, 0xa5, sizeof(stars));
    result = AnalyzeStars(pixels, width, height, background, limit, stars);
    CHECK(result >= 0 && (unsigned)result <= limit);
    CHECK(stars[result].positive_pixel_count == 0);
    return result;
}

static int analyze_limit(uint16_t *pixels, unsigned width, unsigned height, unsigned limit)
{
    return analyze_background(pixels, width, height, 1000, limit);
}

static int analyze(uint16_t *pixels, unsigned width, unsigned height)
{
    return analyze_limit(pixels, width, height, MAX_STARS);
}

static void small_sources(void)
{
    uint16_t pixels[64 * 64];
    StarRecord saved;
    fill(pixels, 64 * 64, 1000);
    CHECK(analyze(pixels, 64, 64) == 0);

    pixels[32 * 64 + 32] = 1100;
    pixels[32 * 64 + 33] = 1200;
    pixels[33 * 64 + 32] = 1100;
    CHECK(analyze(pixels, 64, 64) == 1);
    CHECK(stars[0].positive_pixel_count == 3 && stars[0].max_pixel == 1200);
    CHECK(stars[0].saturated_pixel_count == 0);
    NEAR(stars[0].x, 32.5, 1e-12);
    NEAR(stars[0].y, 32.25, 1e-12);
    NEAR(stars[0].flux, 400, 1e-12);
    NEAR(stars[0].qxx, 0.25, 1e-12);
    NEAR(stars[0].qyy, 0.1875, 1e-12);
    NEAR(stars[0].qxy, -0.125, 1e-12);
    NEAR(stars[0].hfd, 2 * sqrt(0.3125), 0.125);
    saved = stars[0];
    CHECK(analyze(pixels, 64, 64) == 1);
    CHECK(memcmp(&saved, stars, sizeof(saved)) == 0);
    CHECK(analyze_background(pixels, 64, 64, 1100, MAX_STARS) == 0);

    pixels[32 * 64 + 33] = 65535;
    CHECK(analyze(pixels, 64, 64) == 1);
    CHECK(stars[0].saturated_pixel_count == 1 && stars[0].max_pixel == 65535);
    NEAR(stars[0].flux, 64735, 1e-12);

    /* A fixed 16-ADU threshold, minimum three eight-connected pixels. */
    fill(pixels, 64 * 64, 1000);
    pixels[30 * 64 + 30] = 1016;
    pixels[31 * 64 + 31] = 1016;
    pixels[32 * 64 + 32] = 1016;
    CHECK(analyze(pixels, 64, 64) == 0);
    pixels[30 * 64 + 30] = 1017;
    pixels[31 * 64 + 31] = 1017;
    CHECK(analyze(pixels, 64, 64) == 0);
    pixels[32 * 64 + 32] = 1017;
    CHECK(analyze(pixels, 64, 64) == 1);
    NEAR(stars[0].x, 31, 1e-12);
    NEAR(stars[0].qxx, stars[0].qyy, 1e-12);
    CHECK(stars[0].qxy > 0); /* equal diagonals do not imply roundness */

    fill(pixels, 64 * 64, 1000);
    pixels[32 * 64 + 14] = pixels[32 * 64 + 15] = pixels[32 * 64 + 16] = 1100;
    CHECK(analyze(pixels, 64, 64) == 0);
    fill(pixels, 64 * 64, 1000);
    pixels[32 * 64 + 15] = pixels[32 * 64 + 16] = pixels[32 * 64 + 17] = 1100;
    CHECK(analyze(pixels, 64, 64) == 1); /* aperture just touches x=0 */

    fill(pixels, 64 * 64, 65535);
    CHECK(analyze_background(pixels, 64, 64, 65535, MAX_STARS) == 0);
    fill(pixels, 64 * 64, 0);
    CHECK(analyze_background(pixels, 64, 64, 0, MAX_STARS) == 0);
    CHECK(analyze(pixels, 1, 1) == 0);
}

static void gaussian(void)
{
    uint16_t pixels[128 * 128];
    unsigned x, y;
    for (y = 0; y < 128; ++y) {
        for (x = 0; x < 128; ++x) {
            double dx = x - 64.25, dy = y - 63.75;
            pixels[y * 128 + x] = (uint16_t)(1000 + 20000 * exp(-(dx * dx + dy * dy) / 8));
        }
    }
    CHECK(analyze(pixels, 128, 128) == 1);
    NEAR(stars[0].x, 64.25, 0.002);
    NEAR(stars[0].y, 63.75, 0.002);
    NEAR(stars[0].qxx, 4, 0.02);
    NEAR(stars[0].qyy, 4, 0.02);
    NEAR(stars[0].qxy, 0, 0.01);
    NEAR(stars[0].hfd, 4 * sqrt(2 * log(2.0)), 0.25);

    /* Exercise growth of the connected-component workspace beyond 1024 pixels. */
    fill(pixels, 128 * 128, 1000);
    for (y = 40; y < 88; ++y)
        for (x = 40; x < 88; ++x)
            pixels[y * 128 + x] = 1100;
    CHECK(analyze(pixels, 128, 128) == 1);
    NEAR(stars[0].x, 63.5, 1e-12);
    NEAR(stars[0].y, 63.5, 1e-12);
}

static void limit_and_independence(void)
{
    enum { WIDTH = 704 };
    uint16_t *pixels = malloc(WIDTH * WIDTH * sizeof(*pixels));
    StarRecord other[MAX_STARS + 1], guard;
    unsigned x, y;
    CHECK(pixels != NULL);
    fill(pixels, WIDTH * WIDTH, 1000);
    for (y = 24; y <= 664; y += 40) {
        for (x = 24; x <= 664; x += 40) {
            pixels[y * WIDTH + x - 1] = 1100;
            pixels[y * WIDTH + x] = (uint16_t)(1200 + (y - 24) / 40 * 17 + (x - 24) / 40);
            pixels[y * WIDTH + x + 1] = 1100;
        }
    }
    CHECK(analyze(pixels, WIDTH, WIDTH) == 256);
    memset(&guard, 0xa5, sizeof(guard));
    CHECK(memcmp(stars + MAX_STARS + 1, &guard, sizeof(guard)) == 0);
    NEAR(stars[0].x, 664, 1e-12); /* brightest star is the last discovered */
    NEAR(stars[0].y, 664, 1e-12);
    NEAR(stars[255].x, 664, 1e-12);
    NEAR(stars[255].y, 64, 1e-12);
    for (x = 1; x < MAX_STARS; ++x)
        CHECK(stars[x - 1].flux > stars[x].flux);
    fill(pixels, WIDTH * WIDTH, 2000);
    CHECK(AnalyzeStars(pixels, WIDTH, WIDTH, 2000, 10, other) == 0);
    NEAR(stars[255].y, 64, 1e-12); /* no process-global result storage */
    free(pixels);
}

static void source(uint16_t *pixels, unsigned width, unsigned x, unsigned y, unsigned signal)
{
    pixels[y * width + x - 1] = (uint16_t)(1000 + signal);
    pixels[y * width + x] = (uint16_t)(1000 + 2 * signal);
    pixels[y * width + x + 1] = (uint16_t)(1000 + signal);
}

static void selection_and_order(void)
{
    enum { WIDTH = 192, HEIGHT = 160 };
    uint16_t pixels[WIDTH * HEIGHT];
    StarRecord guard;
    unsigned x, y;

    fill(pixels, WIDTH * HEIGHT, 1000);
    source(pixels, WIDTH, 32, 32, 100);
    source(pixels, WIDTH, 112, 32, 200);
    source(pixels, WIDTH, 64, 96, 300);
    source(pixels, WIDTH, 8, 96, 1000); /* bright edge candidate is omitted */
    CHECK(analyze_limit(pixels, WIDTH, HEIGHT, 1) == 1);
    NEAR(stars[0].x, 64, 1e-12);
    NEAR(stars[0].y, 96, 1e-12);
    memset(&guard, 0xa5, sizeof(guard));
    CHECK(memcmp(stars + 2, &guard, sizeof(guard)) == 0);
    CHECK(analyze_limit(pixels, WIDTH, HEIGHT, 2) == 2);
    NEAR(stars[0].flux, 1200, 1e-12);
    NEAR(stars[1].flux, 800, 1e-12);
    CHECK(analyze_limit(pixels, WIDTH, HEIGHT, 10) == 3);
    CHECK(analyze_limit(pixels, WIDTH, HEIGHT, 0) == 0);
    CHECK(AnalyzeStars(pixels, WIDTH, HEIGHT, 1000, 257, stars) == -1);
    CHECK(stars[0].positive_pixel_count == 0);

    /* Equal discovery/measured flux: X precedes Y, even at the selection cut. */
    fill(pixels, WIDTH * HEIGHT, 1000);
    source(pixels, WIDTH, 112, 32, 100);
    source(pixels, WIDTH, 32, 32, 100);
    source(pixels, WIDTH, 32, 96, 100);
    CHECK(analyze_limit(pixels, WIDTH, HEIGHT, 3) == 3);
    NEAR(stars[0].x, 32, 1e-12);
    NEAR(stars[0].y, 32, 1e-12);
    NEAR(stars[1].x, 32, 1e-12);
    NEAR(stars[1].y, 96, 1e-12);
    NEAR(stars[2].x, 112, 1e-12);
    CHECK(analyze_limit(pixels, WIDTH, HEIGHT, 2) == 2);
    NEAR(stars[1].x, 32, 1e-12);
    NEAR(stars[1].y, 96, 1e-12);

    /* Sub-threshold wings can reverse measured flux, but not candidate selection. */
    fill(pixels, WIDTH * HEIGHT, 1000);
    for (y = 28; y <= 36; ++y)
        for (x = 28; x <= 36; ++x)
            pixels[y * WIDTH + x] = 1010;
    source(pixels, WIDTH, 32, 32, 100);
    source(pixels, WIDTH, 112, 32, 110);
    CHECK(analyze_limit(pixels, WIDTH, HEIGHT, 1) == 1);
    NEAR(stars[0].x, 112, 1e-12);
    CHECK(analyze_limit(pixels, WIDTH, HEIGHT, 2) == 2);
    NEAR(stars[0].x, 32, 1e-12);
    NEAR(stars[0].flux, 1180, 1e-12);
    NEAR(stars[1].flux, 440, 1e-12);
}

int main(void)
{
    CHECK(sizeof(StarRecord) == 72 && offsetof(StarRecord, x) == 8);
    CHECK(offsetof(StarRecord, max_pixel) == 64);
    small_sources();
    gaussian();
    limit_and_independence();
    selection_and_order();
    CHECK(AnalyzeStars(NULL, 64, 64, 1000, 10, stars) == -1);
    CHECK(stars[0].positive_pixel_count == 0);
    puts("Star analysis native regression passed");
    return 0;
}
