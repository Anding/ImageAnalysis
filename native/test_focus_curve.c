/* Link with focus_curve.c alone: no bitmap analyzer or equipment dependencies. */
#include "focus_curve.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(test) do { if (!(test)) { \
    fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #test); exit(1); \
} } while (0)

static StarRecord records[MAX_STARS + 1];

static StarRecord star(double x, double y, double hfd, uint16_t peak)
{
    StarRecord result = {0};
    result.positive_pixel_count = 10;
    result.x = x;
    result.y = y;
    result.hfd = hfd;
    result.flux = 100000;
    result.max_pixel = peak;
    return result;
}

static double vcurve(double position, double focus)
{
    return 2 + (position < focus ? 0.02 * (focus - position) :
                                  0.03 * (position - focus));
}

static void failed_finish(void)
{
    int32_t focus = -123;
    uint32_t error = 456;
    CHECK(finish_computation(&focus, &error) == -1);
    CHECK(focus == -123 && error == 456);
}

static void lifecycle(void)
{
    uint32_t exposure = 1000000;
    memset(records, 0, sizeof(records));
    CHECK(add_focus_point(records, 5000, &exposure) == -1);
    failed_finish();
    CHECK(start_computation(0) == -1);
    CHECK(start_computation(MAX_STARS + 1) == -1);
    CHECK(start_computation(1) == 0);
    CHECK(add_focus_point(records, 5000, &exposure) == -1);
    records[0] = star(32, 32, 5, 30000);
    CHECK(add_focus_point(records, 5000, NULL) == -1);
    CHECK(add_focus_point(NULL, 5000, &exposure) == -1);
    exposure = 0;
    CHECK(add_focus_point(records, 5000, &exposure) == -1);
    exposure = 1000000;
    records[0].hfd = NAN;
    CHECK(add_focus_point(records, 5000, &exposure) == -1);
    records[0].hfd = 5;
    CHECK(add_focus_point(records, 5000, &exposure) == 0);
    CHECK(start_computation(1) == 0); /* replaces, rather than appends to, a run */
    failed_finish();
    CHECK(add_focus_point(records, 5000, &exposure) == -1);
    CHECK(start_computation(1) == 0);
    CHECK(finish_computation(NULL, NULL) == -1);
    failed_finish();
}

static void exposure_guidance(void)
{
    uint32_t exposure;
    int i;
    memset(records, 0, sizeof(records));
    CHECK(start_computation(1) == 0);
    records[0] = star(32, 32, 5, 60000);
    exposure = 1000000;
    CHECK(add_focus_point(records, 5000, &exposure) == 0);
    CHECK(exposure == 1000000);
    records[0].max_pixel = 5000;
    exposure = 1000000;
    CHECK(add_focus_point(records, 5000, &exposure) == 0);
    CHECK(exposure == 1000000);
    records[0].max_pixel = 65535;
    exposure = 1000000;
    CHECK(add_focus_point(records, 5000, &exposure) == 0);
    CHECK(exposure == 1000000);
    records[0].max_pixel = 30000;
    records[0].saturated_pixel_count = 1;
    exposure = 1000000;
    CHECK(add_focus_point(records, 5000, &exposure) == 0);
    CHECK(exposure == 1000000);
    records[0].saturated_pixel_count = 0;
    records[0].max_pixel = 5000;
    exposure = UINT32_MAX;
    CHECK(add_focus_point(records, 5000, &exposure) == 0);
    CHECK(exposure == UINT32_MAX);
    records[0].max_pixel = 65535;
    records[0].saturated_pixel_count = 4;
    exposure = 1;
    CHECK(add_focus_point(records, 5000, &exposure) == -1);
    CHECK(exposure == 1);
    exposure = 150000;
    CHECK(add_focus_point(records, 5000, &exposure) == 1);
    CHECK(exposure == 100000);
    CHECK(add_focus_point(records, 5000, &exposure) == -1);
    CHECK(exposure == 100000);
    failed_finish(); /* repeated positions cannot form a focus curve */

    CHECK(start_computation(1) == 0);
    for (i = 0; i < 4; ++i) {
        int32_t position = 4800 + 100 * i;
        records[0] = star(32, 32, vcurve(position, 4950), i % 2 ? 10000 : 50000);
        exposure = 1000000;
        CHECK(add_focus_point(records, position, &exposure) == 0);
        CHECK(exposure == 1000000);
    }
    /* A rejected fifth position still leaves too few accepted positions. */
    records[0] = star(32, 32, 12, 60000);
    records[0].saturated_pixel_count = 4;
    CHECK(add_focus_point(records, 5200, &exposure) == 1);
    failed_finish();
}

static void single_star(double wanted, int32_t centre)
{
    int32_t focus;
    uint32_t error, exposure = 1000000;
    int i;
    CHECK(start_computation(1) == 0);
    for (i = 0; i < 7; ++i) {
        int32_t position = centre - 300 + i * 100;
        memset(records, 0, sizeof(records));
        records[0] = star(32 + i * 0.25, 32, vcurve(position, wanted), 30000);
        if (i == 3) {
            uint32_t retry_exposure = exposure;
            records[0].max_pixel = 65535;
            records[0].saturated_pixel_count = 4;
            records[0].hfd = 12; /* would spoil the curve if stored */
            CHECK(add_focus_point(records, position, &retry_exposure) == 1);
            records[0].max_pixel = 30000;
            records[0].saturated_pixel_count = 0;
            records[0].hfd = vcurve(position, wanted);
        }
        CHECK(add_focus_point(records, position, &exposure) == 0);
        memset(records, 0, sizeof(records)); /* the computation owns its samples */
    }
    CHECK(finish_computation(&focus, &error) == 0);
    CHECK(focus == (int32_t)round(wanted));
    CHECK(error == 1);
    failed_finish();
}

static void reordered_identities(void)
{
    int32_t focus;
    uint32_t error, exposure = 1000000;
    int i, j;
    CHECK(start_computation(3) == 0);
    for (i = 0; i < 9; ++i) {
        int32_t position = 4600 + 100 * i;
        memset(records, 0, sizeof(records));
        for (j = 0; j < 3; ++j) {
            int index = i % 2 ? 2 - j : j;
            records[index] = star(32 + j * 40 + i * 0.2, 40,
                vcurve(position, 4990 + j * 10), 30000);
        }
        CHECK(add_focus_point(records, position, &exposure) == 0);
    }
    CHECK(finish_computation(&focus, &error) == 0);
    CHECK(focus == 5000);
    CHECK(error == 15); /* ceil(1.4826 * 10), not zero from individually exact fits */
}

static void missing_and_ambiguous(void)
{
    int32_t focus;
    uint32_t error, exposure = 1000000;
    int i;
    CHECK(start_computation(2) == 0);
    memset(records, 0, sizeof(records));
    records[0] = star(32, 32, vcurve(4700, 5000), 30000);
    records[1] = star(36, 32, vcurve(4700, 5100), 30000);
    CHECK(add_focus_point(records, 4700, &exposure) == 0);
    memset(records, 0, sizeof(records));
    records[0] = star(34, 32, 5, 5000); /* equally near both known identities */
    CHECK(add_focus_point(records, 4800, &exposure) == -1);
    CHECK(exposure == 1000000); /* ambiguity is not underexposure */
    records[0] = star(100, 32, 5, 5000); /* outside tracking radius */
    CHECK(add_focus_point(records, 4800, &exposure) == -1);
    for (i = 0; i < 7; ++i) {
        int32_t position = 4800 + 100 * i;
        memset(records, 0, sizeof(records));
        records[0] = star(150, 32, 50, 65535); /* untracked: do not use for exposure */
        records[1] = star(36, 32, vcurve(position, 5100), 30000);
        CHECK(add_focus_point(records, position, &exposure) == 0);
    }
    CHECK(finish_computation(&focus, &error) == 0);
    CHECK(focus == 5100 && error == 1); /* star 1 never became missing star 0 */
}

static void retry_tracking(void)
{
    uint32_t exposure = 1000000;
    memset(records, 0, sizeof(records));
    CHECK(start_computation(1) == 0);
    records[0] = star(100, 32, 5, 65535);
    records[0].saturated_pixel_count = 4;
    CHECK(add_focus_point(records, 5000, &exposure) == 1);
    records[0] = star(32, 32, 5, 30000);
    CHECK(add_focus_point(records, 5000, &exposure) == 0);
    records[0] = star(38, 32, 12, 60000);
    records[0].saturated_pixel_count = 4;
    CHECK(add_focus_point(records, 5100, &exposure) == 1);
    records[0] = star(26, 32, 4, 30000);
    CHECK(add_focus_point(records, 5100, &exposure) == 0);
    failed_finish(); /* neither rejected frame established or moved an identity */
}

static void even_combination(void)
{
    int32_t focus;
    uint32_t error, exposure = 1000000;
    int i;
    CHECK(start_computation(2) == 0);
    for (i = 0; i < 7; ++i) {
        int32_t position = 4700 + 100 * i;
        memset(records, 0, sizeof(records));
        records[0] = star(32, 32, vcurve(position, 5000), 30000);
        records[1] = star(80, 32, vcurve(position, 5020), 30000);
        CHECK(add_focus_point(records, position, &exposure) == 0);
    }
    CHECK(finish_computation(&focus, &error) == 0);
    CHECK(focus == 5010 && error == 15);
}

static void late_and_clipped_detections(void)
{
    int32_t focus;
    uint32_t error, exposure = 1000000;
    int i;
    CHECK(start_computation(2) == 0);
    for (i = 0; i < 9; ++i) {
        int32_t position = 4600 + i * 100;
        memset(records, 0, sizeof(records));
        records[0] = star(32, 32, 10, 30000); /* compact, focus-independent object */
        if (i > 0)
            records[1] = star(80, 32, vcurve(position, 5000), 30000);
        CHECK(add_focus_point(records, position, &exposure) == 0);
    }
    CHECK(finish_computation(&focus, &error) == 0);
    CHECK(focus == 5000 && error == 1);

    CHECK(start_computation(1) == 0);
    for (i = 0; i < 7; ++i) {
        int32_t position = 4700 + i * 100;
        memset(records, 0, sizeof(records));
        records[0] = star(33, 32, vcurve(position, 5000), 30000);
        if (i > 0) {
            records[1] = records[0];
            records[1].x = 32;
            records[0].saturated_pixel_count = 4;
        }
        CHECK(add_focus_point(records, position, &exposure) == 0);
        CHECK(exposure == 1000000);
    }
    CHECK(finish_computation(&focus, &error) == 0);
    CHECK(focus == 5000 && error == 1);
}

static void contaminated_field(void)
{
    int32_t focus;
    uint32_t error, exposure = 100000;
    int i;
    CHECK(start_computation(3) == 0);
    for (i = 0; i < 7; ++i) {
        int32_t position = 4700 + i * 100;
        memset(records, 0, sizeof(records));
        records[0] = star(32, 32, 0.1, 65535); /* isolated hot pixel */
        records[0].saturated_pixel_count = 1;
        records[1] = star(80, 32, 5, 65535); /* clipped bright star */
        records[1].saturated_pixel_count = 4;
        records[2] = star(120, 32, vcurve(position, 5000), 1500);
        records[2].saturated_pixel_count = i == 3 ? 3 : 0;
        CHECK(add_focus_point(records, position, &exposure) == 0);
        CHECK(exposure == 100000);
    }
    CHECK(finish_computation(&focus, &error) == 0);
    CHECK(focus == 5000 && error == 1);
}

static void robust_combination(void)
{
    int32_t focus;
    uint32_t error, exposure = 1000000;
    const double centres[] = {5000, 5002, 5500};
    int i, j;
    CHECK(start_computation(3) == 0);
    for (i = 0; i < 11; ++i) {
        int32_t position = 4200 + 200 * i;
        memset(records, 0, sizeof(records));
        for (j = 0; j < 3; ++j)
            records[j] = star(32 + j * 40, 32, vcurve(position, centres[j]) / 2, 30000);
        CHECK(add_focus_point(records, position, &exposure) == 0);
    }
    CHECK(finish_computation(&focus, &error) == 0);
    CHECK(focus == 5002);
    CHECK(error == 3); /* one outlying star does not dominate the median/MAD */
}

static void fit_quality(void)
{
    const double noise[] = {0.3, -0.2, 0.1, -0.1, 0.2, -0.3, 0.1};
    const int32_t repeats[] = {4800, 4800, 4900, 4900, 5200, 5200};
    int32_t focus;
    uint32_t error, exposure = 1000000;
    int i;
    CHECK(start_computation(1) == 0);
    for (i = 0; i < 7; ++i) {
        int32_t position = 4700 + i * 100;
        memset(records, 0, sizeof(records));
        records[0] = star(32, 32, vcurve(position, 5050) + noise[i], 30000);
        CHECK(add_focus_point(records, position, &exposure) == 0);
    }
    CHECK(finish_computation(&focus, &error) == 0);
    CHECK(abs(focus - 5050) < 25 && error > 1); /* residuals contribute to error */

    CHECK(start_computation(1) == 0);
    for (i = 0; i < 6; ++i) {
        records[0] = star(32, 32, 10 - i, 30000); /* one descending branch */
        CHECK(add_focus_point(records, 4700 + i * 100, &exposure) == 0);
    }
    failed_finish();
    CHECK(start_computation(1) == 0);
    for (i = 0; i < 7; ++i) {
        int32_t position = 4700 + i * 100;
        records[0] = star(32, 32, 10 + 0.01 * vcurve(position, 5000), 30000);
        CHECK(add_focus_point(records, position, &exposure) == 0);
    }
    failed_finish(); /* an extended, almost focus-independent source is not a star curve */
    CHECK(start_computation(1) == 0);
    for (i = 0; i < 6; ++i) {
        records[0] = star(32, 32, 5, 30000); /* no focus dependence */
        CHECK(add_focus_point(records, 4700 + i * 100, &exposure) == 0);
    }
    failed_finish();
    CHECK(start_computation(1) == 0);
    for (i = 0; i < 6; ++i) {
        records[0] = star(32, 32, vcurve(repeats[i], 5000), 30000);
        CHECK(add_focus_point(records, repeats[i], &exposure) == 0);
    }
    failed_finish(); /* repeated frames cannot invent five distinct positions */
}

static void limits_and_growth(void)
{
    int32_t focus;
    uint32_t error, exposure = 1000000;
    int i, j;
    CHECK(start_computation(MAX_STARS) == 0);
    for (i = 0; i < 21; ++i) { /* grow sample storage past its initial capacity */
        int32_t position = 4000 + 100 * i;
        memset(records, 0, sizeof(records));
        for (j = 0; j < MAX_STARS; ++j)
            records[j] = star(32 + j * 20, 32, vcurve(position, 5050) / 2, 30000);
        CHECK(add_focus_point(records, position, &exposure) == 0);
    }
    CHECK(finish_computation(&focus, &error) == 0);
    CHECK(focus == 5050 && error == 1);

    CHECK(start_computation(1) == 0);
    for (i = 0; i < 7; ++i) {
        int32_t position = 4700 + 100 * i;
        memset(records, 0, sizeof(records));
        records[0] = star(32, 32, vcurve(position, 5050), 30000);
        records[1] = star(80, 32, 100, 65535); /* beyond the selected limit */
        CHECK(add_focus_point(records, position, &exposure) == 0);
    }
    CHECK(finish_computation(&focus, &error) == 0);
    CHECK(focus == 5050 && error == 1);
}

int main(void)
{
    lifecycle();
    exposure_guidance();
    single_star(5050, 5000);
    single_star(5000.4, 5000);
    single_star(-5000.6, -5000);
    single_star(2147483000.4, 2147483000);
    reordered_identities();
    missing_and_ambiguous();
    retry_tracking();
    even_combination();
    late_and_clipped_detections();
    contaminated_field();
    robust_combination();
    fit_quality();
    limits_and_growth();
    puts("Independent focus-curve native regression passed");
    return 0;
}
