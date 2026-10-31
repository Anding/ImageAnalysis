#include "focus_curve.h"
#include "analysis_parameters.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define CHECK(test) do { if (!(test)) { \
    fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #test); exit(1); \
} } while (0)

static StarRecord records[MAX_STARS + 1];

static StarRecord star(double x, double y, double hfd, double flux,
    double roundness, double over, double under, double scale)
{
    StarRecord result = {0};
    result.positive_pixel_count = 10;
    result.x = x;
    result.y = y;
    result.hfd = hfd;
    result.flux = flux;
    result.roundness = roundness;
    result.over_exposed_risk = over;
    result.under_exposed_risk = under;
    result.exposure_scale = scale;
    return result;
}

static double hyperbola(double position, double focus)
{
    double offset = 0.02 * (position - focus);
    return sqrt(2.0 * 2.0 + offset * offset);
}

static void reset_records(void)
{
    memset(records, 0, sizeof(records));
}

static void select_one(double x, uint32_t *exposure)
{
    reset_records();
    records[0] = star(x, 32, 5, 100000, 1, 0, 0, 1);
    CHECK(start_computation(1) == 0);
    CHECK(select_focus_stars(records, exposure) == 0);
}

static void add_and_commit(int32_t position, double hfd, double over,
    double under, double scale, int expected_status)
{
    uint32_t exposure = 1000000;
    reset_records();
    records[0] = star(32, 32, hfd, 100000, 1, over, under, scale);
    CHECK(add_focus_point(records, position, &exposure) == expected_status);
    CHECK(finish_focus_point() == 0);
}

static void require_failed_fit(void)
{
    int32_t focus = -123456;
    uint32_t error = 0;
    int status = finish_computation(&focus, &error);
    CHECK(status == -1);
    CHECK(focus == -123456);
    CHECK(error == UINT32_MAX);
}

static void layout_and_lifecycle(void)
{
    uint32_t exposure = 1000000;
    int32_t focus = 20;
    uint32_t error = 0;
    CHECK(sizeof(StarRecord) == 104);
    CHECK(offsetof(StarRecord, roundness) == 72);
    CHECK(offsetof(StarRecord, over_exposed_risk) == 80);
    CHECK(offsetof(StarRecord, under_exposed_risk) == 88);
    CHECK(offsetof(StarRecord, exposure_scale) == 96);
    CHECK(add_focus_point(records, 0, &exposure) == -1);
    CHECK(start_computation(0) == -1);
    CHECK(start_computation(MAX_STARS + 1) == -1);
    CHECK(start_computation(1) == 0);
    CHECK(select_focus_stars(records, &exposure) == -1);
    CHECK(select_focus_stars(records, &exposure) == -1);
    CHECK(add_focus_point(records, 0, &exposure) == -1);
    CHECK(finish_computation(&focus, &error) == -1);
    CHECK(focus == 20 && error == UINT32_MAX);
    CHECK(start_computation(1) == 0);
    reset_records();
    records[0] = star(32, 32, 5, 100000, 1, 0, 0, 1);
    CHECK(select_focus_stars(records, NULL) == -1);
    CHECK(select_focus_stars(records, &exposure) == 0);
    CHECK(add_focus_point(NULL, 0, &exposure) == -1);
    CHECK(add_focus_point(records, 0, NULL) == -1);
    CHECK(finish_focus_point() == -1);
    require_failed_fit();
}

static void selection_ranking_and_duplicates(void)
{
    uint32_t exposure = 1000000;
    int i;
    CHECK(start_computation(1) == 0);
    reset_records();
    records[0] = star(10, 32, 5, 1000, 0.1, 0, 0, 1);  /* bright, poor roundness */
    records[1] = star(50, 32, 5, 300, 1.0, 0, 0, 1);   /* weighted winner */
    records[2] = star(90, 32, 5, 100, 0.5, 0, 0, 1);
    records[3] = star(50.5, 32, 5, 100000, 1.0, 0, 0, 1); /* later duplicate */
    records[0].saturated_pixel_count = 100;
    records[1].saturated_pixel_count = 100;
    records[2].saturated_pixel_count = 100;
    records[3].saturated_pixel_count = 100;
    CHECK(select_focus_stars(records, &exposure) == 0);
    for (i = 0; i < 7; ++i) {
        int32_t position = 4700 + i * 100;
        reset_records();
        records[0] = star(10, 32, hyperbola(position, 5000), 1000, 0.1, 0, 0, 1);
        records[1] = star(50, 32, hyperbola(position, 5100), 300, 1.0, 0, 0, 1);
        records[2] = star(90, 32, hyperbola(position, 5200), 100, 0.5, 0, 0, 1);
        records[3] = star(50.5, 32, 18, 100000, 1.0, 0, 0, 1);
        records[0].saturated_pixel_count = 100;
        records[1].saturated_pixel_count = 100;
        records[2].saturated_pixel_count = 100;
        records[3].saturated_pixel_count = 100;
        CHECK(add_focus_point(records, position, &exposure) == 0);
        CHECK(finish_focus_point() == 0);
    }
    {
        int32_t focus = 0;
        uint32_t error = UINT32_MAX;
        CHECK(finish_computation(&focus, &error) == 0);
        CHECK(focus == 5100);
        CHECK(error == 1);
    }
}

static void selection_only_no_recruitment(void)
{
    uint32_t exposure = 1000000;
    int i;
    CHECK(start_computation(1) == 0);
    reset_records();
    records[0] = star(10, 10, 5, 100000, 1, 0, 0, 1);
    CHECK(select_focus_stars(records, &exposure) == 0);
    for (i = 0; i < 7; ++i) {
        int32_t position = 4700 + i * 100;
        reset_records();
        records[0] = star(100, 100, hyperbola(position, 5000),
            100000, 1, 0, 0, 1);
        CHECK(add_focus_point(records, position, &exposure) == -1);
        CHECK(finish_focus_point() == -1);
    }
    require_failed_fit();
}

static void exposure_guidance_and_bounds(void)
{
    uint32_t exposure;
    int status;
    CHECK(start_computation(1) == 0);
    reset_records();
    records[0] = star(32, 32, 5, 100000, 1, 0, 0.8, 2);
    exposure = 500000;
    CHECK(select_focus_stars(records, &exposure) == 0);
    CHECK(exposure == 1000000); /* baseline uses the detector scale */
    reset_records();
    records[0] = star(32, 32, 5, 100000, 1, 0, 0.8, 2);
    status = add_focus_point(records, 5000, &exposure);
    CHECK(status == 1 && exposure == 2000000);
    CHECK(finish_focus_point() == 0);
    require_failed_fit();

    CHECK(start_computation(1) == 0);
    reset_records();
    records[0] = star(32, 32, 5, 100000, 1, 0, 0, 1);
    exposure = MAX_EXPOSURE_US + 1u;
    CHECK(select_focus_stars(records, &exposure) == 0);
    CHECK(exposure == (uint32_t)MAX_EXPOSURE_US);
    reset_records();
    records[0] = star(32, 32, 5, 100000, 1, 0.9, 0, 0.001);
    CHECK(add_focus_point(records, 5000, &exposure) == 1);
    CHECK(exposure == (uint32_t)MIN_EXPOSURE_US);
    CHECK(finish_focus_point() == 0);
    require_failed_fit();

    CHECK(start_computation(1) == 0);
    reset_records();
    records[0] = star(32, 32, 5, 100000, 1, 0, 0, 0.001);
    exposure = MAX_EXPOSURE_US + 1u;
    CHECK(select_focus_stars(records, &exposure) == 0);
    CHECK(exposure == (uint32_t)MIN_EXPOSURE_US);
    reset_records();
    records[0] = star(32, 32, 5, 100000, 1, 0, 0, 1);
    CHECK(add_focus_point(records, 5000, &exposure) == 0);
    CHECK(exposure == (uint32_t)MIN_EXPOSURE_US);
    CHECK(finish_focus_point() == 0);
    require_failed_fit();

    CHECK(start_computation(1) == 0);
    reset_records();
    records[0] = star(32, 32, 5, 100000, 1, 0, 0, 1);
    exposure = 1000000;
    CHECK(select_focus_stars(records, &exposure) == 0);
    reset_records();
    records[0] = star(32, 32, 5, 100000, 1, 0, 0, 0.5);
    CHECK(add_focus_point(records, 5000, &exposure) == 0);
    CHECK(exposure == 1000000); /* healthy stars do not chase scale changes */
    CHECK(finish_focus_point() == 0);
    require_failed_fit();

    CHECK(start_computation(2) == 0);
    reset_records();
    records[0] = star(32, 32, 5, 100000, 1, 0, 0.9, 2);
    records[1] = star(80, 32, 5, 100000, 1, 0.7, 0, 0.5);
    exposure = 1000000;
    CHECK(select_focus_stars(records, &exposure) == 0);
    reset_records();
    records[0] = star(32, 32, 5, 100000, 1, 0, 0.9, 2);
    records[1] = star(80, 32, 5, 100000, 1, 0.7, 0, 0.5);
    CHECK(add_focus_point(records, 5000, &exposure) == 1);
    CHECK(exposure > 1000000); /* mixed advice follows stronger under-risk coverage */
    CHECK(finish_focus_point() == 0);
    require_failed_fit();
}

static void exposure_top_count_tracks_requested_limit(void)
{
    uint32_t exposure = 1000000;
    int i;
    CHECK(start_computation(33) == 0);
    reset_records();
    for (i = 0; i < 40; ++i)
        records[i] = star(20.0 + 20.0 * i, 32, 5, 100000, 1, 0, 0, 1);
    CHECK(select_focus_stars(records, &exposure) == 0);
    CHECK(exposure == 1000000);
    reset_records();
    for (i = 0; i < 33; ++i)
        records[i] = star(20.0 + 20.0 * i, 32, 5, 100000, 1, 0, 0, 1);
    records[32].under_exposed_risk = 0.9;
    records[32].exposure_scale = 2;
    CHECK(add_focus_point(records, 5000, &exposure) == 1);
    CHECK(exposure == 2000000); /* the 33rd necessary star participates */
    CHECK(finish_focus_point() == 0);
    require_failed_fit();
}

static void best_trial_and_high_risk_samples(void)
{
    uint32_t exposure = 1000000;
    int i;
    CHECK(start_computation(1) == 0);
    reset_records();
    records[0] = star(32, 32, 5, 100000, 1, 0, 0, 1);
    CHECK(select_focus_stars(records, &exposure) == 0);
    for (i = 0; i < 7; ++i) {
        int32_t position = 4700 + i * 100;
        double good = hyperbola(position, 5000);
        reset_records();
        records[0] = star(32, 32, good, 100000, 1, 0.6, 0, 0.5);
        CHECK(add_focus_point(records, position, &exposure) >= 0);
        reset_records();
        records[0] = star(32, 32, 18, 100000, 1, 0.9, 0, 0.5);
        CHECK(add_focus_point(records, position, &exposure) >= 0);
        CHECK(finish_focus_point() == 0);
    }
    {
        int32_t focus = 0;
        uint32_t error = UINT32_MAX;
        CHECK(finish_computation(&focus, &error) == 0);
        CHECK(focus == 5000);
        CHECK(error == 1);
    }

    CHECK(start_computation(1) == 0);
    reset_records();
    records[0] = star(32, 32, 5, 100000, 1, 0, 0, 1);
    CHECK(select_focus_stars(records, &exposure) == 0);
    for (i = 0; i < 7; ++i) {
        int32_t position = 4700 + i * 100;
        add_and_commit(position, hyperbola(position, 5100), 0.8, 0, 0.5, 1);
    }
    {
        int32_t focus = 0;
        uint32_t error = UINT32_MAX;
        CHECK(finish_computation(&focus, &error) == 0);
        CHECK(focus == 5100); /* high-risk finite samples are not filtered */
    }
}

static void missing_and_ambiguous_matches(void)
{
    uint32_t exposure = 1000000;
    int i;
    CHECK(start_computation(1) == 0);
    reset_records();
    records[0] = star(20, 32, 5, 100000, 1, 0, 0, 1);
    records[1] = star(40, 32, 5, 90000, 0.9, 0, 0, 1);
    CHECK(select_focus_stars(records, &exposure) == 0);
    reset_records();
    records[0] = star(30, 32, 5, 100000, 1, 0, 0, 1);
    CHECK(add_focus_point(records, 4600, &exposure) == -1);
    CHECK(finish_focus_point() == -1);
    for (i = 0; i < 7; ++i) {
        int32_t position = 4700 + i * 100;
        reset_records();
        records[0] = star(40, 32, hyperbola(position, 5100),
            90000, 0.9, 0, 0, 1);
        if (i != 2)
            records[1] = star(20, 32, hyperbola(position, 5000),
                100000, 1, 0, 0, 1);
        CHECK(add_focus_point(records, position, &exposure) == 0);
        CHECK(finish_focus_point() == 0);
    }
    {
        int32_t focus = 0;
        uint32_t error = UINT32_MAX;
        CHECK(finish_computation(&focus, &error) == 0);
        CHECK(focus == 5000); /* higher-quality identity remains first */
        CHECK(error == 1);
    }
}

static void empty_point_does_not_commit(void)
{
        uint32_t exposure = 1000000;
        int32_t focus = -17;
        uint32_t error = 0;
        CHECK(start_computation(1) == 0);
        reset_records();
        records[0] = star(32, 32, 5, 100000, 1, 0, 0, 1);
        CHECK(select_focus_stars(records, &exposure) == 0);
        reset_records();
        records[0] = star(100, 100, 5, 100000, 1, 0, 0, 1);
        CHECK(add_focus_point(records, 5100, &exposure) == -1);
        CHECK(finish_focus_point() == -1);
        CHECK(finish_computation(&focus, &error) == -1);
        CHECK(focus == -17 && error == UINT32_MAX);

        CHECK(start_computation(1) == 0);
        reset_records();
        records[0] = star(32, 32, 5, 100000, 1, 0, 0, 1);
        CHECK(select_focus_stars(records, &exposure) == 0);
        reset_records();
        records[0] = star(32, 32, hyperbola(5000, 5000), 100000, 1, 0, 0, 1);
        CHECK(add_focus_point(records, 5000, &exposure) == 0);
        reset_records();
        records[0] = star(100, 100, 5, 100000, 1, 0, 0, 1);
        CHECK(add_focus_point(records, 5000, &exposure) == -1);
        CHECK(finish_focus_point() == 0); /* failed retries leave the best staged trial */
        require_failed_fit();
}

static void focus_spread_uncertainty(void)
{
        uint32_t exposure = 1000000, error = 0;
        int32_t focus = 0;
        int i, j;
        CHECK(start_computation(7) == 0);
        reset_records();
        for (j = 0; j < 7; ++j)
            records[j] = star(20 + 20 * j, 32, 5, 100000, 1, 0, 0, 1);
        CHECK(select_focus_stars(records, &exposure) == 0);
        for (i = 0; i < 9; ++i) {
            int32_t position = 4610 + i * 100;
            reset_records();
            for (j = 0; j < 7; ++j) {
                double centre = 5010 + (j - 3) * 20;
                records[j] = star(20 + 20 * j, 32,
                    hyperbola(position, centre), 100000, 1, 0, 0, 1);
            }
            CHECK(add_focus_point(records, position, &exposure) == 0);
            CHECK(finish_focus_point() == 0);
        }
        CHECK(finish_computation(&focus, &error) == 0);
        CHECK(focus == 5010);
        CHECK(error == 60); /* MAD is 40 steps; 1.4826 * 40 rounds up to 60. */
}

static void direct_hyperbola_and_uncertainty(int32_t centre)
{
    static const double noise[] = {0.05, -0.07, 0.02, 0.08, -0.04, 0.03, -0.02, 0.06, -0.05};
    uint32_t exposure = 1000000, error = 0;
    int32_t focus = 0;
    int i;
    select_one(32, &exposure);
    for (i = 0; i < 9; ++i) {
        int32_t position = centre - 400 + i * 100;
        add_and_commit(position, hyperbola(position, centre) + noise[i],
            0, 0, 1, 0);
    }
    CHECK(finish_computation(&focus, &error) == 0);
    CHECK(llabs((long long)focus - centre) <= 10);
    CHECK(error > 1); /* Jacobian covariance includes the noisy residuals. */
}

static void count_r2_and_repeated_positions(void)
{
    uint32_t exposure = 1000000, error = 0;
    int32_t focus = 0;
    int i;
    select_one(32, &exposure);
    for (i = 0; i < 4; ++i)
        add_and_commit(4800 + i * 100, hyperbola(4800 + i * 100, 4950),
            0, 0, 1, 0);
    require_failed_fit();

    select_one(32, &exposure);
    for (i = 0; i < 5; ++i) {
        int32_t position = 4800 + i * 100;
        add_and_commit(position, hyperbola(position, 5000), 0, 0, 1, 0);
    }
    CHECK(finish_computation(&focus, &error) == 0);
    CHECK(focus == 5000);

    select_one(32, &exposure);
    for (i = 0; i < 5; ++i) {
        int32_t position = 4800 + i * 100;
        add_and_commit(position, hyperbola(position, 5000), 0, 0, 1, 0);
    }
    add_and_commit(5000, 2.0, 0, 0, 1, 0);
    add_and_commit(5000, 2.2, 0, 0, 1, 0);
    add_and_commit(5000, 1.8, 0, 0, 1, 0);
    CHECK(finish_computation(&focus, &error) == 0);
    CHECK(focus == 5000); /* repeated center frames average to one position */

    select_one(32, &exposure);
    for (i = 0; i < 7; ++i)
        add_and_commit(4700 + i * 100, 3.0 + ((i * 7) % 4) * 2.0,
            0, 0, 1, 0);
    require_failed_fit(); /* low R-squared */

    select_one(32, &exposure);
    for (i = 0; i < 7; ++i)
        add_and_commit(4700 + i * 100, 5.0, 0, 0, 1, 0);
    require_failed_fit(); /* flat, unidentifiable curve */
}

static void final_limit_and_empty_errors(void)
{
    uint32_t exposure = 1000000, error = 0;
    int32_t focus = 0;
    int i;
    CHECK(start_computation(1) == 0);
    reset_records();
    records[0] = star(20, 32, 5, 100000, 1.0, 0, 0, 1);
    records[1] = star(60, 32, 5, 100000, 0.9, 0, 0, 1);
    records[2] = star(100, 32, 5, 100000, 0.5, 0, 0, 1);
    CHECK(select_focus_stars(records, &exposure) == 0);
    for (i = 0; i < 7; ++i) {
        int32_t position = 4700 + i * 100;
        reset_records();
        records[0] = star(20, 32, 5, 100000, 1.0, 0, 0, 1);
        records[1] = star(60, 32, hyperbola(position, 5100),
            100000, 0.9, 0, 0, 1);
        records[2] = star(100, 32, hyperbola(position, 5500),
            100000, 0.5, 0, 0, 1);
        CHECK(add_focus_point(records, position, &exposure) == 0);
        CHECK(finish_focus_point() == 0);
    }
    CHECK(finish_computation(&focus, &error) == 0);
    CHECK(focus == 5100); /* first valid fits in frozen quality order, limit one */

    select_one(32, &exposure);
    for (i = 0; i < 7; ++i)
        add_and_commit(4700 + i * 100, hyperbola(4700 + i * 100, 0),
            0, 0, 1, 0);
    require_failed_fit();
}

int main(void)
{
    layout_and_lifecycle();
    selection_ranking_and_duplicates();
    selection_only_no_recruitment();
    exposure_guidance_and_bounds();
    exposure_top_count_tracks_requested_limit();
    best_trial_and_high_risk_samples();
    missing_and_ambiguous_matches();
    empty_point_does_not_commit();
    focus_spread_uncertainty();
    direct_hyperbola_and_uncertainty(2147483000);
    direct_hyperbola_and_uncertainty(-2147483000);
    count_r2_and_repeated_positions();
    final_limit_and_empty_errors();
    puts("Redesigned native focus-curve regression passed");
    return 0;
}
