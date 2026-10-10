#include "image_analysis.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define WIDTH 160u
#define HEIGHT 120u
#define STAR_CAPACITY 16u

static int failures;

#define CHECK(condition, message) \
    do { \
        if (!(condition)) { \
            fprintf(stderr, "FAIL: %s\n", message); \
            ++failures; \
        } \
    } while (0)

static void add_gaussian(
    uint16_t *image,
    double center_x,
    double center_y,
    double sigma_x,
    double sigma_y,
    double amplitude)
{
    uint32_t y;
    uint32_t x;
    for (y = 0u; y < HEIGHT; ++y) {
        for (x = 0u; x < WIDTH; ++x) {
            double dx = (x - center_x) / sigma_x;
            double dy = (y - center_y) / sigma_y;
            double value = image[y * WIDTH + x] +
                amplitude * exp(-0.5 * (dx * dx + dy * dy));
            image[y * WIDTH + x] = value >= 65535.0 ? 65535u : (uint16_t)(value + 0.5);
        }
    }
}

static void test_detection_and_measurement(void)
{
    uint16_t *image = (uint16_t *)malloc(WIDTH * HEIGHT * sizeof(uint16_t));
    IAConfig config;
    IAStar stars[STAR_CAPACITY];
    IAFrameSummary summary;
    uint32_t workspace_bytes = IA_WorkspaceBytes(STAR_CAPACITY);
    void *workspace = malloc(workspace_bytes);
    int status;
    uint32_t i;
    IAStar *center_star = NULL;
    IAStar locked;
    IAStar measured;

    CHECK(image != NULL && workspace != NULL, "allocate synthetic test buffers");
    if (image == NULL || workspace == NULL) {
        free(image);
        free(workspace);
        return;
    }

    for (i = 0u; i < WIDTH * HEIGHT; ++i) {
        image[i] = (uint16_t)(1000u + ((i * 17u + i / WIDTH * 13u) % 7u));
    }
    add_gaussian(image, 70.25, 50.75, 2.0, 2.0, 30000.0);
    add_gaussian(image, 120.5, 80.5, 3.0, 2.0, 18000.0);
    add_gaussian(image, 25.0, 25.0, 1.5, 1.5, 70000.0);

    IA_DefaultConfig(&config);
    summary.struct_size = sizeof(summary);
    status = IA_AnalyzeFrame(
        image, WIDTH, HEIGHT, WIDTH, &config, workspace, workspace_bytes,
        stars, STAR_CAPACITY, &summary);

    CHECK(status == IA_OK, "analyze synthetic image");
    CHECK(summary.measured_stars >= 3u, "detect three synthetic stars");
    CHECK(summary.background_milli_adu >= 1000000 &&
          summary.background_milli_adu <= 1006000, "global background");
    CHECK(summary.saturated_stars >= 1u, "report saturated star");

    for (i = 0u; i < summary.measured_stars; ++i) {
        if (abs(stars[i].x_milli - 70250) < 1000 &&
            abs(stars[i].y_milli - 50750) < 1000) {
            center_star = &stars[i];
            break;
        }
    }
    CHECK(center_star != NULL, "find central synthetic star");
    if (center_star != NULL) {
        CHECK(abs(center_star->x_milli - 70250) < 150, "centroid x");
        CHECK(abs(center_star->y_milli - 50750) < 150, "centroid y");
        CHECK(center_star->hfd_milli_pixels > 4300 &&
              center_star->hfd_milli_pixels < 5200, "Gaussian HFD");
        CHECK(center_star->fwhm_milli_pixels > 4300 &&
              center_star->fwhm_milli_pixels < 5200, "Gaussian-equivalent FWHM");
        CHECK(center_star->ellipticity_milli < 50, "round-star ellipticity");
        locked = *center_star;
        locked.x_milli += 2000;
        locked.y_milli -= 2000;
        status = IA_MeasureLockedStars(
            image, WIDTH, HEIGHT, WIDTH, &config, &locked, 1u, 4u, &measured);
        CHECK(status == IA_OK, "remeasure locked star");
        CHECK(abs(measured.x_milli - 70250) < 150, "relocate locked star x");
        CHECK(abs(measured.y_milli - 50750) < 150, "relocate locked star y");
    }

    free(workspace);
    free(image);
}

static void test_calculations(void)
{
    uint32_t exposure = 0u;
    IAStar stars[2];
    IARect roi;
    int status;

    status = IA_RecommendExposure(1000u, 21000u, 1000u, 41000u,
        100u, 10000u, 4000u, &exposure);
    CHECK(status == IA_OK, "recommend exposure");
    CHECK(exposure == 2000u, "linear exposure recommendation");

    memset(stars, 0, sizeof(stars));
    stars[0].x_milli = 20000;
    stars[0].y_milli = 30000;
    stars[1].x_milli = 80000;
    stars[1].y_milli = 60000;
    status = IA_ComputeStarROI(stars, 2u, 100u, 100u, 10u, &roi);
    CHECK(status == IA_OK, "compute star ROI");
    CHECK(roi.x == 10 && roi.y == 20 && roi.width == 81 && roi.height == 51,
        "ideal ROI bounds");
}

static int32_t piecewise_metric(
    int32_t focus,
    int32_t center,
    int32_t baseline,
    int32_t left_slope,
    int32_t right_slope)
{
    int32_t distance = focus - center;
    int32_t slope = distance < 0 ? left_slope : right_slope;
    if (distance < 0) distance = -distance;
    return baseline + (int32_t)((int64_t)distance * slope / IA_MILLI);
}

static int32_t hyperbola_metric(
    int32_t focus,
    int32_t center,
    int32_t baseline,
    int32_t radius,
    int32_t slope)
{
    double distance = focus - center;
    double scaled = distance * slope / IA_MILLI;
    return baseline + (int32_t)(sqrt((double)radius * radius + scaled * scaled) + 0.5);
}

static void test_focus_fitting(void)
{
    IAFocusSample samples[11];
    IAFocusModel model;
    IAFocusFit fits[3];
    IAFocusSummary summary;
    int32_t workspace[3];
    uint32_t i;
    int status;

    memset(&model, 0, sizeof(model));
    model.struct_size = sizeof(model);
    model.kind = IA_FOCUS_PIECEWISE_V;
    model.baseline_milli = 3000;
    model.left_slope_micro_per_step = 38000;
    model.right_slope_micro_per_step = 42000;
    for (i = 0u; i < 11u; ++i) {
        samples[i].focus_position = 4900 + (int32_t)i * 20;
        samples[i].metric_milli = piecewise_metric(
            samples[i].focus_position, 5037, 3250, 38000, 42000);
        samples[i].weight_milli = IA_MILLI;
        samples[i].flags = 0u;
    }
    fits[0].struct_size = sizeof(fits[0]);
    status = IA_FitFocusModel(&model, samples, 11u, 0u, 250u, &fits[0]);
    CHECK(status == IA_OK, "fit piecewise focus model");
    CHECK(abs(fits[0].focus_milli_steps - 5037000) <= 250,
        "piecewise focus position");
    CHECK(fits[0].rms_residual_milli <= 10, "piecewise residual");
    CHECK(fits[0].points_left > 0u && fits[0].points_right > 0u,
        "piecewise bracket geometry");

    model.kind = IA_FOCUS_HYPERBOLA;
    model.baseline_milli = 500;
    model.hyperbola_radius_milli = 2200;
    model.hyperbola_slope_micro_per_step = 45000;
    for (i = 0u; i < 11u; ++i) {
        samples[i].metric_milli = hyperbola_metric(
            samples[i].focus_position, 5026, 700,
            model.hyperbola_radius_milli,
            model.hyperbola_slope_micro_per_step);
    }
    fits[1].struct_size = sizeof(fits[1]);
    status = IA_FitFocusModel(&model, samples, 11u, 0u, 250u, &fits[1]);
    CHECK(status == IA_OK, "fit hyperbolic focus model");
    CHECK(abs(fits[1].focus_milli_steps - 5026000) <= 250,
        "hyperbolic focus position");
    CHECK(fits[1].rms_residual_milli <= 10, "hyperbolic residual");

    fits[2] = fits[1];
    fits[0].focus_milli_steps = 5025000;
    fits[1].focus_milli_steps = 5026000;
    fits[2].focus_milli_steps = 5099000;
    fits[2].flags = IA_FOCUS_AT_BOUNDARY;
    summary.struct_size = sizeof(summary);
    status = IA_CombineFocusFits(
        fits, 3u, IA_FOCUS_AT_BOUNDARY, workspace, sizeof(workspace), &summary);
    CHECK(status == IA_OK, "combine focus fits");
    CHECK(summary.fits_included == 2u, "exclude flagged focus fit");
    CHECK(summary.median_focus_milli_steps == 5025500, "combined focus median");
    CHECK(summary.mad_focus_milli_steps == 500, "combined focus MAD");
}

int main(void)
{
    CHECK(IA_Version() == IA_API_VERSION, "API version");
    test_detection_and_measurement();
    test_calculations();
    test_focus_fitting();
    if (failures != 0) {
        fprintf(stderr, "%d native image-analysis test(s) failed\n", failures);
        return 1;
    }
    puts("Native image-analysis tests passed");
    return 0;
}
