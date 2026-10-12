/*
 * Executable tutorial and regression suite for the native API.
 *
 * Each test constructs a small synthetic image or focus series with a known
 * answer, invokes the same exported functions used by Forth, and checks
 * physically meaningful results rather than implementation details.  Start
 * with test_detection_and_measurement: it demonstrates workspace allocation,
 * default configuration, frame analysis, star locking, and remeasurement.
 * Later tests isolate ROI, exposure, and focus-model edge cases.
 */

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

static long debug_log_size(void)
{
    FILE *file = fopen(IA_DEBUG_LOG_PATH, "rb");
    long size;
    if (file == NULL) {
        return -1;
    }
    fseek(file, 0, SEEK_END);
    size = ftell(file);
    fclose(file);
    return size;
}

static void test_debug_logging(void)
{
    long before;
    long after;

    IA_SetDebug(0u);
    before = debug_log_size();
    IA_Version();
    CHECK(debug_log_size() == before, "native debug disabled by default");

    IA_SetDebug(1u);
    IA_Version();
    IA_SetDebug(0u);
    after = debug_log_size();
    CHECK(after > before, "native debug appends and flushes");
}

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

static void test_locked_noise_rejection(void)
{
    uint16_t image[WIDTH * HEIGHT];
    IAConfig config;
    IAStar locked;
    IAStar measured;
    uint32_t i;
    int status;

    /* Maintenance: rectifying aperture noise gives a large positive flux
       and a plausible HFD even after a locked star has disappeared. */
    for (i = 0u; i < WIDTH * HEIGHT; ++i) {
        image[i] = (uint16_t)(1000u + (i * 17u + i / WIDTH * 13u) % 41u);
    }
    IA_DefaultConfig(&config);
    memset(&locked, 0, sizeof(locked));
    locked.struct_size = sizeof(locked);
    locked.x_milli = 80000;
    locked.y_milli = 60000;
    status = IA_MeasureLockedStars(
        image, WIDTH, HEIGHT, WIDTH, &config, &locked, 1u, 16u, &measured);
    CHECK(status == IA_OK, "measure missing locked star");
    CHECK((measured.flags & IA_STAR_LOW_SNR) != 0u,
        "noise-only aperture must not masquerade as a usable locked star");
    CHECK(measured.snr_milli < (int32_t)config.minimum_snr_milli,
        "noise-only signal quality uses signed background-subtracted flux");
}

static void test_focus_rejected_centroid(void)
{
    uint16_t image[WIDTH * HEIGHT];
    uint32_t run_bytes = IA_FocusRunBytes();
    void *run = malloc(run_bytes);
    IAConfig analysis_config;
    IAFocusConfig focus_config;
    IAFocusRequest request;
    IAFocusDiagnostic diagnostic;
    uint32_t i;
    uint32_t attempt;
    int status;

    CHECK(run != NULL, "allocate rejected-centroid test run");
    if (run == NULL) return;
    IA_DefaultConfig(&analysis_config);
    IA_DefaultFocusConfig(&focus_config);
    focus_config.star_limit = 2u;
    request.struct_size = sizeof(request);
    for (i = 0u; i < WIDTH * HEIGHT; ++i) image[i] = 1000u;
    add_gaussian(image, 50.0, 40.0, 2.0, 2.0, 30000.0);
    add_gaussian(image, 130.0, 80.0, 2.0, 2.0, 25000.0);
    status = IA_BeginFocusRun(
        run, run_bytes, &focus_config, &analysis_config,
        image, WIDTH, HEIGHT, WIDTH, 5000, 1000u, &request);
    CHECK(status == IA_OK, "begin rejected-centroid test");
    if (status != IA_OK) {
        free(run);
        return;
    }
    /* Maintenance: rejected peaks must not walk a locked search center away
       from its last valid star while another usable star advances the run. */
    for (attempt = 0u; attempt < 3u; ++attempt) {
        for (i = 0u; i < WIDTH * HEIGHT; ++i) image[i] = 1000u;
        add_gaussian(image, 130.0, 80.0, 2.0, 2.0, 25000.0);
        image[40u * WIDTH + 64u + attempt * 14u] = 65535u;
        status = IA_SubmitFocusFrame(
            run, image, WIDTH, HEIGHT, WIDTH, 0, 0,
            request.focus_position, 1000u, &request);
        CHECK(status == IA_OK && request.action == IA_FOCUS_CAPTURE,
            "another usable star advances past a rejected centroid");
    }
    for (i = 0u; i < WIDTH * HEIGHT; ++i) image[i] = 1000u;
    add_gaussian(image, 50.0, 40.0, 2.0, 2.0, 30000.0);
    add_gaussian(image, 130.0, 80.0, 2.0, 2.0, 25000.0);
    diagnostic.struct_size = sizeof(diagnostic);
    status = IA_DiagnoseFocusFrame(
        run, image, WIDTH, HEIGHT, WIDTH, 0, 0, &diagnostic);
    CHECK(status == IA_OK && diagnostic.stars_measured == 2u,
        "returning star is recovered at its last valid locked position");
    free(run);
}

static void make_focus_frame(
    uint16_t *image,
    int32_t focus,
    uint32_t exposure_milliseconds)
{
    double sigma = 1.5 + fabs(focus - 5030) * 0.02;
    double exposure_scale = exposure_milliseconds / 1000.0;
    double focus_scale = 2.25 / (sigma * sigma);
    uint32_t i;

    for (i = 0u; i < WIDTH * HEIGHT; ++i) {
        image[i] = (uint16_t)(1000u + (i * 17u) % 5u);
    }
    add_gaussian(
        image, 50.0, 40.0, sigma, sigma,
        30000.0 * exposure_scale * focus_scale);
    add_gaussian(
        image, 75.0, 55.0, sigma * 1.05, sigma * 1.05,
        26000.0 * exposure_scale * focus_scale);
    add_gaussian(
        image, 100.0, 70.0, sigma * 0.95, sigma * 0.95,
        22000.0 * exposure_scale * focus_scale);
}

static void test_focus_run(void)
{
    uint16_t *image =
        (uint16_t *)malloc(WIDTH * HEIGHT * sizeof(uint16_t));
    uint32_t run_bytes = IA_FocusRunBytes();
    void *run = malloc(run_bytes);
    IAConfig analysis_config;
    IAFocusConfig focus_config;
    IAFocusRequest request;
    IAFocusResult result;
    IAFocusDiagnostic diagnostic;
    uint32_t observations = 0u;
    int status;

    CHECK(image != NULL && run != NULL, "allocate focus-run buffers");
    if (image == NULL || run == NULL) {
        free(image);
        free(run);
        return;
    }
    IA_DefaultConfig(&analysis_config);
    IA_DefaultFocusConfig(&focus_config);
    focus_config.star_limit = 3u;
    focus_config.target_minimum_peak_adu = 5000u;
    focus_config.target_maximum_peak_adu = 60000u;
    focus_config.maximum_residual_milli = 1500u;
    make_focus_frame(image, 5030, 1000u);
    request.struct_size = sizeof(request);
    status = IA_BeginFocusRun(
        run, run_bytes, &focus_config, &analysis_config,
        image, WIDTH, HEIGHT, WIDTH, 5030, 1000u, &request);
    CHECK(status == IA_OK, "begin bitmap-driven focus run");
    CHECK(request.action == IA_FOCUS_CAPTURE, "request first focus frame");
    CHECK(request.roi.width > 0 && request.roi.height > 0,
        "recommend compact focus ROI");

    while (status == IA_OK &&
           (request.action == IA_FOCUS_CAPTURE ||
            request.action == IA_FOCUS_RETAKE) &&
           observations < 30u) {
        make_focus_frame(
            image, request.focus_position,
            request.exposure_milliseconds);
        status = IA_SubmitFocusFrame(
            run, image, WIDTH, HEIGHT, WIDTH, 0, 0,
            request.focus_position, request.exposure_milliseconds,
            &request);
        ++observations;
    }
    CHECK(status == IA_OK, "submit focus frames");
    CHECK(request.action == IA_FOCUS_COMPLETE, "complete focus run");
    result.struct_size = sizeof(result);
    status = IA_EndFocusRun(run, &result);
    CHECK(status == IA_OK, "return focus result");
    CHECK(result.stars_locked == 3u && result.stars_fitted >= 1u,
        "fit independently locked stars");
    CHECK(abs(result.focus_milli_steps - 5030000) < 5000,
        "intersecting lines recover focus");
    CHECK(result.recommended_exposure_milliseconds > 0u,
        "recommend verification exposure");

    make_focus_frame(
        image, (result.focus_milli_steps + 500) / IA_MILLI,
        result.recommended_exposure_milliseconds);
    diagnostic.struct_size = sizeof(diagnostic);
    status = IA_DiagnoseFocusFrame(
        run, image, WIDTH, HEIGHT, WIDTH, 0, 0, &diagnostic);
    CHECK(status == IA_OK, "diagnose verification frame");
    CHECK(diagnostic.stars_measured >= 1u &&
          diagnostic.median_hfd_milli_pixels > 0,
        "return HFD observing diagnostic");

    image[40u * WIDTH + 50u] = UINT16_MAX;
    diagnostic.struct_size = sizeof(diagnostic);
    status = IA_DiagnoseFocusFrame(
        run, image, WIDTH, HEIGHT, WIDTH, 0, 0, &diagnostic);
    CHECK(status == IA_OK, "diagnose partially saturated frame");
    CHECK(diagnostic.saturated_stars == 1u &&
          diagnostic.stars_measured >= 1u,
        "retain unsaturated diagnostic stars");
    free(run);
    free(image);
}

static void test_focus_diagnostic_exposure(void)
{
    uint16_t *image =
        (uint16_t *)malloc(WIDTH * HEIGHT * sizeof(uint16_t));
    uint32_t run_bytes = IA_FocusRunBytes();
    void *run = malloc(run_bytes);
    IAConfig analysis_config;
    IAFocusConfig focus_config;
    IAFocusRequest request;
    IAFocusResult result;
    IAFocusDiagnostic diagnostic;
    uint32_t variant;

    CHECK(image != NULL && run != NULL, "allocate diagnostic exposure buffers");
    if (image == NULL || run == NULL) {
        free(image);
        free(run);
        return;
    }
    IA_DefaultConfig(&analysis_config);
    IA_DefaultFocusConfig(&focus_config);
    focus_config.star_limit = 3u;
    focus_config.target_minimum_peak_adu = 5000u;
    focus_config.target_maximum_peak_adu = 60000u;
    focus_config.maximum_residual_milli = 1500u;
    focus_config.search_radius_pixels = 4u;

    for (variant = 0u; variant < 3u; ++variant) {
        uint32_t observations = 0u;
        int status;
        make_focus_frame(image, 5030, 997u);
        request.struct_size = sizeof(request);
        status = IA_BeginFocusRun(
            run, run_bytes, &focus_config, &analysis_config,
            image, WIDTH, HEIGHT, WIDTH, 5030, 997u, &request);
        while (status == IA_OK &&
               (request.action == IA_FOCUS_CAPTURE ||
                request.action == IA_FOCUS_RETAKE) &&
               observations < 30u) {
            int32_t position = request.focus_position;
            uint32_t actual_exposure = position == 5030 ? 1234u : 1500u;
            make_focus_frame(image, position, actual_exposure);
            if (position == 5030 && variant == 1u) {
                image[40u * WIDTH + 50u] = UINT16_MAX;
            } else if (position == 5030 && variant == 2u) {
                uint32_t x;
                uint32_t y;
                for (y = 0u; y < HEIGHT; ++y) {
                    for (x = 0u; x < 65u; ++x) {
                        image[y * WIDTH + x] = 1000u;
                    }
                }
            }
            status = IA_SubmitFocusFrame(
                run, image, WIDTH, HEIGHT, WIDTH, 0, 0,
                position, actual_exposure, &request);
            if (position == 5030 && variant != 0u) {
                uint32_t expected_reason = variant == 1u ?
                    IA_FOCUS_REASON_SATURATED : IA_FOCUS_REASON_TOO_FAINT;
                CHECK((request.reason_flags & expected_reason) != 0u,
                    "mark unusable nearest sweep exposure");
            }
            ++observations;
        }
        CHECK(status == IA_OK && request.action == IA_FOCUS_COMPLETE,
            "complete diagnostic exposure sweep");
        result.struct_size = sizeof(result);
        status = IA_EndFocusRun(run, &result);
        CHECK(status == IA_OK &&
              abs(result.focus_milli_steps - 5030000) < 5000,
            "fit diagnostic exposure sweep");
        CHECK(result.recommended_exposure_milliseconds ==
              (variant == 0u ? 1234u : 997u),
            "reuse exact nearest successful exposure, latest on ties");
        CHECK(request.exposure_milliseconds ==
              result.recommended_exposure_milliseconds,
            "completion request carries diagnostic exposure");
        make_focus_frame(
            image, (result.focus_milli_steps + 500) / IA_MILLI,
            result.recommended_exposure_milliseconds);
        diagnostic.struct_size = sizeof(diagnostic);
        status = IA_DiagnoseFocusFrame(
            run, image, WIDTH, HEIGHT, WIDTH, 0, 0, &diagnostic);
        CHECK(status == IA_OK && diagnostic.stars_measured == 3u &&
              diagnostic.saturated_stars == 0u,
            "nearest successful exposure gives unsaturated diagnostic");
    }
    free(run);
    free(image);
}

static void test_focus_failures(void)
{
    uint16_t *image =
        (uint16_t *)malloc(WIDTH * HEIGHT * sizeof(uint16_t));
    uint32_t run_bytes = IA_FocusRunBytes();
    void *run = malloc(run_bytes);
    IAConfig analysis_config;
    IAFocusConfig focus_config;
    IAFocusRequest request;
    uint32_t submissions;
    int saw_retake;
    int status;

    CHECK(image != NULL && run != NULL, "allocate focus failure buffers");
    if (image == NULL || run == NULL) {
        free(image);
        free(run);
        return;
    }
    IA_DefaultConfig(&analysis_config);
    IA_DefaultFocusConfig(&focus_config);
    focus_config.star_limit = 3u;
    focus_config.target_minimum_peak_adu = 5000u;
    focus_config.target_maximum_peak_adu = 60000u;

    make_focus_frame(image, 5030, 1000u);
    request.struct_size = sizeof(request);
    status = IA_BeginFocusRun(
        run, run_bytes, &focus_config, &analysis_config,
        image, WIDTH, HEIGHT, WIDTH, 5030, 1000u, &request);
    saw_retake = 0;
    for (submissions = 0u;
         status == IA_OK && request.action != IA_FOCUS_FAILED &&
         submissions < 30u;
         ++submissions) {
        make_focus_frame(image, request.focus_position, 10000u);
        status = IA_SubmitFocusFrame(
            run, image, WIDTH, HEIGHT, WIDTH, 0, 0,
            request.focus_position, 10000u,
            &request);
        if (request.action == IA_FOCUS_RETAKE) saw_retake = 1;
    }
    CHECK(status == IA_OK && saw_retake, "retake saturated focus frame");
    CHECK(request.action == IA_FOCUS_FAILED &&
          (request.reason_flags &
           (IA_FOCUS_REASON_SATURATED | IA_FOCUS_REASON_NO_STARS)) ==
          (IA_FOCUS_REASON_SATURATED | IA_FOCUS_REASON_NO_STARS),
        "fail after saturated-frame retries");

    make_focus_frame(image, 5030, 1000u);
    request.struct_size = sizeof(request);
    status = IA_BeginFocusRun(
        run, run_bytes, &focus_config, &analysis_config,
        image, WIDTH, HEIGHT, WIDTH, 5030, 1000u, &request);
    for (submissions = 0u;
         status == IA_OK &&
         (request.action == IA_FOCUS_CAPTURE ||
          request.action == IA_FOCUS_RETAKE) &&
         submissions < 30u;
         ++submissions) {
        make_focus_frame(image, 4900, request.exposure_milliseconds);
        status = IA_SubmitFocusFrame(
            run, image, WIDTH, HEIGHT, WIDTH, 0, 0,
            request.focus_position, request.exposure_milliseconds,
            &request);
    }
    CHECK(status == IA_OK && request.action == IA_FOCUS_FAILED &&
          (request.reason_flags & IA_FOCUS_REASON_NO_CURVE) != 0u,
        "reject focus series without opposing branches");
    free(run);
    free(image);
}

int main(void)
{
    CHECK(IA_Version() == IA_API_VERSION, "API version");
    test_debug_logging();
    test_detection_and_measurement();
    test_calculations();
    test_locked_noise_rejection();
    test_focus_rejected_centroid();
    test_focus_run();
    test_focus_diagnostic_exposure();
    test_focus_failures();
    if (failures != 0) {
        fprintf(stderr, "%d native image-analysis test(s) failed\n", failures);
        return 1;
    }
    puts("Native image-analysis tests passed");
    return 0;
}
