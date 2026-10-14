/*
 * Public C ABI for native astronomical image analysis.
 *
 * Tutorial for callers:
 *   - Initialize every structure's struct_size field.
 *   - Obtain scratch size with IA_WorkspaceBytes and retain that buffer.
 *   - Analyze a full frame, lock a small set of good stars, then measure the
 *     same stars in each focus frame.
 *   - Store physical quantities in the named milli- or micro-unit fields;
 *     these scaled integers are part of the ABI, not display formatting.
 *
 * Structures are append-only ABI records.  Reordering fields would also
 * invalidate the direct Forth structure declarations in ImageAnalysis.f.
 */

#ifndef IMAGE_ANALYSIS_H
#define IMAGE_ANALYSIS_H

#include <stdint.h>

#ifdef _WIN32
#define IA_EXPORT __declspec(dllexport)
#define IA_CALL __cdecl
#else
#define IA_EXPORT
#define IA_CALL
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define IA_API_VERSION 2u
#define IA_MILLI 1000
#define IA_FOCUS_MAX_STARS 7u
#define IA_FOCUS_MAX_POINTS 15u
#ifndef IA_DEBUG_LOG_PATH
#define IA_DEBUG_LOG_PATH "E:\\images\\logs\\image-analysis-native.log"
#endif

enum IAStatus {
    IA_OK = 0,
    IA_INVALID_ARGUMENT = -1,
    IA_WORKSPACE_TOO_SMALL = -2,
    IA_NO_STARS = -3
};

enum IAStarFlags {
    IA_STAR_SATURATED = 1u << 0,
    IA_STAR_EDGE = 1u << 1,
    IA_STAR_BLENDED = 1u << 2,
    IA_STAR_TOO_SMALL = 1u << 3,
    IA_STAR_LOW_SNR = 1u << 4
};

enum IAFocusAction {
    IA_FOCUS_CAPTURE = 1,
    IA_FOCUS_RETAKE = 2,
    IA_FOCUS_COMPLETE = 3,
    IA_FOCUS_FAILED = 4
};

enum IAFocusReason {
    IA_FOCUS_REASON_NONE = 0,
    IA_FOCUS_REASON_SATURATED = 1u << 0,
    IA_FOCUS_REASON_TOO_FAINT = 1u << 1,
    IA_FOCUS_REASON_NO_STARS = 1u << 2,
    IA_FOCUS_REASON_NO_CURVE = 1u << 3
};

typedef struct IAConfig {
    uint32_t struct_size;
    uint32_t detection_sigma_milli;
    uint32_t minimum_threshold_adu;
    uint32_t saturation_adu;
    uint32_t minimum_area_pixels;
    uint32_t minimum_separation_pixels;
    uint32_t aperture_radius_pixels;
    uint32_t annulus_inner_radius_pixels;
    uint32_t annulus_outer_radius_pixels;
    uint32_t minimum_snr_milli;
    uint32_t blend_minimum_contrast_milli;
} IAConfig;

typedef struct IAStar {
    uint32_t struct_size;
    int32_t x_milli;
    int32_t y_milli;
    uint32_t flux_adu;
    uint32_t peak_adu;
    uint32_t saturated_pixels;
    int32_t background_milli_adu;
    int32_t noise_milli_adu;
    int32_t snr_milli;
    int32_t hfr_milli_pixels;
    int32_t hfd_milli_pixels;
    int32_t fwhm_milli_pixels;
    int32_t sigma_major_milli_pixels;
    int32_t sigma_minor_milli_pixels;
    int32_t ellipticity_milli;
    int32_t theta_millidegrees;
    uint32_t positive_pixels;
    int32_t left;
    int32_t top;
    int32_t right;
    int32_t bottom;
    uint32_t flags;
} IAStar;

typedef struct IAFrameSummary {
    uint32_t struct_size;
    uint32_t detected_stars;
    uint32_t measured_stars;
    uint32_t included_stars;
    uint32_t saturated_stars;
    uint32_t edge_stars;
    uint32_t blended_stars;
    uint32_t low_snr_stars;
    int32_t background_milli_adu;
    int32_t noise_milli_adu;
    uint32_t detection_threshold_adu;
    int32_t median_hfd_milli_pixels;
    int32_t mad_hfd_milli_pixels;
    int32_t median_fwhm_milli_pixels;
    int32_t mad_fwhm_milli_pixels;
    int32_t median_ellipticity_milli;
    int32_t median_snr_milli;
} IAFrameSummary;

typedef struct IARect {
    int32_t x;
    int32_t y;
    int32_t width;
    int32_t height;
} IARect;

typedef struct IAFocusConfig {
    uint32_t struct_size;
    uint32_t star_limit;
    uint32_t point_count;
    uint32_t focus_increment;
    uint32_t target_minimum_peak_adu;
    uint32_t target_maximum_peak_adu;
    uint32_t minimum_exposure_ms;
    uint32_t maximum_exposure_ms;
    uint32_t initial_maximum_change_milli;
    uint32_t step_maximum_change_milli;
    uint32_t maximum_retries;
    uint32_t roi_margin_pixels;
    uint32_t search_radius_pixels;
    uint32_t reject_star_flags;
    uint32_t minimum_branch_points;
    uint32_t maximum_residual_milli;
} IAFocusConfig;

typedef struct IAFocusRequest {
    uint32_t struct_size;
    uint32_t action;
    uint32_t reason_flags;
    int32_t focus_position;
    uint32_t exposure_milliseconds;
    IARect roi;
} IAFocusRequest;

typedef struct IAFocusResult {
    uint32_t struct_size;
    uint32_t stars_locked;
    uint32_t stars_fitted;
    int32_t focus_milli_steps;
    int32_t mad_focus_milli_steps;
    int32_t median_rms_residual_milli;
    uint32_t recommended_exposure_milliseconds;
} IAFocusResult;

typedef struct IAFocusDiagnostic {
    uint32_t struct_size;
    uint32_t stars_measured;
    uint32_t saturated_stars;
    int32_t median_hfd_milli_pixels;
    int32_t mad_hfd_milli_pixels;
} IAFocusDiagnostic;

IA_EXPORT uint32_t IA_CALL IA_Version(void);
IA_EXPORT void IA_CALL IA_SetDebug(uint32_t enabled);
IA_EXPORT void IA_CALL IA_DefaultConfig(IAConfig *config);
IA_EXPORT uint32_t IA_CALL IA_WorkspaceBytes(uint32_t maximum_stars);

IA_EXPORT int IA_CALL IA_AnalyzeFrame(
    const uint16_t *pixels,
    uint32_t frame_size_bytes,
    uint32_t width,
    uint32_t height,
    uint32_t stride_pixels,
    const IAConfig *config,
    void *workspace,
    uint32_t workspace_bytes,
    IAStar *stars,
    uint32_t star_capacity,
    IAFrameSummary *summary);

IA_EXPORT int IA_CALL IA_SummarizeStars(
    const IAStar *stars,
    uint32_t star_count,
    uint32_t reject_flags,
    void *workspace,
    uint32_t workspace_bytes,
    IAFrameSummary *summary);

IA_EXPORT int IA_CALL IA_MeasureLockedStars(
    const uint16_t *pixels,
    uint32_t width,
    uint32_t height,
    uint32_t stride_pixels,
    const IAConfig *config,
    const IAStar *locked_stars,
    uint32_t star_count,
    uint32_t search_radius_pixels,
    IAStar *measurements);

IA_EXPORT int IA_CALL IA_RecommendExposure(
    uint32_t current_milliseconds,
    uint32_t measured_peak_adu,
    uint32_t background_adu,
    uint32_t target_peak_adu,
    uint32_t minimum_milliseconds,
    uint32_t maximum_milliseconds,
    uint32_t maximum_change_milli,
    uint32_t *recommended_milliseconds);

IA_EXPORT int IA_CALL IA_ComputeStarROI(
    const IAStar *stars,
    uint32_t star_count,
    uint32_t image_width,
    uint32_t image_height,
    uint32_t margin_pixels,
    IARect *roi);

IA_EXPORT uint32_t IA_CALL IA_FocusRunBytes(void);
IA_EXPORT void IA_CALL IA_DefaultFocusConfig(IAFocusConfig *config);

IA_EXPORT int IA_CALL IA_BeginFocusRun(
    void *run,
    uint32_t run_bytes,
    const IAFocusConfig *focus_config,
    const IAConfig *analysis_config,
    const uint16_t *pixels,
    uint32_t frame_size_bytes,
    uint32_t width,
    uint32_t height,
    uint32_t stride_pixels,
    int32_t center_focus_position,
    uint32_t exposure_milliseconds,
    IAFocusRequest *request);

IA_EXPORT int IA_CALL IA_SubmitFocusFrame(
    void *run,
    const uint16_t *pixels,
    uint32_t width,
    uint32_t height,
    uint32_t stride_pixels,
    int32_t origin_x,
    int32_t origin_y,
    int32_t focus_position,
    uint32_t exposure_milliseconds,
    IAFocusRequest *request);

IA_EXPORT int IA_CALL IA_EndFocusRun(
    void *run,
    IAFocusResult *result);

IA_EXPORT int IA_CALL IA_DiagnoseFocusFrame(
    void *run,
    const uint16_t *pixels,
    uint32_t width,
    uint32_t height,
    uint32_t stride_pixels,
    int32_t origin_x,
    int32_t origin_y,
    IAFocusDiagnostic *diagnostic);

#ifdef __cplusplus
}
#endif

#endif
