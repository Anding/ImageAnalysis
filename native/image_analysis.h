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

#define IA_API_VERSION 1u
#define IA_MILLI 1000

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

enum IAFocusModelKind {
    IA_FOCUS_PIECEWISE_V = 1,
    IA_FOCUS_HYPERBOLA = 2
};

enum IAFocusFitFlags {
    IA_FOCUS_INSUFFICIENT_POINTS = 1u << 0,
    IA_FOCUS_ONE_SIDED = 1u << 1,
    IA_FOCUS_AT_BOUNDARY = 1u << 2,
    IA_FOCUS_INVALID_MODEL = 1u << 3
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

typedef struct IAFocusModel {
    uint32_t struct_size;
    uint32_t kind;
    int32_t baseline_milli;
    int32_t left_slope_micro_per_step;
    int32_t right_slope_micro_per_step;
    int32_t hyperbola_radius_milli;
    int32_t hyperbola_slope_micro_per_step;
} IAFocusModel;

typedef struct IAFocusSample {
    int32_t focus_position;
    int32_t metric_milli;
    uint32_t weight_milli;
    uint32_t flags;
} IAFocusSample;

typedef struct IAFocusFit {
    uint32_t struct_size;
    int32_t focus_milli_steps;
    int32_t vertical_offset_milli;
    int32_t rms_residual_milli;
    uint32_t points_used;
    uint32_t points_left;
    uint32_t points_right;
    uint32_t flags;
} IAFocusFit;

typedef struct IAFocusSummary {
    uint32_t struct_size;
    uint32_t fits_included;
    int32_t median_focus_milli_steps;
    int32_t mad_focus_milli_steps;
    int32_t median_rms_residual_milli;
} IAFocusSummary;

IA_EXPORT uint32_t IA_CALL IA_Version(void);
IA_EXPORT void IA_CALL IA_DefaultConfig(IAConfig *config);
IA_EXPORT uint32_t IA_CALL IA_WorkspaceBytes(uint32_t maximum_stars);

IA_EXPORT int IA_CALL IA_AnalyzeFrame(
    const uint16_t *pixels,
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

IA_EXPORT int IA_CALL IA_FitFocusModel(
    const IAFocusModel *model,
    const IAFocusSample *samples,
    uint32_t sample_count,
    uint32_t reject_sample_flags,
    uint32_t resolution_milli_steps,
    IAFocusFit *fit);

IA_EXPORT int IA_CALL IA_CombineFocusFits(
    const IAFocusFit *fits,
    uint32_t fit_count,
    uint32_t reject_fit_flags,
    void *workspace,
    uint32_t workspace_bytes,
    IAFocusSummary *summary);

#ifdef __cplusplus
}
#endif

#endif
