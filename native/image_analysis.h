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

#ifdef __cplusplus
}
#endif

#endif
