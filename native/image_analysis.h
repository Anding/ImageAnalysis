/* Stateless star measurements. Keep this layout aligned with Frame.f. */
#ifndef IMAGE_ANALYSIS_H
#define IMAGE_ANALYSIS_H

#include <stdint.h>

#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#define CALL __cdecl
#else
#define EXPORT
#define CALL
#endif

#ifdef __cplusplus
extern "C" {
#endif

enum { MAX_STARS = 256 };

#pragma pack(push, 8)
typedef struct StarRecord {
    uint32_t positive_pixel_count; /* zero terminates the list */
    uint32_t saturated_pixel_count;
    double x, y;                  /* bitmap-local pixel centres */
    double flux;                  /* positive background-subtracted ADU */
    double hfd;                   /* pixels */
    double qxx, qyy, qxy;          /* central second moments, pixels squared */
    uint16_t max_pixel;            /* raw sample; six bytes of tail padding */
    double roundness;              /* minor/major moment-axis width, 0..1 */
    double over_exposed_risk;       /* soft peak-headroom risk, 0..1 */
    double under_exposed_risk;      /* soft peak-significance risk, 0..1 */
    double exposure_scale;          /* calculated multiplier; clipping only gives a trial */
} StarRecord;
#pragma pack(pop)

/*
 * pixels: contiguous width*height unsigned 16-bit samples.
 * median, mad: previously calculated image-wide statistics, in ADU.
 * C calculates the signal threshold from analysis_parameters.h:
 * max(MIN_SIGNAL_INCREMENT_ADU, MAD_MULTIPLIER * mad), above median.
 * The fifth argument is raw MAD, not a precomputed threshold.
 * Return up to MAX_STARS; retain by discovery flux on overflow, with no clipping
 * cutoff. Reject invalid apertures and HFD outside the configured range.
 * No deblending or galaxy classification:
 * compact galaxy nuclei can resemble stars and require focus-curve validation.
 * stars: space for MAX_STARS records plus a zero terminator.
 * Records sort by measured flux descending, then x ascending, then y ascending.
 * Returns the count (0..MAX_STARS), -1 for invalid inputs, -2 for allocation failure.
 * On error the output list is empty. Inputs are never modified.
 */
EXPORT int CALL AnalyzeStars(const uint16_t *pixels, uint32_t width,
    uint32_t height, uint32_t median, uint32_t mad,
    StarRecord *stars);

#ifdef __cplusplus
}
#endif
#endif
