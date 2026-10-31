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
} StarRecord;
#pragma pack(pop)

/*
 * pixels: contiguous width*height unsigned 16-bit samples.
 * median, mad: previously calculated image-wide statistics, in ADU.
 * C calculates the signal threshold from analysis_parameters.h:
 * max(MIN_SIGNAL_INCREMENT_ADU, MAD_MULTIPLIER * mad), above median.
 * The fifth argument is raw MAD, not a precomputed threshold.
 * max_stars: 0..MAX_STARS; with default parameters, prefer candidates with at most
 * three saturated pixels, then select by discovery flux after rejecting invalid apertures and HFD
 * outside 1..16 pixels. No deblending or galaxy classification:
 * compact galaxy nuclei can resemble stars and require focus-curve validation.
 * stars: space for MAX_STARS records plus a zero terminator.
 * Records sort by measured flux descending, then x ascending, then y ascending.
 * Returns the count (0..max_stars), -1 for invalid inputs, -2 for allocation failure.
 * On error the output list is empty. Inputs are never modified.
 */
EXPORT int CALL AnalyzeStars(const uint16_t *pixels, uint32_t width,
    uint32_t height, uint32_t median, uint32_t mad, uint32_t max_stars,
    StarRecord *stars);

#ifdef __cplusplus
}
#endif
#endif
