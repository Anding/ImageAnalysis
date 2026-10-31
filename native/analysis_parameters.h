/* Compile-time tuning; rebuild ImageAnalysisNative.dll after changing values.
 * Distances are in bitmap pixels, not arcseconds; keep ROI/binning fixed.
 * Mathematical constants and the StarRecord/buffer layout are not settings.
 */
#ifndef ANALYSIS_PARAMETERS_H
#define ANALYSIS_PARAMETERS_H

/* Discovery: signal must exceed max(this floor, MAD_MULTIPLIER * image MAD). */
#define MIN_SIGNAL_INCREMENT_ADU 16
#define MAD_MULTIPLIER 6
#define MIN_COMPONENT_PIXELS 3          /* Eight-connected pixels above threshold. */

/* Measurement: aperture radius and accepted HFD are separate choices. */
#define APERTURE_RADIUS_PIXELS 16       /* Positive integer; larger admits more neighbours. */
#define RADIAL_BINS_PER_PIXEL 16        /* Positive integer; HFD radial resolution. */
#define MIN_HFD_PIXELS 1.0              /* Reject single-pixel/hot-pixel profiles. */
#define MAX_HFD_PIXELS 16.0             /* Retain defocused stars, reject broad sources. */
#define MAX_SATURATED_PIXELS 3          /* Raw pixels equal to 65535 per aperture. */

/* Autofocus: centroid tracking and minimum evidence for a bracketed V curve. */
#define MATCH_RADIUS_PIXELS 8.0
#define MIN_FOCUS_POSITIONS 5           /* At least five, for positive fit degrees of freedom. */
#define MIN_BRANCH_POINTS 2            /* At least two observations per fitted line. */
#define MIN_BRANCH_RISE_PIXELS 0.5
#define MIN_BRANCH_RISE_FRACTION 0.1    /* Fraction of fitted minimum HFD; 0.1 is 10%. */

/* Exposure guidance applies only when no usable match survives clipping. */
#define EXPOSURE_REDUCTION_DIVISOR 2    /* Integer greater than one. */
#define MIN_EXPOSURE_US 100000          /* Recommended exposure floor, in microseconds. */

#if MIN_FOCUS_POSITIONS < 5 || MIN_BRANCH_POINTS < 2
#error Focus fitting requires at least five positions and two points per branch.
#endif

#endif
