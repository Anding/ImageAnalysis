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
/* Selection weights are preferences, not probabilities or shape rejection gates. */
#define ROUNDNESS_WEIGHT 0.5
#define BRIGHTNESS_WEIGHT 0.5
#define DUPLICATE_RADIUS_PIXELS 1.0     /* Suppress coincident measured identities. */
#define MATCH_RADIUS_PIXELS 8.0

/* Exposure heuristics: peak significance is not a calibrated stellar SNR.
 * MAD -> Gaussian sigma; with sky-dominated noise significance grows as sqrt(t).
 * Linear raw-level scaling assumes negligible fixed bias. Clipping requires trials.
 */
#define TARGET_PEAK_SIGMA 20.0
#define MIN_NOISE_ADU 1.0
#define SAFE_PEAK_FRACTION 0.75         /* Raw full-scale headroom starts the soft risk. */
#define TARGET_PEAK_FRACTION 0.5        /* Selection-frame baseline target. */
#define AF_EXPOSURE_TOP_STARS 32
#define MAX_EXPOSURE_RISK 0.5
#define MIN_EXPOSURE_CHANGE_FRACTION 0.05
#define MIN_EXPOSURE_US 100000
#define MAX_EXPOSURE_US 60000000

/* Direct three-parameter hyperbola fit; only distinct count and R² filter quality. */
#define MIN_FOCUS_POSITIONS 5
#define MIN_FIT_R_SQUARED 0.8
#define FIT_MAX_ITERATIONS 100
#define FIT_RELATIVE_TOLERANCE 1e-9
#define FIT_INITIAL_DAMPING 1e-3

#if MIN_FOCUS_POSITIONS < 4
#error Three-parameter fit needs positive residual degrees of freedom.
#endif

#endif
