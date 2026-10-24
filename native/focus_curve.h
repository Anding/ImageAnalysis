/* One active computation; star records are the only image-analysis boundary. */
#ifndef FOCUS_CURVE_H
#define FOCUS_CURVE_H

#include "image_analysis.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Reset the previous run. Follow up to max_stars (1..256); return 0 or -1. */
EXPORT int CALL start_computation(uint32_t max_stars);

/*
 * Copy matched HFDs from a zero-terminated StarRecord list.
 * focus_position: signed integer focuser steps. exposure: unsigned microseconds.
 * 0: accepted; 1: overexposed; 2: underexposed; -1: other failure.
 * 1/2 rewrite exposure and do not add a sample or update tracking: retry at
 * the same focus. Other returns leave exposure unchanged.
 * Input records can be reused/freed immediately. Keep ROI and binning fixed.
 */
EXPORT int CALL add_focus_point(const StarRecord *stars, int32_t focus_position,
    uint32_t *exposure);

/*
 * End the run and release its samples, on success or failure.
 * 0: write best_focus and error in integer focuser steps.
 * -1: insufficient/invalid fit or another failure; outputs remain unchanged.
 * error is a fit-quality diagnostic, not a calibrated confidence interval.
 */
EXPORT int CALL finish_computation(int32_t *best_focus, uint32_t *error);

#ifdef __cplusplus
}
#endif
#endif
