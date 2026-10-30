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
 * 0: accepted; 1: all usable matches clipped; -1: no usable sample/other failure.
 * A star needs HFD in 1..16 pixels and at most three saturated pixels. Clipped
 * stars are omitted when other stars remain. No raw-peak brightness target.
 * Unused tracking slots may acquire later detections, without reassigning
 * missing identities. Clipped detections do not hide nearby usable matches.
 * 1 halves exposure (minimum 100 ms), without adding a sample or updating
 * tracking: retry at the same focus. Other returns leave exposure unchanged.
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
