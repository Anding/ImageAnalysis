/* One active computation; star records are the only image-analysis boundary. */
#ifndef FOCUS_CURVE_H
#define FOCUS_CURVE_H

#include "image_analysis.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Reset the previous run. max_stars is the requested final fit count (1..256). */
EXPORT int CALL start_computation(uint32_t max_stars);

/*
 * Freeze the ranked identity pool from the first frame and set its exposure
 * baseline. Saturation is not an admission gate. Returns 0 on selection and
 * -1 for invalid inputs or an empty selection. Called once per run.
 */
EXPORT int CALL select_focus_stars(const StarRecord *stars, uint32_t *exposure);

/*
 * Stage this exposure trial for a focus position. High-risk positive finite
 * HFDs remain eligible samples; risks only guide exposure and trial choice.
 * 0: no useful exposure change; 1: wrote a recommended exposure; -1: invalid
 * input or no matched sample. Better trials replace worse trials until the
 * point is finished. Exposure is bounded by the configured minimum/maximum.
 * Inputs may be reused immediately. Keep ROI and binning fixed.
 */
EXPORT int CALL add_focus_point(const StarRecord *stars, int32_t focus_position,
    uint32_t *exposure);

/* Commit the best staged trial and clear its state; -1 if no trial has evidence. */
EXPORT int CALL finish_focus_point(void);

/*
 * End the run and release its samples, on success or failure.
 * Commits a pending trial. 0 writes the best_focus and fit uncertainty in
 * integer focuser steps. -1 leaves best_focus unchanged and writes UINT32_MAX
 * to error when available. Error is not a calibrated real-world guarantee.
 */
EXPORT int CALL finish_computation(int32_t *best_focus, uint32_t *error);

#ifdef __cplusplus
}
#endif
#endif
