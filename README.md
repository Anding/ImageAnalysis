# ImageAnalysis
Astronomical image analysis in Forth

`compute-imageStats` builds the full 16-bit histogram and publishes integer
`MEAN`, `MEDIAN`, `MEDIANAD`, and `SATPIX` values into the frame's ordered
FITS metadata map. `SATPIX` is the exact number of pixels in histogram bin
65535; a saturation fraction remains derivable from the image dimensions.

## Native star analysis

`native\image_analysis.c` builds as the 32-bit `ImageAnalysisNative.dll` used
by VFX Forth. It operates directly on an in-memory unsigned 16-bit frame; file
loading, acquisition policy, and hardware control remain in Forth.

The native API:

- estimates a global median/MAD background;
- detects and measures stars;
- returns fixed-point per-star and frame statistics;
- summarizes any Forth-selected subset using a caller-supplied flag mask;
- runs stateful bitmap-driven autofocus with compact star selection,
  adaptive exposure, and intersecting-line fits;
- remeasures locked stars without changing their identities;
- robustly combines independently fitted stellar focus positions;
- returns an optional aggregate HFD diagnostic at the accepted focus.

The completed autofocus result reuses the exact exposure from the successful
discovery or sweep frame nearest the fitted focus, preferring the latest on
equal distances. Frames with rejected, saturated, or low-SNR locked-star
measurements are excluded; sweep exposure adaptation is unchanged.

All arrays and workspaces are caller-owned. Metric values use milli-units so
the ABI contains only 32-bit fields and is directly accessible from 32-bit
VFX Forth. HFD is twice HFR. `FWHM` is the Gaussian-equivalent width derived
from second moments and is identified as such because defocused stars need
not have Gaussian profiles.

Signal-to-noise uses signed background-subtracted aperture flux and all
aperture pixels. HFD still uses positive flux as defined above. Do not reuse
that rectified flux for SNR: noise-only apertures then appear significant and
can masquerade as surviving locked stars. The synthetic locked-noise
regression protects this distinction.

Autofocus updates a locked search center only from an accepted measurement.
Rejected/noise centroids must not walk a pending retake away from the star's
last valid position. A transient-loss/return regression protects recovery.

Overlapping locked searches can land on the same peak. A shared centroid
retains the first identity and flags subsequent identities as blended using
the supplied minimum separation, preventing duplicate stellar fit votes.

Build `ImageAnalysis_project\ImageAnalysis_project.vcxproj` for
`Release|Win32`. `native\test_image_analysis.c` contains dependency-free
synthetic tests for detection, measurement, exposure calculation, and ROI
calculation. Native autofocus receives in-memory bitmaps only; it contains no
image-file loader.

Native diagnostics are independent of Forth logging and are off by default.
`IA_SetDebug(1)` (or `ia.debug-on`) appends immediately flushed calculation
checkpoints to `E:\images\logs\image-analysis-native.log`;
`IA_SetDebug(0)` (or `ia.debug-off`) disables them.
Locked-star result checkpoints include SNR, local background/noise and
positive-aperture flux/area, so unexpectedly large faint-star HFD can be
investigated without changing acquisition or final-diagnostic exposure policy.
