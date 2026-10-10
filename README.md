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
- calculates a revised exposure and an ideal multi-star subframe.
- remeasures locked stars without changing their identities;
- fits saved piecewise-V or hyperbolic focus models;
- robustly combines independently fitted stellar focus positions.

All arrays and workspaces are caller-owned. Metric values use milli-units so
the ABI contains only 32-bit fields and is directly accessible from 32-bit
VFX Forth. HFD is twice HFR. `FWHM` is the Gaussian-equivalent width derived
from second moments and is identified as such because defocused stars need
not have Gaussian profiles.

Build `ImageAnalysis_project\ImageAnalysis_project.vcxproj` for
`Release|Win32`. `native\test_image_analysis.c` contains dependency-free
synthetic tests for detection, measurement, exposure calculation, and ROI
calculation. `native\focus_corpus.c` is an offline XISF corpus tool; file
loading exists only in that test utility, not in the DLL.

Initial measurements of the 2026-10-10 focus corpus remeasured seven locked
stars in roughly 0-1 ms per 2504x2504 frame. Four clear brackets produced
median per-star model residuals of 0.36-0.58 pixels and focus estimates from
5234 to 5253. The cloudy bracket produced residuals above 2.3 pixels, giving
Forth policy a strong numerical rejection signal without embedding that
decision in C.
