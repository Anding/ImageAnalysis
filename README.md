# ImageAnalysis

`NEED ImageAnalysis` loads the complete repository through `Library.f`:
pixel/statistics, display functions, stateless stars and the stateful focus-curve
ABI. It depends on Formats, not on the application. These remain separate source
files, not independently loadable packages.

`ImageAnalysis.f` contains the Forth histogram, pixel-statistics and image
combination code. `compute-imageStats` writes only the existing pixel FITS
keywords `MEAN`, `MEDIAN`, `MEDIANAD` and `SATPIX` (plus their section marker).
It does not invoke star analysis.

`StarAnalysis.f` is the separate native star-analysis binding and caller.
Pass a caller-owned `frame` explicitly; no global `image` is defined:

```forth
need ImageAnalysis
frame compute-imageStats
frame analyze-stars
frame FRAME_STARS @  \ Address of the headerless, zero-terminated StarRecord list.
```

`analyze-stars ( frame -- )` reads `MEDIAN` from the frame metadata,
converts it to an unsigned 16-bit integer, and supplies that background to C.
Missing, malformed or out-of-range metadata supplies zero. It neither passes
nor scans the histogram. The caller owns consistent stage selection.
It allocates no result/statistics storage and writes **no FITS
metadata**. Native failures leave an empty list and are logged without aborting.
See `StarAnalysis_test1.f` for an executable synthetic-frame tutorial.
`DisplayFunction.f` remains the histogram-based display stretch.

## Native contract

`AnalyzeStars(pixels, width, height, median, mad, max_stars, stars)` receives a
contiguous unsigned 16-bit bitmap, its logical dimensions, the already
calculated image-wide median and raw MAD in ADU, a per-call measurement limit, and the output
address. `StarAnalysis.f` defines `starAnalysis.star-limit` with default 10.
Callers can change it directly; `analyze-stars` reads its current value and
passes it to C. Temporary overrides are saved/restored by their callers.
Neither `<expose>` nor `pipeline-process` takes an analysis-limit argument.
C stores no configured limit.
`max_stars` accepts 0..256; zero returns an empty list. It returns the number
of measured stars, or a negative error (`-1` invalid input, including a limit
above 256; `-2` allocation failure).
There is no native background calculation, SNR, frame summary, equipment
context or run state.

The Forth binding reads both values directly from `FRAME_STATISTICS`, without
reading FITS metadata or recalculating a histogram. C calculates
`max(MIN_SIGNAL_INCREMENT_ADU, MAD_MULTIPLIER * mad)` above the median.
The seven-argument layout and `StarRecord` are unchanged, but the fifth argument
now means **raw MAD**, not a precomputed threshold. Update the DLL and Forth
binding together, and restart Forth after deployment.

`native\analysis_parameters.h` contains the compile-time detection, measurement,
tracking, curve-quality and exposure-guidance settings, with units and comments.
Rebuild the DLL after changing them; Forth and the C implementation do not duplicate
the tuning values. Mathematical constants and ABI capacities are not tuning
parameters. The following description uses the unchanged default values.

Detection uses eight-connected pixels **greater than median + max(16, 6 * MAD) ADU**, at
least three pixels per component. Each component supplies one candidate;
there is no deblending. Discovery scans the **whole bitmap**, accumulating
background-subtracted flux over each component's above-threshold pixels.
After measurement rejects invalid apertures and HFD outside 1..16 pixels,
selection prefers candidates with at most three saturated pixels, then
retains the brightest `max_stars` candidates by this **discovery flux**,
breaking equal-flux ties by increasing X, then increasing Y. Candidates whose
discovery-centred aperture crosses the edge are excluded before selection.

Measurement precedes selection so invalid bright objects cannot consume slots.
A positive-flux centroid
refinement sets the circular 16-pixel measurement aperture. Candidates whose
aperture crosses an image edge are omitted. Nearby candidates can contaminate
each other's apertures; these initial rules do not classify blends.

Return at most **max_stars** records, ordered by **measured aperture flux
descending, X ascending, Y ascending**, then a zero record. Clipped candidates
remain available as fallbacks for exposure guidance. Discovery flux is a selection proxy: this does
not guarantee the globally brightest stars by aperture flux. Sub-threshold
wings can change the measured order.

The caller always
provides **257 records**, including the terminator; there is no buffer header
or separate capacity argument. Scratch storage is private to each call and freed before
return. On error the list starts with a zero record.

Brightness selection helps retain useful stars across focus frames, but list
rank is not identity: the independent focus curve matches stars by position.

The 72-byte, eight-byte-aligned `StarRecord` layout is declared in
`native\image_analysis.h` and in `ForthAstroFormats\Frame.f`:

| Offset | Field | Representation |
| --- | --- | --- |
| 0 | `positive_pixel_count` | `uint32_t`; zero terminates the list |
| 4 | `saturated_pixel_count` | `uint32_t`; raw samples equal to 65535 |
| 8, 16 | `x`, `y` | `double`; positive-flux centroid, zero-based bitmap pixel centres |
| 24 | `flux` | `double`; sum of positive background-subtracted ADU |
| 32 | `hfd` | `double`; half-flux diameter in bitmap pixels |
| 40, 48, 56 | `qxx`, `qyy`, `qxy` | `double`; positive-flux central second moments in pixels squared |
| 64 | `max_pixel` | `uint16_t`; maximum raw aperture sample; six tail-padding bytes |

Use `DF@` for C doubles in VFX Forth, **not `F@`** (80-bit extended).
Qxx=Qyy alone does not imply roundness: Qxy must also vanish. HFD comes from
positive pixel flux accumulated in 1/16-pixel radial bins about the returned
centroid, with linear interpolation in the half-flux bin, not from Q.
The flux centroid and central moments describe the pixels selected by the
refined aperture; the centroid can move slightly as its boundary changes.

Every frame always allocates its own zeroed statistics block and initially
empty star list, and frees both with `free-frame`. `compute-imageStats` reuses
that statistics block rather than reallocating it.
`FRAME_STARS` never refers to shared native result memory;
a later analysis of another frame does not overwrite it. ROI dimensions do
not change the frame's bitmap allocation capacity or its star-list allocation.

Build `ImageAnalysis_project\ImageAnalysis_project.vcxproj` as `Release|Win32`
and deploy `ImageAnalysisNative.dll` to MPE's Bin directory. Compile
`native\image_analysis.c` with `native\test_image_analysis.c` to run the
dependency-free numerical regression.

## Independent native focus curve

The same DLL exports the API in `native\focus_curve.h`, bound directly by
`StarAnalysis.f`:

```c
int start_computation(uint32_t max_stars);
int add_focus_point(const StarRecord *stars, int32_t focus_position,
                    uint32_t *exposure);
int finish_computation(int32_t *best_focus, uint32_t *error);
```

There is one active computation, not thread-safe. `start_computation` discards
any previous run and accepts 1..256 identities to follow; it returns 0 or -1.
The first accepted list establishes up to that many usable identities, in list
order. Later lists match mutual nearest centroids within eight bitmap pixels;
ties and competing matches are omitted. Keep ROI/binning fixed. List rank can
change, and missing identities are never replaced by newly discovered stars.
Unused slots may acquire later detections, with earlier samples left missing.
Clipped detections do not displace nearby usable matches.

`add_focus_point` copies matched HFDs and the actual signed focuser-step
position; the caller may immediately reuse or free the input list. Partial
matches are accepted; no matches returns -1, not an underexposure diagnosis.
Exposure is unsigned integer microseconds. Return values are:

| Return | Meaning |
| --- | --- |
| 0 | Sample accepted; exposure unchanged |
| 1 | Only clipped matches available; exposure halved; retry at the same focus |
| -1 | Other failure; exposure unchanged |

Each usable identity needs HFD in 1..16 pixels and at most three saturated
pixels. Dim stars and isolated saturated pixels do not trigger exposure changes.
Other usable stars allow a sample to proceed despite clipped detections.
Only when no usable match remains does clipping recommend half the exposure,
with a 100 ms floor; unchanged or increasing recommendations fail. Rejected samples
neither enter the fit nor establish/update tracking.

`finish_computation` ends and frees the run, returning 0 or -1. Each usable
star needs at least five distinct focus positions. Fit two least-squares
linear branches, with at least two observations per branch, negative/positive
slopes and their intersection bracketed between the inner branch samples.
Each branch must rise by at least 0.5 pixel or 10% of the fitted minimum,
whichever is larger, to reject nearly focus-independent sources.
Choose the valid split with lowest residual sum of squares; do not extrapolate.
Combine usable stellar focus estimates by their median. Return the rounded
signed integer best focus and an unsigned integer uncertainty in focuser steps:
the ceiling of the largest of one step, median propagated fit error, and
1.4826 times the median absolute deviation between stellar focus estimates.
This is a diagnostic, not a calibrated confidence interval. No usable fit
returns -1 and leaves both outputs unchanged.
These are aperture and curve-quality filters, not a galaxy classifier.
Compact nuclei can resemble stellar profiles; offline checks do not establish
reliable autofocus under every sky condition.

The curve module consumes only star records: it never invokes bitmap analysis,
controls equipment, writes FITS metadata or decides whether to apply focus.
Compile `native\focus_curve.c` with `native\test_focus_curve.c` alone for its
synthetic regression. `AstroImagingInForth\capabilities\Autofocus.f` is the
Forth caller: acquire through the existing pipeline, submit `FRAME_STARS`,
retry exposure at the same position, and finish/apply or restore prior focus.

This is a replacement API, not an adapter. The previous `IA_*` analysis and
autofocus exports are removed. The legacy application autofocus script is no
longer automatically loaded.
