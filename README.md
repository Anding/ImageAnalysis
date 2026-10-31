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

`AnalyzeStars(pixels, width, height, median, mad, stars)` receives a
contiguous unsigned 16-bit bitmap, its logical dimensions, the already
calculated image-wide median and raw MAD in ADU, and the output address.
There is no detector limit argument or `starAnalysis.star-limit`: every call
returns up to 256 measured candidates. Negative results are -1 for invalid input
and -2 for allocation failure.
There is no native background calculation, calibrated SNR, frame summary, equipment
context or run state.

The Forth binding reads both values directly from `FRAME_STATISTICS`, without
reading FITS metadata or recalculating a histogram. C calculates
`max(MIN_SIGNAL_INCREMENT_ADU, MAD_MULTIPLIER * mad)` above the median.
This is a coordinated ABI change: six arguments and 104-byte records.
Update the DLL, ImageAnalysis bindings and ForthAstroFormats frame layout together,
and restart Forth after deployment. Old bindings or record allocations are unsafe.

`native\analysis_parameters.h` contains the compile-time detection, measurement,
tracking, curve-quality and exposure-guidance settings, with units and comments.
Rebuild the DLL after changing them; Forth and the C implementation do not duplicate
the tuning values. Mathematical constants and ABI capacities are not tuning
parameters. The following description uses the header's default values.

Detection uses eight-connected pixels **greater than median + max(16, 6 * MAD) ADU**, at
least three pixels per component. Each component supplies one candidate;
there is no deblending. Discovery scans the **whole bitmap**, accumulating
background-subtracted flux over each component's above-threshold pixels.
After measurement rejects invalid apertures and HFD outside 1..16 pixels,
overflow retains the brightest 256 candidates by this **discovery flux**,
breaking equal-flux ties by increasing X, then increasing Y. Candidates whose
discovery-centred aperture crosses the edge are excluded before selection.

Measurement precedes selection so invalid bright objects cannot consume slots.
A positive-flux centroid
refinement sets the circular 16-pixel measurement aperture. Candidates whose
aperture crosses an image edge are omitted. Nearby candidates can contaminate
each other's apertures; these initial rules do not classify blends.

Return at most **256** records, ordered by **measured aperture flux
descending, X ascending, Y ascending**, then a zero record. Saturation is not a
detection cutoff. Discovery flux is an overflow proxy: this does
not guarantee the globally brightest stars by aperture flux. Sub-threshold
wings can change the measured order.

The caller always
provides **257 records**, including the terminator; there is no buffer header
or separate capacity argument. Scratch storage is private to each call and freed before
return. On error the list starts with a zero record.

Brightness selection helps retain useful stars across focus frames, but list
rank is not identity: the independent focus curve matches stars by position.

The 104-byte, eight-byte-aligned `StarRecord` layout is declared in
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
| 72 | `roundness` | `double`; minor/major moment-axis width, 0..1 |
| 80 | `over_exposed_risk` | `double`; soft raw-peak headroom risk, 0..1 |
| 88 | `under_exposed_risk` | `double`; soft peak-significance risk, 0..1 |
| 96 | `exposure_scale` | `double`; calculated exposure multiplier |

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
int select_focus_stars(const StarRecord *stars, uint32_t *exposure);
int add_focus_point(const StarRecord *stars, int32_t focus_position,
                    uint32_t *exposure);
int finish_focus_point(void);
int finish_computation(int32_t *best_focus, uint32_t *error);
```

There is one active computation, not thread-safe. `start_computation` discards
the previous run and accepts the desired final survivor limit, 1..256.
`select_focus_stars` uses the entry-focus frame only for selection, not curve
samples. It ranks all candidates by equal default weights of roundness and
`min(1, flux / (2 * median_candidate_flux))`; coincident centroids are suppressed.
It freezes up to 256 identities and recommends the run's baseline exposure.
Later frames never populate new slots. Mutual nearest matching within eight
pixels skips ambiguous or competing detections. Keep ROI/binning fixed.

`add_focus_point` copies the best trial's complete matched HFD set at the current
position. Trial quality uses exposure risks in the top 32 selected identities,
including a penalty for missing detections. No clipping/roundness gate discards
curve samples. `finish_focus_point` commits the best available trial after
bounded retries, even when exposure targets could not be met. The input list
may immediately be reused or freed. Exposure is unsigned microseconds.
`add_focus_point` return values are:

| Return | Meaning |
| --- | --- |
| 0 | Desired low-risk coverage met, or no useful exposure change |
| 1 | Calculated exposure written; retake locally at the same focus |
| -1 | No measurement evidence or invalid input/allocation failure |

Exposure risks are heuristics, not calibrated probabilities: overexposure uses
raw ADC headroom; underexposure uses peak signal divided by `max(1, 1.4826*MAD)`.
Positive-only aperture flux is not treated as a calibrated SNR. The target peak
significance is 20 and peak headroom begins at 75% of 65535. Recommendations
use approximately linear raw-level scaling and squared significance scaling
for sky-dominated noise. Fixed camera bias, read noise and clipping weaken these
approximations: clipping gives only a shorter trial, not recovered true flux.
Exposure is bounded to 0.1..60 seconds. Forth accepts the selection recommendation
as the baseline and resets to it at every bracket position; local trials do not
change that baseline. These defaults require validation under the actual sky.

`finish_computation` commits pending samples, ends and frees the run. Each star
needs five distinct positions and whole-fit R² >= 0.8. The direct-HFD model is
`sqrt(h² + [s*(position-f)]²)`, with positive h/s and f inside the sampled bracket.
Repeated positions are averaged before fitting. No exposure-risk, saturation,
roundness or empirical contrast filter enters this stage. Choose up to the
requested number of valid curves in frozen selection rank, then combine their
focuses by the median.

E is the ceiling of the largest of one step, median individual focus standard
errors, and `1.4826 * MAD` of the stellar focus estimates. Individual errors use
the nonlinear-fit Jacobian covariance and residual variance with three fitted
parameters. E is a local diagnostic, not a calibrated confidence interval or
guarantee against systematic/model bias. A mathematically unavailable estimate
returns -1, leaves focus unchanged and reports E=`UINT32_MAX`; Forth retains
entry focus and continues. Only Forth gates an estimate by E (default <=50 steps).

Numerical background: [NIST least squares](https://www.itl.nist.gov/div898/handbook/pmd/section4/pmd431.htm)
and [Hocus Focus hyperbolic fitting](https://ghilios.github.io/hocus-focus/overview/hyperbola-fitting/).
The agreed three-parameter symmetric model is simpler than Hocus Focus's
four/five-parameter families; asymmetry can bias it. We do not adopt their
extra robust weighting/outlier gates, preserving the sample-count/R² boundary.

The curve module consumes only star records: it never invokes bitmap analysis,
controls equipment, writes FITS metadata or decides whether to apply focus.
Compile `native\focus_curve.c` with `native\test_focus_curve.c` alone for its
synthetic regression. `AstroImagingInForth\capabilities\Autofocus.f` is the
Forth caller: acquire through the existing pipeline, submit `FRAME_STARS`,
select once, retain best bounded trials at each position, and finish/apply or
retain entry focus according to E.

This is a replacement API, not an adapter. The previous `IA_*` analysis and
autofocus exports are removed. The legacy application autofocus script is no
longer automatically loaded.
