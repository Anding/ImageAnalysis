# ImageAnalysis

`ImageAnalysis.f` contains the Forth histogram, pixel-statistics and image
combination code. `compute-imageStats` writes only the existing pixel FITS
keywords `MEAN`, `MEDIAN`, `MEDIANAD` and `SATPIX` (plus their section marker).
It does not invoke star analysis.

`StarAnalysis.f` is the separate native star-analysis binding and caller:

```forth
need StarAnalysis
image compute-imageStats
image 10 analyze-stars
image FRAME_STARS @  \ Address of the headerless, zero-terminated StarRecord list.
```

`analyze-stars ( frame max_stars -- )` reads `MEDIAN` from the frame metadata,
converts it to an unsigned 16-bit integer, and supplies that background to C.
Missing, malformed or out-of-range metadata supplies zero. It neither passes
nor scans the histogram. The caller owns consistent stage selection.
It allocates no result/statistics storage and writes **no FITS
metadata**. Native failures leave an empty list and are logged without aborting.
See `StarAnalysis_test1.f` for an executable synthetic-frame tutorial.
`DisplayFunction.f` remains the histogram-based display stretch.

## Native contract

`AnalyzeStars(pixels, width, height, background, max_stars, stars)` receives a
contiguous unsigned 16-bit bitmap, its logical dimensions, the already
calculated background in ADU, a per-call measurement limit, and the output
address. The application owns `starAnalysis.default.maxStars` (environment
value 10). `expose` reads it and passes it through `<expose>` and
`pipeline-process` to `analyze-stars`; explicit callers supply their own limit.
C stores no configured limit.
`max_stars` accepts 0..256; zero returns an empty list. It returns the number
of measured stars, or a negative error (`-1` invalid input, including a limit
above 256; `-2` allocation failure).
There is no native background calculation, SNR, frame summary, equipment
context or run state.

Detection uses eight-connected pixels **greater than median + 16 ADU**, at
least three pixels per component. Each component supplies one candidate;
there is no deblending. Discovery scans the **whole bitmap**, accumulating
background-subtracted flux over each component's above-threshold pixels.
It retains the brightest `max_stars` candidates by this **discovery flux**,
breaking equal-flux ties by increasing X, then increasing Y. Candidates whose
discovery-centred aperture crosses the edge are excluded before selection.

Measurement is a separate phase, applied only to the selected candidates.
A positive-flux centroid
refinement sets the circular 16-pixel measurement aperture. Candidates whose
aperture crosses an image edge are omitted. Nearby candidates can contaminate
each other's apertures; these initial rules do not classify blends.

Return at most **max_stars** records, ordered by **measured aperture flux
descending, X ascending, Y ascending**, then a zero record. A candidate rejected
during measurement is omitted without measuring a replacement, so fewer
records can be returned. Discovery flux is a fast selection proxy: this does
not guarantee the globally brightest stars by aperture flux. Sub-threshold
wings can change the measured order.

The caller always
provides **257 records**, including the terminator; there is no buffer header
or separate capacity argument. Scratch storage is private to each call and freed before
return. On error the list starts with a zero record.

Brightness selection helps retain useful stars across focus frames, but list
rank is not identity: a future autofocus curve still needs positional matching.

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

This is a replacement API, not an adapter. The previous `IA_*` analysis and
autofocus exports are removed. The legacy application autofocus script is no
longer automatically loaded. A separate stateful curve modeller and its
Forth autofocus orchestration are not implemented here.
