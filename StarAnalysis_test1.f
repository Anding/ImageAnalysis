\ Executable tutorial: frame-owned native records, independent of FITS output.
need ImageAnalysis
need simple-tester

64 64 1 allocate-frame constant stars.test.frame
64 64 1 allocate-frame constant stars.test.other
stars.test.frame FRAME_STARS @ constant stars.test.list
stars.test.frame FRAME_STATISTICS @ constant stars.test.stats

: stars.test.background { frame -- }
\ Fill a frame with a constant 1000-ADU sky.
    4096 0 do 1000 frame FRAME_BITMAP i 2* + w! loop
;

: stars.test.pixel { value x y -- }
\ Write one synthetic pixel into the first frame.
    value stars.test.frame FRAME_BITMAP y 64 * x + 2* + w!
;

: stars.test.triangle ( -- )
\ Three connected pixels with an exactly known fractional centroid and moments.
    1100 32 32 stars.test.pixel
    1200 33 32 stars.test.pixel
    1100 32 33 stars.test.pixel
;

Tstart
T{ StarRecord star-buffer-bytes }T 104 26728 ==
\ Every frame has its own empty list, even if star analysis is never requested.
T{ stars.test.list 0<> }T -1 ==
T{ stars.test.other FRAME_STARS @ stars.test.list <> }T -1 ==
T{ stars.test.list star.positive-pixel-count @ }T 0 ==
T{ stars.test.other FRAME_STARS @ star.positive-pixel-count @ }T 0 ==

stars.test.frame stars.test.background
stars.test.triangle
\ Pixel statistics alone must not invoke star analysis.
stars.test.frame compute-imageStats
T{ stars.test.frame FRAME_STATISTICS @ }T stars.test.stats ==
T{ stars.test.list star.positive-pixel-count @ }T 0 ==
T{ stars.test.frame analyze-stars }T ==
T{ stars.test.frame FRAME_STARS @ }T stars.test.list ==
T{ stars.test.list star.positive-pixel-count @ }T 3 ==
T{ stars.test.list StarRecord + star.positive-pixel-count @ }T 0 ==
T{ stars.test.list star.max-pixel w@ }T 1200 ==
T{ stars.test.list star.x df@ f>s }T 32 ==
T{ stars.test.list star.y df@ f>s }T 32 ==
T{ stars.test.list star.flux df@ f>s }T 400 ==
T{ stars.test.list star.qxx df@ 1000e f* f>s }T 250 ==
T{ stars.test.list star.qyy df@ 10000e f* f>s }T 1875 ==
T{ stars.test.list star.qxy df@ 1000e f* f>s }T -125 ==
T{ stars.test.list star.hfd df@ 1e f> }T -1 ==
T{ stars.test.list star.saturated-pixel-count @ }T 0 ==
T{ stars.test.list star.roundness df@ 0e f> }T -1 ==
T{ stars.test.list star.over-exposed-risk df@ f>s }T 0 ==
T{ stars.test.list star.under-exposed-risk df@ f>s }T 0 ==
T{ stars.test.list star.exposure-scale df@ 0e f> }T -1 ==
\ Existing pixel FITS keys are unchanged; no star keys are added.
T{ stars.test.frame FRAME_METADATA @ map-count }T 5 ==
T{ s" NSTARS" stars.test.frame FRAME_METADATA @ map@ }T 0 0 ==
T{ s" HFD" stars.test.frame FRAME_METADATA @ map@ }T 0 0 ==

\ Existing statistics suffice after erasing the FITS map; no metadata is required.
stars.test.other stars.test.background
stars.test.other compute-imageStats
stars.test.other FRAME_METADATA @ erase-map
T{ stars.test.other analyze-stars }T ==
T{ stars.test.other FRAME_METADATA @ map-count }T 0 ==
T{ stars.test.other FRAME_STARS @ star.positive-pixel-count @ }T 0 ==
T{ stars.test.list star.positive-pixel-count @ }T 3 ==

\ Even misleading metadata cannot replace the frame's statistical median.
s" invalid" stars.test.frame FRAME_METADATA @ =>" MEDIAN"
T{ stars.test.frame analyze-stars }T ==
T{ stars.test.list star.positive-pixel-count @ }T 3 ==
\ Pass raw MAD: 10 admits this triangle; 20 leaves only one bright pixel.
10 stars.test.stats MEDIAN_ABSOLUTE_DEVIATION !
T{ stars.test.frame analyze-stars }T ==
T{ stars.test.list star.positive-pixel-count @ }T 3 ==
20 stars.test.stats MEDIAN_ABSOLUTE_DEVIATION !
T{ stars.test.frame analyze-stars }T ==
T{ stars.test.list star.positive-pixel-count @ }T 0 ==
0 stars.test.stats MEDIAN_ABSOLUTE_DEVIATION !

\ The existing hot-pixel rejection also excludes a peak with two weak neighbours.
65535 33 32 stars.test.pixel
stars.test.frame analyze-stars
T{ stars.test.list star.positive-pixel-count @ }T 0 ==
\ A native failure empties the list and returns without aborting the run.
0 stars.test.frame FRAME_WIDTH !
T{ stars.test.frame analyze-stars }T ==
T{ stars.test.list star.positive-pixel-count @ }T 0 ==
64 stars.test.frame FRAME_WIDTH !
stars.test.frame stars.test.background
stars.test.frame compute-imageStats
T{ stars.test.frame FRAME_STATISTICS @ }T stars.test.stats ==
stars.test.frame analyze-stars
T{ stars.test.list star.positive-pixel-count @ }T 0 ==
\ Detection keeps both candidates; autofocus owns its independent survivor count.
1100 18 18 stars.test.pixel
1200 19 18 stars.test.pixel
1100 18 19 stars.test.pixel
2100 45 45 stars.test.pixel
3200 46 45 stars.test.pixel
2100 45 46 stars.test.pixel
T{ stars.test.frame analyze-stars }T ==
T{ stars.test.list star.x df@ f>s }T 45 ==
T{ stars.test.list StarRecord + star.x df@ f>s }T 18 ==
T{ stars.test.list StarRecord 2* + star.positive-pixel-count @ }T 0 ==
T{ stars.test.frame free-frame stars.test.other free-frame }T ==
Tend
