\ Executable tutorial: frame-owned native records, independent of FITS output.
need StarAnalysis
need simple-tester

64 64 1 allocate-frame constant stars.test.frame
64 64 1 allocate-frame constant stars.test.other
stars.test.frame FRAME_STARS @ constant stars.test.list

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
\ Every frame has its own empty list, even if star analysis is never requested.
T{ stars.test.list 0<> }T -1 ==
T{ stars.test.other FRAME_STARS @ stars.test.list <> }T -1 ==
T{ stars.test.list star.positive-pixel-count @ }T 0 ==
T{ stars.test.other FRAME_STARS @ star.positive-pixel-count @ }T 0 ==

stars.test.frame stars.test.background
stars.test.triangle
\ Pixel statistics alone must not invoke star analysis.
stars.test.frame compute-imageStats
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
\ Existing pixel FITS keys are unchanged; no star keys are added.
T{ stars.test.frame FRAME_METADATA @ map-count }T 5 ==
T{ s" NSTARS" stars.test.frame FRAME_METADATA @ map@ }T 0 0 ==
T{ s" HFD" stars.test.frame FRAME_METADATA @ map@ }T 0 0 ==

\ The caller prepares a histogram even without compute-imageStats.
stars.test.other stars.test.background
T{ stars.test.other analyze-stars }T ==
T{ stars.test.other FRAME_METADATA @ map-count }T 0 ==
T{ stars.test.other FRAME_STARS @ star.positive-pixel-count @ }T 0 ==
T{ stars.test.list star.positive-pixel-count @ }T 3 ==

65535 33 32 stars.test.pixel
stars.test.frame analyze-stars
T{ stars.test.list star.max-pixel w@ }T 65535 ==
T{ stars.test.list star.saturated-pixel-count @ }T 1 ==
stars.test.frame stars.test.background
stars.test.frame analyze-stars
T{ stars.test.list star.positive-pixel-count @ }T 0 ==
T{ stars.test.frame free-frame stars.test.other free-frame }T ==
Tend
