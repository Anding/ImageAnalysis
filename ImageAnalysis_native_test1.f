need ImageAnalysis
need simple-tester

64 64 1 allocate-frame constant native.test.frame

: native.test.background ( -- )
    64 64 * 0 do
        1000 native.test.frame FRAME_BITMAP i 2* + w!
    loop
;

: native.test.pixel { value x y -- }
    value native.test.frame FRAME_BITMAP y 64 * x + 2* + w!
;

: native.test.star { cx cy -- }
    12000 cx 2 - cy native.test.pixel
    12000 cx 2 + cy native.test.pixel
    12000 cx cy 2 - native.test.pixel
    12000 cx cy 2 + native.test.pixel
    24000 cx 1 - cy native.test.pixel
    24000 cx 1 + cy native.test.pixel
    24000 cx cy 1 - native.test.pixel
    24000 cx cy 1 + native.test.pixel
    50000 cx cy native.test.pixel
;

native.test.background
32 31 native.test.star
native.test.frame compute-imageStats

create native.test.model <IA_FOCUS_MODEL> allot
<IA_FOCUS_MODEL> native.test.model IA_FOCUS_MODEL.struct_size !
IA.FOCUS.PIECEWISE_V native.test.model IA_FOCUS_MODEL.kind !
38000 native.test.model IA_FOCUS_MODEL.left_slope_micro_per_step !
42000 native.test.model IA_FOCUS_MODEL.right_slope_micro_per_step !

create native.test.samples <IA_FOCUS_SAMPLE> 3 * allot

: native.test.sample { focus metric index -- }
    focus index <IA_FOCUS_SAMPLE> * native.test.samples +
        IA_FOCUS_SAMPLE.focus_position !
    metric index <IA_FOCUS_SAMPLE> * native.test.samples +
        IA_FOCUS_SAMPLE.metric_milli !
    1000 index <IA_FOCUS_SAMPLE> * native.test.samples +
        IA_FOCUS_SAMPLE.weight_milli !
;

4980 3520 0 native.test.sample
5020 2000 1 native.test.sample
5060 3680 2 native.test.sample
create native.test.fit <IA_FOCUS_FIT> allot

Tstart
T{ IA_Version }T 1 ==
T{ ia.summary IA_FRAME_SUMMARY.detected_stars @ 0> }T -1 ==
T{ ia.summary IA_FRAME_SUMMARY.included_stars @ }T 1 ==
T{ 0 ia.star IA_STAR.x_milli @ 31000 33000 within }T -1 ==
T{ 0 ia.star IA_STAR.y_milli @ 30000 32000 within }T -1 ==
T{ s" NSTARS" native.test.frame FRAME_METADATA @ >string hashS }T s" 1" hashS ==
T{ s" HFD" native.test.frame FRAME_METADATA @ >string nip 0> }T -1 ==
T{ 1000 21000 1000 41000 100 10000 4000 ia.recommend-exposure }T 2000 0 ==
T{ 0 ia.star 1 64 64 5 ia.compute-roi
    swap IA_RECT.width @ swap }T 11 0 ==
T{ native.test.model native.test.samples 3 0 250 native.test.fit
    ia.fit-focus }T 0 ==
T{ native.test.fit IA_FOCUS_FIT.focus_milli_steps @
    5019000 5021001 within }T -1 ==
Tend

native.test.frame free-frame
