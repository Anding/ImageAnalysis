\ Native-binding regression using a deterministic synthetic star field.
\ Maintenance: exercise the public Forth wrappers here; numerical edge cases
\ that do not cross the ABI belong in native\test_image_analysis.c.

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

create native.test.measurement <IA_STAR> allot
create native.test.focus-config <IA_FOCUS_CONFIG> allot
native.test.focus-config IA_DefaultFocusConfig
1 native.test.focus-config IA_FOCUS_CONFIG.star_limit !
create native.test.focus-request <IA_FOCUS_REQUEST> allot
IA_FocusRunBytes constant native.test.focus-run-bytes
native.test.focus-run-bytes allocate throw constant native.test.focus-run

Tstart
T{ IA_Version }T 2 ==
T{ ia.debug-on ia.debug-off IA_Version }T 2 ==
T{ ia.summary IA_FRAME_SUMMARY.detected_stars @ 0> }T -1 ==
T{ ia.summary IA_FRAME_SUMMARY.included_stars @ }T 1 ==
T{ 0 ia.star IA_STAR.x_milli @ 31000 33000 within }T -1 ==
T{ 0 ia.star IA_STAR.y_milli @ 30000 32000 within }T -1 ==
T{ s" NSTARS" native.test.frame FRAME_METADATA @ >string hashS }T s" 1" hashS ==
T{ s" HFD" native.test.frame FRAME_METADATA @ >string nip 0> }T -1 ==
T{ 1000 21000 1000 41000 100 10000 4000 ia.recommend-exposure }T 2000 0 ==
T{ 0 ia.star 1 64 64 5 ia.compute-roi
    swap IA_RECT.width @ swap }T 11 0 ==
T{ native.test.frame 0 ia.star 1 5 native.test.measurement
    ia.measure-locked }T 0 ==
T{ native.test.measurement IA_STAR.hfd_milli_pixels @ 0> }T -1 ==
T{ native.test.focus-run native.test.focus-run-bytes
    native.test.focus-config native.test.frame 5000 1000
    native.test.focus-request ia.begin-focus }T 0 ==
T{ native.test.focus-request IA_FOCUS_REQUEST.action @
    }T IA.FOCUS.CAPTURE ==
T{ native.test.focus-request IA_FOCUS_REQUEST.roi_width @ 0>
    }T -1 ==
Tend

native.test.focus-run free throw
native.test.frame free-frame
