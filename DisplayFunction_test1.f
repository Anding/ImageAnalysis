\ Unit regression for the display-function parameter calculation.
\ Maintenance: these exact integer results protect the Forth fixed-point
\ contract independently of image loading and output codecs, including the
\ bright-branch stats-field loads and zero-signal endpoint protection.

need ImageAnalysis
need simple-tester
0 value imageStats

allocate-frameStats -> imageStats

32768 imageStats MEDIAN !
4096  imageStats MEDIAN_ABSOLUTE_DEVIATION !

imageStats compute-displayParameters

: display.test-flat { value | frame -- first last }
    32 16 1 allocate-frame -> frame
    512 0 do value frame FRAME_BITMAP i 2* + w! loop
    frame compute-imageStats
    frame frame FRAME_BITMAP apply-displayFunction-to
    frame FRAME_BITMAP w@
    frame FRAME_BITMAP 511 2* + w@
    frame free-frame
;

Tstart
T{ imageStats df.MADN @ }T 6072 ==
T{ imageStats df.B @ }T 16384 ==
T{ imageStats df.s @ }T 15767 ==
T{ imageStats df.h @ }T 65536 ==

49152 imageStats MEDIAN !
T{ imageStats compute-displayParameters }T ==
T{ imageStats df.a @ }T 65536 ==
T{ imageStats df.s @ }T 0 ==
T{ imageStats df.h @ }T 65536 ==
T{ imageStats df.t @ }T 16384 ==
T{ imageStats df.m @ }T 32768 ==

0 imageStats MEDIAN_ABSOLUTE_DEVIATION !
T{ imageStats compute-displayParameters }T ==
T{ imageStats df.h @ }T 65536 ==
T{ imageStats df.t @ }T 16384 ==
T{ imageStats df.m @ }T 32768 ==
T{ 0 0 <MID> }T 0 ==
T{ 32768 0 <MID> }T 0 ==
T{ 16384 16384 <MID> }T 32768 ==
T{ 16384 65536 <MID> }T 65536 ==
T{ 0 display.test-flat }T 0 0 ==
T{ 49152 display.test-flat }T 49152 49152 ==
Tend

imageStats free drop
