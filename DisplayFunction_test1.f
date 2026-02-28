need ImageAnalysis
include "%idir%/DisplayFunction.f"
0 value imageStats

allocate-imageStats -> imageStats

32768 imageStats MEDIAN !
4096  imageStats MEDIAN_ABSOLUTE_DEVIATION !

imageStats compute-display_parameters
cr
imageStats .display_parameters
cr

\ 0 dup cr . displayScale .
\ 8192 dup cr . displayScale .
\ 16384 dup cr . displayScale .
\ 24576 dup cr . displayScale .
\ 32768 dup cr . displayScale .
\ 40960 dup cr . displayScale .
\ 49152 dup cr . displayScale .
\ 57344 dup cr . displayScale .
\ 65535 dup cr . displayScale .
