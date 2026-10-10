need ImageAnalysis
need simple-tester
0 value imageStats

allocate-frameStats -> imageStats

32768 imageStats MEDIAN !
4096  imageStats MEDIAN_ABSOLUTE_DEVIATION !

imageStats compute-displayParameters

Tstart
T{ imageStats df.MADN @ }T 6072 ==
T{ imageStats df.B @ }T 16384 ==
T{ imageStats df.s @ }T 15767 ==
T{ imageStats df.h @ }T 65536 ==
Tend

imageStats free drop
