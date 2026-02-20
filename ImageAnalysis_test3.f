need imageAnalysis

0 value image
0 value imagestats
cr

cr ." load XISF text image"
s" E:\testdata\images\LUM-E155-F5100-f7843758a3f5.xisf" xisf.load-file drop -> image

cr ." compute image statistics"
allocate-imageStats -> imagestats
image imagestats compute-imageStats
imagestats .imageStats

cr