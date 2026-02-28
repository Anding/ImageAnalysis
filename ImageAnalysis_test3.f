need imageAnalysis

0 value image
cr

cr ." load XISF text image"
s" E:\testdata\images\LUM-E155-F5100-f7843758a3f5.xisf" xisf.load-file drop -> image

cr ." compute image statistics"
image compute-imageStats
image .imageStats

cr