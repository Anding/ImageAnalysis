need imageAnalysis
need ForthAstroFormats

cr

cr ." load XISF text image"
s" E:\testdata\images\LUM-E155-F5100-f7843758a3f5.xisf" xisf.load-file drop constant stats.example-frame

cr ." compute image statistics"
stats.example-frame compute-imageStats
stats.example-frame .imageStats

cr