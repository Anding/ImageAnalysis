# ImageAnalysis
Astronomical image analysis in Forth

`compute-imageStats` builds the full 16-bit histogram and publishes integer
`MEAN`, `MEDIAN`, `MEDIANAD`, and `SATPIX` values into the frame's ordered
FITS metadata map. `SATPIX` is the exact number of pixels in histogram bin
65535; a saturation fraction remains derivable from the image dimensions.
