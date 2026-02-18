include "%idir%\ImageAnalysis.f"

0 value image
0 value imagestats

: make-testXISF { | map img -- img }
    640 480 1 allocate-image -> img
    img FITS_MAP @ -> map
    s" 16" map =>" BITPIX"	
    s" 2"	map =>" NAXIS"	
    s" 640" map =>" NAXIS1"
    s" 480" map =>" NAXIS2" 
    640 480 * 0 do
        0x10000 choose img IMAGE_BITMAP i 2* + w!   \ random 16 bit words
    loop   
    img
;

make-testXISF -> image
allocate-imageStats -> imagestats
image imagestats compute-histogram
