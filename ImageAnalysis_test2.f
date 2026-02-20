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
    img
; 

: make-random ( image --)
    >R
    640 480 * 0 do
        0x10000 choose j ( loop obscures R@) IMAGE_BITMAP i 2* + w!   \ random 16 bit words
    loop   
    R> drop
;

: make-constant ( image --)
    >R
    640 480 * 0 do
        0x8000 j ( loop obscures R@) IMAGE_BITMAP i 2* + w!   \ random 16 bit words
    loop   
    R> drop
;   

: make-binary ( image --)
    >R
    640 480 * 0 do
        i 1 and if 0x5000 else 0xb000 then  \ alternate 0 and -1
        j ( loop obscures R@) IMAGE_BITMAP i 2* + w!   \ random 16 bit words
    loop   
    R> drop
;   

make-testXISF -> image
allocate-imageStats -> imagestats

cr
cr ." zero image"
image imagestats compute-imageStats
imagestats .imageStats
cr

cr ." constant image"
image make-constant
image imagestats compute-imageStats
imagestats .imageStats
cr

cr ." binary image"
image make-binary
image imagestats compute-imageStats
imagestats .imageStats
cr

cr ." random image"
image make-random
image imagestats compute-imageStats
imagestats .imageStats
cr