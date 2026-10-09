need ImageAnalysis
need simple-tester

0 value image
FILEPATH_SIZE allocate-buffer constant histogram.test.path

: histogram.test-write-filepath
    { frame suffix-addr suffix-u filepath-buffer -- }
    frame drop
    filepath-buffer reset-buffer
    s" E:\Coding\ImageAnalysis\testdata\histogram"
        filepath-buffer write-buffer drop
    suffix-addr suffix-u filepath-buffer write-buffer drop
;

ASSIGN histogram.test-write-filepath TO-DO write-filepath

: make-testXISF { | map img -- img }
    640 480 1 allocate-frame -> img
    img FRAME_METADATA @ -> map
    s" 16" map =>" BITPIX"	
    s" 2"	map =>" NAXIS"	
    s" 640" map =>" NAXIS1"
    s" 480" map =>" NAXIS2" 
    img
; 

: make-random ( image --)
    >R
    640 480 * 0 do
        0x10000 choose j ( loop obscures R@) FRAME_BITMAP i 2* + w!   \ random 16 bit words
    loop   
    R> drop
;

: make-constant ( image --)
    >R
    640 480 * 0 do
        0x8000 j ( loop obscures R@) FRAME_BITMAP i 2* + w!   \ random 16 bit words
    loop   
    0xffff R@ FRAME_BITMAP 640 480 * 1- 2* + w!
    R> drop
;   

: make-binary ( image --)
    >R
    640 480 * 0 do
        i 1 and if 0x5000 else 0xb000 then  \ alternate 0 and -1
        j ( loop obscures R@) FRAME_BITMAP i 2* + w!   \ random 16 bit words
    loop   
    R> drop
;   

make-testXISF -> image

image compute-imageStats

image make-constant
image compute-imageStats

image histogram.test.path save-Histogram

Tstart
T{ image FRAME_STATISTICS @ TOTAL_PIXELS @ }T 640 480 * ==
T{ image FRAME_STATISTICS @ MEAN @ }T 32768 ==
T{ image FRAME_STATISTICS @ MEDIAN @ }T 32768 ==
T{ image FRAME_STATISTICS @ histogram.saturated }T 1 ==
T{ s" SATPIX" image FRAME_METADATA @ >string hashS }T s" 1" hashS ==
T{ histogram.test.path buffer-to-string FileExists? }T -1 ==
Tend

image free-frame
bye