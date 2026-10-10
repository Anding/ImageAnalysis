\ Full-frame display-function integration with explicit XISF output.

NEED ImageAnalysis
NEED ForthImageLoaders
NEED simple-tester

0 value img1
0 value img2
FILEPATH_SIZE allocate-buffer constant display.test.path

: display.make-frame { | frame map -- frame }
    640 480 1 allocate-frame -> frame
    frame FRAME_METADATA @ -> map
    s" 16" map =>" BITPIX"
    s" 2" map =>" NAXIS"
    s" 640" map =>" NAXIS1"
    s" 480" map =>" NAXIS2"
    s" UInt16" map =>" SMPLFRMT"
    s" Gray" map =>" COLORSPC"
    s" Light" map =>" IMAGETYP"
    s" 0" map =>" OFFSET"
    s" display-test" map =>" UUID"
    640 480 * 0 do
        i 257 * 0x10000 mod frame FRAME_BITMAP i 2* + w!
    loop
    frame
;

display.make-frame -> img1
img1 FRAME_WIDTH @ img1 FRAME_HEIGHT @ img1 FRAME_DEPTH @ allocate-frame -> img2
s" UInt16" img2 FRAME_METADATA @ =>" SMPLFRMT"
s" Gray" img2 FRAME_METADATA @ =>" COLORSPC"
s" Light" img2 FRAME_METADATA @ =>" IMAGETYP"
s" 0" img2 FRAME_METADATA @ =>" OFFSET"
s" display-test-stretched" img2 FRAME_METADATA @ =>" UUID"

img1 compute-imageStats
img1 img2 apply-displayFunction

display.test.path reset-buffer
s" E:\testdata\images\" display.test.path write-buffer drop
display.test.path buffer-punctuate-filepath
s" displayfunction1.xisf" display.test.path write-buffer drop
img2 display.test.path save-XISFimage-to

Tstart
T{ display.test.path buffer-to-string FileExists? }T -1 ==
Tend

img1 free-frame
img2 free-frame
