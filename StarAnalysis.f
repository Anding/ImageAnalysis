\ Stateless native star analysis; records and storage belong to the frame.
need ImageAnalysis

LIBRARY: ImageAnalysisNative.dll
Extern: int "C" AnalyzeStars(
    int * pixels, int width, int height, int * histogram, int * stars
) ;

: analyze-stars { frame -- }
\ Replace FRAME_STARS with up to 256 records and a zero terminator; no FITS.
\ Prepare this frame's histogram independently of other pipeline requests.
    frame FRAME_STATISTICS @ 0= if
        allocate-frameStats frame FRAME_STATISTICS !
    then
    frame compute-histogram
    frame FRAME_BITMAP
    frame FRAME_WIDTH @
    frame FRAME_HEIGHT @
    frame FRAME_STATISTICS @ HISTOGRAM
    frame FRAME_STARS @
    AnalyzeStars dup 0< abort" Native star analysis failed" drop
;
