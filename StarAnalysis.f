\ Stateless native star analysis; records and storage belong to the frame.
need ImageAnalysis
need ForthVT100

LIBRARY: ImageAnalysisNative.dll
Extern: int "C" AnalyzeStars(
    int * pixels, int width, int height, int * histogram, int * stars
) ;

: analyze-stars { frame -- }
\ Replace FRAME_STARS with up to 256 records and a zero terminator; no FITS.
\ Caller must have prepared the histogram; both buffers belong to the frame.
    frame FRAME_BITMAP
    frame FRAME_WIDTH @
    frame FRAME_HEIGHT @
    frame FRAME_STATISTICS @ HISTOGRAM
    frame FRAME_STARS @
    AnalyzeStars 0< IF s" Native star analysis failed" .E> THEN
    \ A negative result leaves an empty list; log the failure and continue.
;
