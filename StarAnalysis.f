\ Stateless native star analysis; records and storage belong to the frame.
need ImageAnalysis
need ForthVT100

LIBRARY: ImageAnalysisNative.dll
Extern: int "C" AnalyzeStars(
    int * pixels, int width, int height, int * histogram,
    int max_stars, int * stars
) ;

10 value default.maxstars \ Environment-selectable measurement limit, 0..256.

: analyze-stars { frame -- }
\ Read default.maxstars on each call; measure that many brightest candidates.
\ Return a flux-ordered list, with a zero terminator and no FITS output.
\ Caller must have prepared the histogram; both buffers belong to the frame.
    frame FRAME_BITMAP
    frame FRAME_WIDTH @
    frame FRAME_HEIGHT @
    frame FRAME_STATISTICS @ HISTOGRAM
    default.maxstars
    frame FRAME_STARS @
    AnalyzeStars 0< IF s" Native star analysis failed" .E> THEN
    \ A negative result leaves an empty list; log the failure and continue.
;
