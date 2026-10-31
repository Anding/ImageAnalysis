\ Stateless native star analysis; records belong to frames.
need ForthVT100

LIBRARY: ImageAnalysisNative.dll
Extern: int "C" AnalyzeStars(
    int * pixels, int width, int height, int median, int mad,
    int * stars
) ;

: analyze-stars { frame | stats -- }
\ Measure up to 256 candidates with shape and exposure-risk metrics.
\ Return a flux-ordered list, with a zero terminator and no FITS output.
\ Pass existing median/MAD statistics; C owns the configured detection threshold.
    frame FRAME_STATISTICS @ -> stats
    frame FRAME_BITMAP
    frame FRAME_WIDTH @
    frame FRAME_HEIGHT @
    stats MEDIAN @
    stats MEDIAN_ABSOLUTE_DEVIATION @
    frame FRAME_STARS @
    AnalyzeStars 0< IF s" Native star analysis failed" .E> THEN
    \ A negative result leaves an empty list; log the failure and continue.
;
