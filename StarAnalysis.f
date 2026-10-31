\ Stateless native star analysis; records belong to frames.
need ForthVT100

10 value starAnalysis.star-limit \ Brightest candidates to measure, from 0 to 256.

LIBRARY: ImageAnalysisNative.dll
Extern: int "C" AnalyzeStars(
    int * pixels, int width, int height, int median, int mad,
    int max_stars, int * stars
) ;

: analyze-stars { frame | stats -- }
\ Measure the brightest candidates using the current starAnalysis.star-limit.
\ Return a flux-ordered list, with a zero terminator and no FITS output.
\ Pass existing median/MAD statistics; C owns the configured detection threshold.
    frame FRAME_STATISTICS @ -> stats
    frame FRAME_BITMAP
    frame FRAME_WIDTH @
    frame FRAME_HEIGHT @
    stats MEDIAN @
    stats MEDIAN_ABSOLUTE_DEVIATION @
    starAnalysis.star-limit
    frame FRAME_STARS @
    AnalyzeStars 0< IF s" Native star analysis failed" .E> THEN
    \ A negative result leaves an empty list; log the failure and continue.
;
