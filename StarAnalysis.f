\ Stateless native star analysis; records and storage belong to the frame.
need ImageAnalysis
need ForthVT100

LIBRARY: ImageAnalysisNative.dll
Extern: int "C" AnalyzeStars(
    int * pixels, int width, int height, int background,
    int max_stars, int * stars
) ;

: star-background { frame | low high rest remaining -- background }
\ Read integer MEDIAN metadata; missing, malformed or out-of-range means zero.
    0 0 s" MEDIAN" frame FRAME_METADATA @ map@ >number
    -> remaining -> rest -> high -> low
    remaining high or 0= low 65536 u< and if low else 0 then
;

: analyze-stars { frame max_stars -- }
\ Measure the requested number of brightest candidates, supplied on each call.
\ Return a flux-ordered list, with a zero terminator and no FITS output.
\ Background comes from MEDIAN metadata, not a histogram calculation.
    frame FRAME_BITMAP
    frame FRAME_WIDTH @
    frame FRAME_HEIGHT @
    frame star-background
    max_stars
    frame FRAME_STARS @
    AnalyzeStars 0< IF s" Native star analysis failed" .E> THEN
    \ A negative result leaves an empty list; log the failure and continue.
;
