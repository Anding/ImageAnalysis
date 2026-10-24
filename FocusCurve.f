\ Numerical focus-curve binding; equipment and acceptance policy belong to callers.
LIBRARY: ImageAnalysisNative.dll

\ ( max_stars -- status ) Reset the curve run and select identities to follow.
Extern: int "C" start_computation( int max_stars ) ;
\ ( stars position exposure-addr -- status ) Copy a sample or recommend exposure.
Extern: int "C" add_focus_point(
    int * stars, int focus_position, int * exposure
) ;
\ ( focus-addr error-addr -- status ) Finish the run; outputs are integer steps.
Extern: int "C" finish_computation( int * best_focus, int * error ) ;
