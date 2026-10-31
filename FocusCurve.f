\ Numerical focus-curve binding; equipment and acceptance policy belong to callers.
LIBRARY: ImageAnalysisNative.dll

\ ( max_stars -- status ) Reset the run; limit final survivors, not tracked identities.
Extern: int "C" start_computation( int max_stars ) ;
\ ( stars exposure-addr -- status ) Select the fixed pool and recommend run baseline.
Extern: int "C" select_focus_stars( int * stars, int * exposure ) ;
\ ( stars position exposure-addr -- status ) Retain best trial; recommend local exposure.
Extern: int "C" add_focus_point(
    int * stars, int focus_position, int * exposure
) ;
\ ( -- status ) Commit the best available trial at this focus position.
Extern: int "C" finish_focus_point( ) ;
\ ( focus-addr error-addr -- status ) Finish the run; outputs are integer steps.
Extern: int "C" finish_computation( int * best_focus, int * error ) ;
