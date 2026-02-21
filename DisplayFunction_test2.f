\ test the equivalence of CLIP and <CLIP> 
need simple-tester 
include "%idir%/DisplayFunction.f"    
  
0 value x1
0 value x2
0 value x3  

: choose-range ( a b -- x)
\ choose a random number a <= x < b
    over - choose +
;    

: random3 ( -- )
\ set the parameters for CLIP at random, respecting x1 < x3
    0xFFFF  choose -> x1     
    0x10000 choose -> x2
    x1 0x10000 choose-range -> x3   
;  

: do-tests ( n --)
    0 do
     T{ random3 x1 x2 x3 <CLIP> }T x1 x2 x3 CLIP ==
    loop
;
   
CR
Tstart 
\ corner cases  
T{ 0 0 65535 <CLIP> }T 0 0 65535 CLIP == 
T{ 0 65535 65535 <CLIP> }T 0 65535 65535 CLIP ==   
\ spot checks
98 do-tests
CR
Tend

: time-CLIP
    ticks
    1000 0 do
        65535 0 do i 0x2000 0xa000 CLIP drop loop
    loop
    ticks swap -
    cr ." CLIP took (ms) " . 
;

: time-<CLIP>
    ticks
    1000 0 do
        65535 0 do i 0x2000 0xa000 <CLIP> drop loop
    loop
    ticks swap -
    cr ." <CLIP> took (ms) " .     
;

time-CLIP
time-<CLIP>


