\ test the equivalence of MID and <MID> 
need simple-tester 
include "%idir%/DisplayFunction.f"    
  
0 value x1
0 value x2

: random2 ( -- )
\ set the parameters for MID
    0x10000  choose -> x1     
    0x10000 choose -> x2
;  

: test ( x1 x2 -- x)
\ take the difference between x1 and x2 and allow a rounding difference of 7
    - 7 - 0 max
;

: do-tests ( n --)
    0 do \ cr
        T{ random2 \ 9 emit x1 . x2 .
        x1 x2 <MID> x1 x2 MID test }T 0 ==
    loop
;
   
CR
Tstart 
\ corner cases  
T{ 0 65535 <MID> 0 65535 MID test }T 0 == 
T{ 32768 32768 <MID> 32768 32768 MID test }T 0 ==
T{ 65535 65535 <MID> 65535 65535 MID test }T 0 ==   
T{ 65535 0 <MID> 65535 0 MID test }T 0 ==   

\ spot checks
98 do-tests
CR
Tend

: time-MID
    ticks
    1000 0 do
        65535 0 do i 0x2000 MID drop loop
    loop
    ticks swap -
    cr ." MID took (ms) " . 
;

: time-<MID>
    ticks
    1000 0 do
        65535 0 do i 0x2000 <MID> drop loop
    loop
    ticks swap -
    cr ." <MID> took (ms) " .     
;

time-MID
time-<MID>


