\ XISF display function implementation

\ Adopt "integer format" where 1.0 is represented as 0x10000
\ This is convenient for adapting algorithms intended to handle pixel values in the range 0..1
\ that are repesented as 16 bit unsigned number in the range 0..65535

65536.0E fconstant scaleFactor
65536   constant 1.00I
131072  constant 2.00I
32768   constant 0.50I
16384   constant 0.25I

: FtoI ( f -- i)
\ convert a floating point number to integer format
    scaleFactor f* f>s
;

: ItoF ( i -- f)
\ convert an integer format number to floating point
    s>f scaleFactor f/
;

: I* ( i1 i2 -- i3)
\ multiply two integer format numbers, i3 = i1 * i2
    1.00I */  
;
  
: I/ ( i1 i2 -- i3)
\ divide two integer format numbers, i3 = i1 / i2
    1.00I swap */
;

: I+ ( i1 i2 -- i3)
\ add two integer format numbers, i3 = i1 + i2
    +
;

: I- ( i1 i2 -- i3)
\ subtract two integer format numbers, i3 = i1 - i2
    -
;  

\ *******************************************************
\ XISF display functions, see DisplayFunction.md 

: CLIP { x s h -- xc }
    x s < if 
        0
    else 
        x h > if 
            1.00I
        else
            x s I- 
            h s I-
            I/
         then
     then
;

CODE <CLIP> ( x s h -- xc )
\ assembly language version of CLIP
    push    esi                 \ callee save
    push    edi                 \ callee save
                                \ ebx = h
    mov     esi, 0 [ebp]        \ esi = s
    mov     edi, 4 [ebp]        \ edi = x
    cmp     edi, esi            \ fall through if x < s
    jae     L$1
    xor     ebx, ebx            \ ebx = 0.0
    jmp     L$3
L$1:                            
    cmp     edi, ebx            \ fall through if x >= h
    jb      L$2                   
    mov     ebx, 0x10000        \ ebx = 1.0
    jmp     L$3                
L$2:                             
    mov     eax, edi            \ eax = x
    sub     eax, esi            \ eax = x-s
    shl     eax, 16             \ eax = (x-s)*0x10000  (left shift by 16 = multiply by 0x10000)
                                \ valid since x-s fits in 16 bits (s <= x < h <= 0xffff)
    xor     edx, edx            \ clear EDX: div uses EDX:EAX as 64-bit dividend
    sub     ebx, esi            \ ebx = h-s
    div     ebx                 \ eax = (x-s)/(h-s)
    mov     ebx, eax            \ ebx = (x-s)/(h-s)
L$3:                         
    pop     edi                 \ callee restore
    pop     esi                 \ callee restore
    lea ebp, 08 [ebp]           \ move the stack pointer up by 2 cells
    NEXT,    
END-CODE

: MID { x m -- xm }
    x 0= if 
        0
    else
        x m = if
            0.50I
        else
            x 1.00I = if
                1.00I
            else
                m 1.00I I- x I*
                m 2.00I I* 1.00I I-
                x I*
                m I-
                I/
             then
         then
     then
;

1.4826E    FtoI constant  1.48I
-2.80E     FtoI constant -2.80I
0.25I   value df.B
-2.80I  value df.C
0       value df.MADN
0       value df.a
0       value df.t
0       value df.s
0       value df.h
0       value df.m

: compute-display_parameters { M MAD | -- }
    MAD 1.48I I* -> df.MADN
    M 0.50I > if 1.00I else 0 then -> df.a
    
    df.a 1.00I = if 
        0 
    else
        df.MADN 0= if
            0
        else
            M df.C df.MADN I* I+
            0 max
            1.00I min
         then
     then -> df.s
     
     df.a 0= if
        1.00I
     else
        df.MADN 0= if
            1.00I
        else
            M df.C df.MADN I* I-
            0 max
            1.00I min
        then
     then -> df.h
     
     df.a 0= if
        M df.s I- -> df.t
        df.t df.B MID
     else
        df.h M I- -> df.t
        df.B df.t MID
     then -> df.m   
;

: displayScale ( x -- x1)
    ( x) df.s df.h <CLIP> ( xc)
    ( xc) df.m MID ( xm)
;

\ ****************************************************
\ utility functions

: .display_parameters
    ." MADN " df.MADN . cr 
    ." B    " df.B . cr
    ." C    " df.C . cr
    ." a    " df.a . cr
    ." s    " df.s . cr
    ." h    " df.h . cr
    ." t    " df.t . cr
    ." m    " df.m . cr
;
    
    