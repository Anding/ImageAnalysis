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

: CLIP { s h x -- xc }
    x s < if 
        0
    else 
        x h > if 
            1.00I 1- \ ~1.0
        else
            x s I- 
            h s I-
            I/
         then
     then
;

LBL: <<CLIP>> ( ebx = x, esi = h edi = s, return in ebx)
\ assembly language subroutine for CLIP
    cmp     ebx, edi            \ fall through if x < s
    jae     L$1
    xor     ebx, ebx            \ ebx = 0.0
    jmp     L$3
L$1:                            
    cmp     ebx, esi            \ fall through if x > h
    jbe      L$2                   
    mov     ebx, 0xffff         \ ebx ~ 1.0
    jmp     L$3                
L$2:                             
    mov     eax, ebx            \ eax = x
    sub     eax, edi            \ eax = x-s
    shl     eax, 16             \ eax = (x-s)*0x10000  (left shift by 16 = multiply by 0x10000)
                                \ valid since x-s fits in 16 bits (s <= x < h <= 0xffff)
    xor     edx, edx            \ clear EDX: div uses EDX:EAX as 64-bit dividend
    sub     esi, edi            \ ebx = h-s
    div     esi                 \ eax = (x-s)/(h-s)
    mov     ebx, eax            \ ebx = (x-s)/(h-s)
L$3:     
    ret
END-CODE

CODE <CLIP> ( s h x -- xc )
\ assembly language version of CLIP
    push    esi                 \ callee save
    push    edi                 \ callee save
                                \ ebx = x
    mov     esi, 0 [ebp]        \ esi = h
    mov     edi, 4 [ebp]        \ edi = s
    call     <<CLIP>>                 
    pop     edi                 \ callee restore
    pop     esi                 \ callee restore
    lea     ebp, 08 [ebp]       \ move the stack pointer up by 2 cells, return value in ebx
    NEXT,    
END-CODE

: MID { m x -- xm }
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

LBL: <<MID>> ( ebx = x, edi = m, return in ebx)
    \ compute the denominator
    mov     ecx, edi            \ ecx = m
    shl     ecx, 1              \ ecx = 2*m
    sub     ecx, 0x10000        \ ecx = 2*m - 1.0, recall 1.0 = 0x10000
    mov     eax, ebx            \ eax = x
    imul    ecx                 \ edx:eax = (2m-1) * x  (signed: 2m-1 is negative when m < 0.5)
    shrd    eax, edx, 16        \ remove the scale factor
    sub     eax, edi            \ eax = (2m-1)*x - m
    mov     ecx, eax            \ ecx = (2m-1)*x - m
    \ compute the numerator     \
    mov     esi, edi            \ esi = m
    sub     esi, 0x10000        \ esi = m - 0x10000  (m-1 in fixed point; negative since m < 1)
    mov     eax, ebx            \ eax = x
    imul    esi                 \ edx:eax = (m-1)*x
    \ shr/shl by 0x10000 cancel: keep raw product as numerator for the division
    \ perform the division
    idiv    ecx                 \ eax = (m-1)*x / (2m-1)*x - m
    mov     ebx, eax            \ ebx = result          
    ret
END-CODE

CODE <MID> ( m x -- xm )
\ assembly language version of MID
    push    esi                 \ callee save
    push    edi                 \ callee save
                                \ ebx = x
    mov     edi, 0 [ebp]        \ edi = m
    call     <<MID>>
    pop     edi                 \ callee restore
    pop     esi                 \ callee restore
    lea     ebp, 04 [ebp]       \ move the stack pointer up by 1 cell, return value in ebx
    NEXT,    
END-CODE

1.4826E    FtoI constant  1.48I
-2.80E     FtoI constant -2.80I

: compute-displayParameters { ist -- }
\ take an imageStats buffer and complete the display function parameters
    0.25I ist df.B !
    -2.80I ist df.C !
    ist MEDIAN_ABSOLUTE_DEVIATION @ 1.48I I* ist df.MADN !
    ist MEDIAN @ 0.50I > if 1.00I else 0 then ist df.a !
    
    ist df.a @ 1.00I = if 
        0 
    else
        ist df.MADN @ 0= if
            0
        else
            ist MEDIAN @ ist df.C @ ist df.MADN @ I* I+
            0 max
            1.00I min
         then
     then ist df.s !
     
     ist df.a @ 0= if
        1.00I
     else
        ist df.MADN @ 0= if
            1.00I
        else
            ist MEDIAN @ ist df.C @ ist df.MADN @ I* I-
            0 max
            1.00I min
        then
     then ist df.h !
     
     ist df.a @ 0= if
        ist MEDIAN @ ist df.s @ I- ist df.t !
        ist df.B @ ist df.t @ MID
     else
        ist df.h MEDIAN @ I- ist df.t !
        ist df.t @ ist df.B @ MID
     then ist df.m !  
;

CODE <apply-displayFunction> ( s h m src dest pixels --)
    push    esi
    push    edi
    mov     ecx, ebx                \ ecx contains the pixel count
    mov     edx,  0 [ebp]           \ edx contains the dest
    mov     eax,  4 [ebp]           \ eax contains the src
L$1:
    test    ecx, ecx                \ check if byte count is zero
    jz      L$2
    movzx   ebx, word 0 [eax]       \ ebx contains x0
    push    edx
    push    eax
    push    ecx    
    mov     esi, 12 [ebp]           \ esi contains h
    mov     edi, 16 [ebp]           \ edi contains s 
    call    <<CLIP>>                \ ebx contains x1
    mov     edi,  8 [ebp]           \ edi contains m
    call    <<MID>>                 \ ebx contains x2
    pop     ecx
    pop     eax
    pop     edx
    mov     word 0 [edx], ebx   
    dec     ecx
    add     edx, 2
    add     eax, 2  
    jmp     L$1
L$2:
    pop     edi
    pop     esi
    mov     ebx, 20 [ebp]           \ move the below stack item to TOS register since there is no return value
    lea     ebp, 24 [ebp]           \ move the stack pointer up by 6 cells
    NEXT,
END-CODE

: apply-displayFunction { imgSrc imgDst -- }
\ apply the display function
\ assumes that image stistics are already computed for imgSrc and that imgDst is already allocated
    imgSrc IMAGE_STATISTICS @ compute-displayParameters
    imgSrc IMAGE_STATISTICS @ df.s @ 
    imgSrc IMAGE_STATISTICS @ df.h @
    imgSrc IMAGE_STATISTICS @ df.m @
    imgSrc IMAGE_BITMAP
    imgDst IMAGE_BITMAP
    imgSrc IMAGE_STATISTICS @ TOTAL_PIXELS @
    ( s h m src dest pixels --) <apply-displayFunction>
;

\ ****************************************************
\ utility functions

: .display_parameters ( imageStats)
    ." MADN " dup df.MADN @ . cr 
    ." B    " dup df.B @ . cr
    ." C    " dup df.C @ . cr
    ." a    " dup df.a @ . cr
    ." s    " dup df.s @ . cr
    ." h    " dup df.h @ . cr
    ." t    " dup df.t @ . cr
    ." m    " dup df.m @ . cr  
    drop
;


    
    