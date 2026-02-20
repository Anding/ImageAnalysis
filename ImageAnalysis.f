\ astronomical image analysis in Forth
need ForthXISF

BEGIN-STRUCTURE IMAGE_STATISTICS
    0x40000 +FIELD HISTOGRAM                        \ one 32 bit cell for each 16 bit brightness value
    0x40000 +FIELD HISTOGRAM_ABSOLUTE_DEVIATION     \ histogram of the absolute deviations of the pixels from the median
          4 +FIELD TOTAL_PIXELS
          4 +FIELD MEAN
          4 +FIELD MEDIAN
          4 +FIELD MEDIAN_ABSOLUTE_DEVIATION
END-STRUCTURE

: allocate-imageStats ( -- imageStats)
    IMAGE_STATISTICS allocate abort" unable to allocate image statistics" 
;

\ internal words in assembly language

CODE <histogram> ( bitmap histogram total_pixels  -- )
    mov     edx, 4 [ebp]            \ source pointer
    mov     ecx, 0 [ebp]            \ histogram pointer    
    test    ebx, ebx                \ exit if no pixels to process
    jz      L$2
L$1:
    movzx   eax, word 0 [edx]       \ load 16-bit word with zero extend
    inc     dword 0 [ecx] [eax*4]   \ increment the longword at address = ax + ecx
    add     edx, 2                  \ move source pointer forward 2 bytes
    dec     ebx                     \ decrement pixel count by 1
    jnz     L$1                     \ continue loop
L$2:
    mov ebx, 08 [ebp]               \ move the 2nd stack item to the cached TOS
    lea ebp, 12 [ebp]               \ move the stack pointer up by 3 cells
    NEXT,    
END-CODE

CODE <histogram-ad> ( median bitmap histogram total_pixels  -- )
    push    edi                     \ callee save
    push    esi                     \ callee save
    mov     edi, 8 [ebp]            \ median pixel value
    mov     esi, 4 [ebp]            \ address of the image array (free EDX for CDQ)
    mov     ecx, 0 [ebp]            \ pointer to the base of the histogram buffer   
    test    ebx, ebx                \ exit if no pixels to process
    jz      L$2
L$1:
    movzx   eax, word 0 [esi]       \ load 16-bit word with zero extend
    sub     eax, edi                \ diff = pixel - median; EDX free because ESI holds image pointer
    cdq                             \ EDX = 0x00000000 if EAX >= 0, 0xFFFFFFFF if EAX < 0
    xor     eax, edx                \ if negative: flip all EAX bits then add 1, which is binary negate
    sub     eax, edx                \ if positive these instructions do nothing to eax
    inc     dword 0 [ecx] [eax*4]   \ increment the histogram bin for this absolute deviation
    add     esi, 2                  \ move source pointer forward 2 bytes
    dec     ebx                     \ decrement pixel count by 1
    jnz     L$1                     \ continue loop
L$2:
    pop esi                         \ callee restore
    pop edi                         \ callee restore
    mov ebx, 12 [ebp]               \ move the 3rd stack item to the cached TOS
    lea ebp, 16 [ebp]               \ move the stack pointer up by 4 cells
    NEXT,    
END-CODE

CODE <median> ( histogram total_pixels -- m )
    mov     ecx, 0 [ebp]            \ ecx will be the address of the current histogram bin
    shr     ebx, 1                  \ divide the number of pixels by 2, this is the target we need to reach
    xor     edx, edx                \ set edx=0, edx will count through the number of histogram bins
L$1:
    sub     ebx, dword 0 [ecx]      \ subtract the number of pixels in the current bin from the remaining target
    jb      L$2                     \ edx contains the median pixel value
    inc     edx                     \ advance to the next bin
    add     ecx, 4                  \ advance to the address of the next bin
    cmp     edx, 0x10000            \ test if edx < 0x10000, i.e. there are still bins remaining
    jb      L$1           
    or      edx, -1                 \ set edx=-1 as the return value since the target number of pixels was not reached
L$2:   
    mov     ebx, edx                \ return the median on TOS
    lea     ebp, 4 [ebp]            \ move the stack pointer up by 1 cells
    NEXT,     
END-CODE

CODE <mean> ( histogram total_pixels -- m )
    push    esi                     \ callee save
    push    edi                     \ callee save
    push    ebx                     \ save the number of pixels on the stack for the final division operation
    xor     esi, esi                \ esi will be the hi 32 bits of a 64-bit accumulator
    xor     edi, edi                \ edi will be the lo 32 bits of a 64-bit accumulator
    xor     ecx, ecx                \ ecx will be the number of the current histogram bin 0..0xffff
    mov     ebx, 0 [ebp]            \ ebx will be the address of the current histogram bin

L$1:
    mov     eax, ecx                \ update eax to the number of the current histogram bin
    mul     dword 0 [ebx]           \ edx:eax = eax * [ebx], the total pixel intensity represented by this bin
    add     edi, eax                \ lo 32 bits of a 64-bit addition to the accumulator
    adc     esi, edx                \ hi 32 bits
    inc     ecx                     \ advance to the next bin
    add     ebx, 4                  \ advance to the address of the next bin
    cmp     ecx, 0x10000            \ test if ecx < 0x10000, i.e. there are still bins remaining
    jb      L$1             
    mov     edx, esi                \ move the accumulated pixel intensity to edx:eax
    mov     eax, edi
    pop     ebx                     \ reload ebx with the number of pixels
    div     ebx                     \ after the division eax contains the mean pixel value
    mov     ebx, eax                \ return the mean on TOS
    lea     ebp, 4 [ebp]            \ move the stack pointer up by 1 cells  
    pop     edi                     \ callee restore    
    pop     esi                     \ callee restore    
    NEXT,     
END-CODE

: compute-histogram { image imageStats -- }
\ prepare a full-resolution histogram for an image 
    image IMAGE_BITMAP
	imageStats HISTOGRAM dup 0x40000 erase
	image IMAGE_SIZE_BYTES @ 2/ dup imageStats TOTAL_PIXELS !       \ each pixel is 2 bytes
    ( bitmap histogram n ) <histogram> 
;

: compute-ASBDhistogram { image imageStats -- }
\ prepare a full-resolution histogram of the absolute deviation values of image 
    imageStats MEDIAN @
    image IMAGE_BITMAP   
 	imageStats HISTOGRAM_ABSOLUTE_DEVIATION dup 0x40000 erase   
 	imageStats TOTAL_PIXELS @
 	( median bitmap histogram n) <histogram-ad>
;   
    
: compute-median ( imageStats -- )
\ compute the median and update imageStats
    >R R@ HISTOGRAM R@ TOTAL_PIXELS @ <median>
    R> MEDIAN !
;

: compute-median_absolute_deviation ( imageStats -- )
\ compute the mediam absolute deviation and update imageStats
    >R R@ HISTOGRAM_ABSOLUTE_DEVIATION R@ TOTAL_PIXELS @ <median> 
    R> MEDIAN_ABSOLUTE_DEVIATION !
;

: compute-mean ( imageStats -- )
\ compute the mean and update imageStats
    >R R@ HISTOGRAM R@ TOTAL_PIXELS @ <mean>
	R> MEAN !
;

: compute-imageStats { image imageStats }
    image imageStats compute-histogram 
    imageStats compute-mean 
    imageStats compute-median
    image imagestats compute-ASBDhistogram              \ must compute the median first
    imagestats compute-median_absolute_deviation
;
    
: histogram.saturated ( imageStats -- )
\ count the number of saturated pixels based on the histogram
	HISTOGRAM 0x3fffc + @
;

: combine-images { n x y addr0 | half-n size dest -- }	\ VFX locals
\ combine n sequential x * y * 16bit monochrome images at located at addr0  
\ space for the combines image must already be allocated at the end of the set of images
	n 2/ -> half-n					\ for rounding
	x y * 2* -> size				\ size in bytes of each image
	size n * addr0 + -> dest		\ storage address of the combined image
	size 0 DO
			0 							\ cumulative count across bins
			n 0 DO
				addr0 i size * + j + w@ 
				+						\ update the cumulative
			LOOP
			half-n + n /			\ adding half-n rounds rather than truncates using integer arithmetic
			dup .
			( mean) dest i + w!
	2 +LOOP
;


\ utility functions

: .imageStats ( imageStats --)
    cr ." Mean      " dup mean ?
    cr ." Median    " dup median ? 
    cr ." Median AD " dup median_absolute_deviation ?
    drop
;
