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
    test    ebx, ebx                \ check if byte count is zero
    jz      L$2
L$1:
    cmp     ebx, 2                  \ check if at least 2 bytes remain
    jb      L$2
    movzx   eax, word 0 [edx]       \ load 16-bit word with zero extend
    inc     dword 0 [ecx] [eax*4]   \ increment the longword at address = ax + ecx
    add     edx, 2                  \ move source pointer forward 2 bytes
    dec     ebx                     \ decrement pixel count by 1
    jmp     L$1                     \ continue loop
L$2:
    mov ebx, 08 [ebp]               \ move the 3rd stack item to the cached TOS
    lea ebp, 12 [ebp]               \ move the stack pointer up by 3 cells
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
    cmp     edx, 0xffff             \ test if edx < 0xffff, i.e. there are still bins remaining
    jb      L$1           
    or      edx, -1                 \ set edx=-1 as the return value since the target number of pixels was not reached
L$2:   
    mov     ebx, edx                \ return the median on the TOS
    lea     ebp, 4 [ebp]            \ move the stack pointer up by 1 cells
    NEXT,     
END-CODE

: compute-histogram { image imageStats -- }
\ prepare a full-resolution histogram for an image 
    image IMAGE_BITMAP
	imageStats HISTOGRAM dup 0x40000 erase
	image IMAGE_SIZE_BYTES @ 2/ dup imageStats TOTAL_PIXELS !       \ each pixel is 2 bytes
    ( bitmap histogram n ) <histogram> 
;

: compute-ASBDhistogram { image imageStats | _median -- }
\ prepare a full-resolution histogram of the absolute deviation values of image 
    imageStats MEDIAN @ -> _median
	imageStats 0x40000 HISTOGRAM_ABSOLUTE_DEVIATION erase
	image IMAGE_BITMAP dup imageStats TOTAL_PIXELS @ + swap ( end start)
	DO
		i w@ _median - abs
		4* imageStats HISTOGRAM_ABSOLUTE_DEVIATION + incr
	2 +LOOP
;  
    
: compute-mean { imageStats -- }
\ compute the mean pixel level based on the histogram
	0 0 ( cumulative_intensity as a double number)
	0x10000 0 	( end start)
	DO
		imageStats i 4* + @	    					    \ num of pixels at this value of the histogram
		( cumulative_intensity pixels) i um* d+	        \ update cumulative intensity
	LOOP	
	( cumulative_intensity) imageStats TOTAL_PIXELS @ UM/MOD nip ( mean)
	imageStats MEAN !
;

: compute-median ( imageStats -- )
    >R R@ HISTOGRAM R@ TOTAL_PIXELS @ <median>
    R> MEDIAN !
;

: compute-median_absolute_deviation ( imageStats -- )
    >R R@ HISTOGRAM_ABSOLUTE_DEVIATION R@ TOTAL_PIXELS @ <median> 
    R> MEDIAN_ABSOLUTE_DEVIATION !
;

: initialize-imageStats { image imageStats }
    image imageStats compute-histogram 
    imageStats compute-mean 
    imageStats compute-median
;
    
: histogram.saturated ( imageStats -- )
\ count the number of saturated pixels based on the histogram
	HISTOGRAM 0x3fffc + @
;

\ : histogram-down ( bits -- )
\ \ downsample the histogram by bits, 0 <= bits < 16
\ \ downsampled histogram is stored in hist.buffer2
\ 	hist.buffer2 0x40000 erase
\ 	0x10000 0	( bits end start)
\ 	DO
\ 		i over rshift 4* hist.buffer2 +
\ 		i 4* hist.buffer @
\ 		swap +!
\ 	LOOP
\ ;

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
 