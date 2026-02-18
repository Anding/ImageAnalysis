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

CODE <histogram> ( bitmap histogram n  -- )
    \ Load pointers: EDX = source, ECX = destination. EBX contains byte count.
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
    sub     ebx, 2                  \ decrement byte counter by 2
    jmp     L$1                     \ continue loop
L$2:
    mov ebx, 08 [ebp]               \ move the 3rd stack item to the cached TOS
    lea ebp, 12 [ebp]               \ move the stack pointer up by 3 cells
    NEXT,    
END-CODE

: compute-histogram { image imageStats -- }
\ prepare a full-resolution histogram for an image 
    image IMAGE_BITMAP
	imageStats HISTOGRAM dup 0x40000 erase
	image IMAGE_SIZE_BYTES @ dup imageStats TOTAL_PIXELS !
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

: histogram.median { histbuffer total_pixels | half_total_pixels cumulative_pixels -- median }
\ compute the median pixel level based on the histogram
	total_pixels @ 2/ -> half_total_pixels 
	zero total_pixels	    				\ now iterate until half the total number of pixels have been counted
	0 begin
		dup 4* histbuffer + @				\ num of pixels at this value of the histogram
		add cumulative_pixels				\ update cumulative pixels
		cumulative_pixels half_total_pixels <
	while
		1+
	repeat	( median)
;

: compute-median ( imageStats -- )
    R> R@ HISTOGRAM R@ TOTAL_PIXELS @ histogram.median 
    R> MEDIAN !
;

: compute-median_absolute_deviation ( imageStats -- )
    R> R@ HISTOGRAM_ABSOLUTE_DEVIATION R@ TOTAL_PIXELS @ histogram.median 
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
 