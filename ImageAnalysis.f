\ astronomical image analysis in Forth

0x40000 buffer: hist.buffer
0x40000 buffer: hist.buffer2
\ one cell for each 16 bit brightness value

: make-histogram ( x y addr --)
\ prepare a full-resolution histogram for an x * y * 16bit monochrome image at addr
	hist.buffer 0x40000 erase
	>R * 2* R@ + R> ( end start)
	DO
		i w@ 4* hist.buffer + incr
	2 +LOOP
;

: histogram.mean ( -- mean) { | cumulative_pixels }		\ VFX locals
\ compute the mean pixel level based on the histogram
	0 0 ( cumulative_intensity)
	zero cumulative_pixels
	0x10000 0 	( end start)
	DO
		i 4* hist.buffer + @								\ num of pixels at this value of the histogram
		dup add cumulative_pixels						\ update cumulative pixels; expect x * y at the end
		( cumulative_intensity pixels) i um* d+	\ update cumulative intensity
	LOOP	
	( cumulative_intensity) cumulative_pixels .s  UM/MOD nip ( mean)
;

: histogram.median ( -- x) { | cumulative_pixels half_total_pixels }		\ VFX locals
\ compute the median pixel level based on the histogram
	zero cumulative_pixels					\ first find the total number of pixels in the image
	0x10000 0	( end start)
	DO
		i 4* hist.buffer + @					\ num of pixels at this value of the histogram
		add cumulative_pixels				\ update cumulative pixels; expect x * y at the end
	LOOP	
	cumulative_pixels 2/ -> half_total_pixels 

	zero cumulative_pixels					\ now iterate until half the total number of pixels have been counted
	0 begin
		dup 4* hist.buffer + @				\ num of pixels at this value of the histogram
		add cumulative_pixels				\ update cumulative pixels
		cumulative_pixels half_total_pixels <
	while
		1+
	repeat	( median)
;

: histogram.saturated ( hist -- x)
\ count the number of saturated pixels based on the histogram
	hist.buffer 0x3fffc + @
;

: histogram-down ( bits -- )
\ downsample the histogram by bits, 0 <= bits < 16
\ downsampled histogram is stored in hist.buffer2
	hist.buffer2 0x40000 erase
	0x10000 0	( bits end start)
	DO
		i over rshift 4* hist.buffer2 +
		i 4* hist.buffer @
		swap +!
	LOOP
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
 