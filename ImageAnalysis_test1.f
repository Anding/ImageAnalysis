\ unit test for ImageAnalyis.f

include %idir%\ImageAnalysis.f

: report
	base @ hex
	CR
	hist.buffer 64 dump	
	CR
	." saturated " image.saturated .

	CR
	." mean      " image.mean .

	CR
	." median    " image.median .
	base !
;

0 value image1
0 value image2
0 value image3
0 value image4


here -> image1
	 0x0001 w,  0x0002 w,  0x0003 w,  0x0004 w,
	 0x0005 w,  0x0006 w,  0x0007 w,  0x0008 w,
	 0x0009 w,  0x000a w,  0x000b w,  0x000c w,
	 0x000d w,  0x000e w,  0x000f w,  0x0010 w,

4 4 image1 make-histogram report

here -> image2
	 0x0001 w,  0x0002 w,  0x0003 w,  0x0004 w,
	 0x0005 w,  0x0006 w,  0x0007 w,  0x0008 w,
	 0x0009 w,  0x000a w,  0x000b w,  0x000c w,
	 0x000d w,  0x000e w,  0x000f w,  0xffff w,

4 4 image2 make-histogram report

here -> image3
	 0xffff w,  0x0002 w,  0x0003 w,  0x0004 w,
	 0x0005 w,  0x0006 w,  0x0007 w,  0x0008 w,
	 0x0009 w,  0x000a w,  0x000b w,  0x000c w,
	 0x000d w,  0x000e w,  0x000f w,  0x0010 w,

4 4 image3 make-histogram report

here -> image4
	0x20 allot

3 4 4 image1 combine-images

CR image4 0x20 dump
