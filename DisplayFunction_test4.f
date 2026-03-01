need ImageAnalysis

0 value img1
0 value img2

: test_write-XISFfilepath { map buf -- }
	s" E:\testdata\images\" buf write-buffer drop
	buf buffer-punctuate-filepath
	s" displayfunction1.xisf" buf write-buffer drop 
;

cr cr ." load XISF text image"
s" E:\testdata\images\LUM-E14-F5100-900080d4354b.xisf" xisf.load-file drop -> img1
img1 xisf.spawn -> img2

cr ." compute image statistics and display parameters"
img1 compute-imageStats

cr ." run display function... "
ticks
img1 img2 apply-displayFunction    \ 141ms in pure assembly; 516 ms in all Forth, 218 ms with assembly subroutines in a Forth loop
ticks swap -
. ."  ms "               

img1 .imageStats

ASSIGN test_write-XISFfilepath TO-DO write-XISFfilepath
img2 save-XISFimage
