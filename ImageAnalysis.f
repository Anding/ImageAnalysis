\ astronomical image analysis in Forth

: histogram ( x y addr -- hist)
\ prepare a full-resolution histogram for an x * y * 16bit monochrome image at addr
;

: hist.mean ( hist -- x)
\ compute the mean pixel level based on the histogram
;

: hist.median ( hist -- x)
\ compute the median pixel level based on the histogram
;

: hist.saturated ( hist -- x)
\ count the number of saturated pixels based on the histogram
;

: histogram-down ( bits hist -- hist-lores)
\ downsample the histogram, 0 < bits < 16
;

: combine-images ( mode n x y addr -- addr')
\ combine n sequential x * y * 16bit monochrome images at located addr 
\ space (at addr') for the combines image must already be allocated at the end of the set of images
\ mode 0 = mean, other modes reserved for expansion
;
 