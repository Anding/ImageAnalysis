\ Export a histogram at an explicit destination.

: SaveHistogramAsBinary { histogram zaddr | fileid -- IOR }
    zaddr zcount delete-file drop
    zaddr zcount w/o create-file if -1 exit then -> fileid
    histogram 0x40000 fileid write-file drop
    fileid close-file ( IOR)
;

: save-Histogram-to { frame filepath-buffer -- }
    filepath-buffer create-imageDirectory
    frame FRAME_STATISTICS @ HISTOGRAM
    filepath-buffer buffer-to-string drop
    ( bitmap width height caddr) SaveHistogramAsBinary abort" Error writing histogram file"
;