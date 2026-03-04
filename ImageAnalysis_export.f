\ export histogram as a binary file

DEFER write-HISTfilepath ( map buf --)

: default_write-HISTfilepath { map buf -- }
	s" e:\images\snapshot\" buf write-buffer drop	
	buf buffer-punctuate-filepath
	s" hist.raw" buf write-buffer drop
	0 buf echo-buffer drop                                   \ zero terminated string
;

    ASSIGN default_write-HISTfilepath TO-DO write-HISTfilepath
    
: initialize-HISTfilepath ( img --)
	>R
	R@ FITS_MAP @ ( map)
	R> RAW_FILEPATH_BUFFER                                   \ reuse the RAW filpath
	FILEPATH_SIZE over ( map buf FILEPATH_SIZE buf) declare-buffer
	( map buf) write-HISTfilepath
;

: SaveHistogramAsBinary { histogram zaddr | fileid -- IOR }
    zaddr zcount delete-file drop
    zaddr zcount w/o create-file if -1 exit then -> fileid
    histogram 0x40000 fileid write-file drop
    fileid close-file ( IOR)
;

: save-Histogram { img -- }
	img initialize-HISTfilepath
	img RAW_FILEPATH_BUFFER create-imageDirectory
    img IMAGE_STATISTICS @ HISTOGRAM
    img RAW_FILEPATH_BUFFER buffer-to-string drop
    ( bitmap width height caddr) SaveHistogramAsBinary abort" Error writing histogram file"
;