\ Export an image histogram through the active binary filepath policy.

NEED ForthPublication

: save-Histogram { frame filepath-buffer -- }
\ Histograms contain 65536 four-byte bins, hence the fixed 0x40000-byte file.
    frame frame FRAME_STATISTICS @ HISTOGRAM 0x40000 filepath-buffer
        save-binary-file
;