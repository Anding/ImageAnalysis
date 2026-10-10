\ astronomical image analysis in Forth
need ForthAstroFormats

LIBRARY: ImageAnalysisNative.dll
Extern: int "C" IA_Version() ;
Extern: void "C" IA_DefaultConfig( int * config ) ;
Extern: int "C" IA_WorkspaceBytes( int maximum_stars ) ;
Extern: int "C" IA_AnalyzeFrame(
    int * pixels, int width, int height, int stride_pixels,
    int * config, int * workspace, int workspace_bytes,
    int * stars, int star_capacity, int * summary
) ;
Extern: int "C" IA_SummarizeStars(
    int * stars, int star_count, int reject_flags,
    int * workspace, int workspace_bytes, int * summary
) ;
Extern: int "C" IA_MeasureLockedStars(
    int * pixels, int width, int height, int stride_pixels,
    int * config, int * locked_stars, int star_count,
    int search_radius_pixels, int * measurements
) ;
Extern: int "C" IA_RecommendExposure(
    int current_milliseconds, int measured_peak_adu, int background_adu,
    int target_peak_adu, int minimum_milliseconds, int maximum_milliseconds,
    int maximum_change_milli, int * recommended_milliseconds
) ;
Extern: int "C" IA_ComputeStarROI(
    int * stars, int star_count, int image_width, int image_height,
    int margin_pixels, int * roi
) ;
Extern: int "C" IA_FitFocusModel(
    int * model, int * samples, int sample_count, int reject_sample_flags,
    int resolution_milli_steps, int * fit
) ;
Extern: int "C" IA_CombineFocusFits(
    int * fits, int fit_count, int reject_fit_flags,
    int * workspace, int workspace_bytes, int * summary
) ;

0  constant IA.OK
-3 constant IA.NO_STARS

1  constant IA.STAR.SATURATED
2  constant IA.STAR.EDGE
4  constant IA.STAR.BLENDED
8  constant IA.STAR.TOO_SMALL
16 constant IA.STAR.LOW_SNR

BEGIN-STRUCTURE <IA_CONFIG>
    4 +FIELD IA_CONFIG.struct_size
    4 +FIELD IA_CONFIG.detection_sigma_milli
    4 +FIELD IA_CONFIG.minimum_threshold_adu
    4 +FIELD IA_CONFIG.saturation_adu
    4 +FIELD IA_CONFIG.minimum_area_pixels
    4 +FIELD IA_CONFIG.minimum_separation_pixels
    4 +FIELD IA_CONFIG.aperture_radius_pixels
    4 +FIELD IA_CONFIG.annulus_inner_radius_pixels
    4 +FIELD IA_CONFIG.annulus_outer_radius_pixels
    4 +FIELD IA_CONFIG.minimum_snr_milli
    4 +FIELD IA_CONFIG.blend_minimum_contrast_milli
END-STRUCTURE

BEGIN-STRUCTURE <IA_STAR>
    4 +FIELD IA_STAR.struct_size
    4 +FIELD IA_STAR.x_milli
    4 +FIELD IA_STAR.y_milli
    4 +FIELD IA_STAR.flux_adu
    4 +FIELD IA_STAR.peak_adu
    4 +FIELD IA_STAR.saturated_pixels
    4 +FIELD IA_STAR.background_milli_adu
    4 +FIELD IA_STAR.noise_milli_adu
    4 +FIELD IA_STAR.snr_milli
    4 +FIELD IA_STAR.hfr_milli_pixels
    4 +FIELD IA_STAR.hfd_milli_pixels
    4 +FIELD IA_STAR.fwhm_milli_pixels
    4 +FIELD IA_STAR.sigma_major_milli_pixels
    4 +FIELD IA_STAR.sigma_minor_milli_pixels
    4 +FIELD IA_STAR.ellipticity_milli
    4 +FIELD IA_STAR.theta_millidegrees
    4 +FIELD IA_STAR.positive_pixels
    4 +FIELD IA_STAR.left
    4 +FIELD IA_STAR.top
    4 +FIELD IA_STAR.right
    4 +FIELD IA_STAR.bottom
    4 +FIELD IA_STAR.flags
END-STRUCTURE

BEGIN-STRUCTURE <IA_FRAME_SUMMARY>
    4 +FIELD IA_FRAME_SUMMARY.struct_size
    4 +FIELD IA_FRAME_SUMMARY.detected_stars
    4 +FIELD IA_FRAME_SUMMARY.measured_stars
    4 +FIELD IA_FRAME_SUMMARY.included_stars
    4 +FIELD IA_FRAME_SUMMARY.saturated_stars
    4 +FIELD IA_FRAME_SUMMARY.edge_stars
    4 +FIELD IA_FRAME_SUMMARY.blended_stars
    4 +FIELD IA_FRAME_SUMMARY.low_snr_stars
    4 +FIELD IA_FRAME_SUMMARY.background_milli_adu
    4 +FIELD IA_FRAME_SUMMARY.noise_milli_adu
    4 +FIELD IA_FRAME_SUMMARY.detection_threshold_adu
    4 +FIELD IA_FRAME_SUMMARY.median_hfd_milli_pixels
    4 +FIELD IA_FRAME_SUMMARY.mad_hfd_milli_pixels
    4 +FIELD IA_FRAME_SUMMARY.median_fwhm_milli_pixels
    4 +FIELD IA_FRAME_SUMMARY.mad_fwhm_milli_pixels
    4 +FIELD IA_FRAME_SUMMARY.median_ellipticity_milli
    4 +FIELD IA_FRAME_SUMMARY.median_snr_milli
END-STRUCTURE

BEGIN-STRUCTURE <IA_RECT>
    4 +FIELD IA_RECT.x
    4 +FIELD IA_RECT.y
    4 +FIELD IA_RECT.width
    4 +FIELD IA_RECT.height
END-STRUCTURE

1 constant IA.FOCUS.PIECEWISE_V
2 constant IA.FOCUS.HYPERBOLA

1 constant IA.FOCUS.INSUFFICIENT_POINTS
2 constant IA.FOCUS.ONE_SIDED
4 constant IA.FOCUS.AT_BOUNDARY
8 constant IA.FOCUS.INVALID_MODEL

BEGIN-STRUCTURE <IA_FOCUS_MODEL>
    4 +FIELD IA_FOCUS_MODEL.struct_size
    4 +FIELD IA_FOCUS_MODEL.kind
    4 +FIELD IA_FOCUS_MODEL.baseline_milli
    4 +FIELD IA_FOCUS_MODEL.left_slope_micro_per_step
    4 +FIELD IA_FOCUS_MODEL.right_slope_micro_per_step
    4 +FIELD IA_FOCUS_MODEL.hyperbola_radius_milli
    4 +FIELD IA_FOCUS_MODEL.hyperbola_slope_micro_per_step
END-STRUCTURE

BEGIN-STRUCTURE <IA_FOCUS_SAMPLE>
    4 +FIELD IA_FOCUS_SAMPLE.focus_position
    4 +FIELD IA_FOCUS_SAMPLE.metric_milli
    4 +FIELD IA_FOCUS_SAMPLE.weight_milli
    4 +FIELD IA_FOCUS_SAMPLE.flags
END-STRUCTURE

BEGIN-STRUCTURE <IA_FOCUS_FIT>
    4 +FIELD IA_FOCUS_FIT.struct_size
    4 +FIELD IA_FOCUS_FIT.focus_milli_steps
    4 +FIELD IA_FOCUS_FIT.vertical_offset_milli
    4 +FIELD IA_FOCUS_FIT.rms_residual_milli
    4 +FIELD IA_FOCUS_FIT.points_used
    4 +FIELD IA_FOCUS_FIT.points_left
    4 +FIELD IA_FOCUS_FIT.points_right
    4 +FIELD IA_FOCUS_FIT.flags
END-STRUCTURE

BEGIN-STRUCTURE <IA_FOCUS_SUMMARY>
    4 +FIELD IA_FOCUS_SUMMARY.struct_size
    4 +FIELD IA_FOCUS_SUMMARY.fits_included
    4 +FIELD IA_FOCUS_SUMMARY.median_focus_milli_steps
    4 +FIELD IA_FOCUS_SUMMARY.mad_focus_milli_steps
    4 +FIELD IA_FOCUS_SUMMARY.median_rms_residual_milli
END-STRUCTURE

256 constant IA.MAX_STARS
create ia.config <IA_CONFIG> allot
create ia.summary <IA_FRAME_SUMMARY> allot
IA.MAX_STARS IA_WorkspaceBytes constant ia.workspace.bytes
ia.workspace.bytes allocate throw constant ia.workspace
<IA_STAR> IA.MAX_STARS * allocate throw constant ia.stars

IA.STAR.SATURATED IA.STAR.EDGE or
IA.STAR.BLENDED or IA.STAR.TOO_SMALL or IA.STAR.LOW_SNR or
value ia.reject-flags
32 value ia.summary-stars

ia.config IA_DefaultConfig

: ia.star ( index -- star )
    <IA_STAR> * ia.stars +
;

: ia.milli$ ( n -- caddr u )
\ Format a signed value scaled by 1000 with exactly three decimal places.
    dup >R abs 0 <# # # # '.' hold #s R> sign #>
;

: ia.summarize ( reject-flags -- status )
    ia.stars
    ia.summary IA_FRAME_SUMMARY.measured_stars @ ia.summary-stars min
    rot
    ia.workspace
    ia.workspace.bytes
    ia.summary
    IA_SummarizeStars
;

: ia.measure-locked { img locked count search-radius results -- status }
    img FRAME_BITMAP
    img FRAME_WIDTH @
    img FRAME_HEIGHT @
    img FRAME_WIDTH @
    ia.config
    locked
    count
    search-radius
    results
    IA_MeasureLockedStars
;

variable ia.recommended-exposure

: ia.recommend-exposure
    { current-ms peak background target min-ms max-ms max-change-milli
      -- recommended-ms status }
    current-ms peak background target min-ms max-ms max-change-milli
    ia.recommended-exposure
    IA_RecommendExposure
    ia.recommended-exposure @ swap
;

create ia.roi <IA_RECT> allot

: ia.compute-roi
    { stars count image-width image-height margin -- roi status }
    stars count image-width image-height margin ia.roi
    IA_ComputeStarROI
    ia.roi swap
;

: ia.fit-focus
    { model samples count reject-flags resolution fit -- status }
    <IA_FOCUS_FIT> fit IA_FOCUS_FIT.struct_size !
    model samples count reject-flags resolution fit IA_FitFocusModel
;

: ia.combine-focus
    { fits count reject-flags summary -- status }
    <IA_FOCUS_SUMMARY> summary IA_FOCUS_SUMMARY.struct_size !
    fits count reject-flags ia.workspace ia.workspace.bytes summary
    IA_CombineFocusFits
;

: ia.add-FITS { img | map -- }
    img FRAME_METADATA @ -> map
    s"  " map =>" #STARS"
    ia.summary IA_FRAME_SUMMARY.detected_stars @ (.) map =>" NDETECT"
    ia.summary IA_FRAME_SUMMARY.included_stars @ (.) map =>" NSTARS"
    ia.summary IA_FRAME_SUMMARY.included_stars @ if
        ia.summary IA_FRAME_SUMMARY.median_hfd_milli_pixels @
            ia.milli$ map =>" HFD"
        ia.summary IA_FRAME_SUMMARY.mad_hfd_milli_pixels @
            ia.milli$ map =>" HFDMAD"
        ia.summary IA_FRAME_SUMMARY.median_fwhm_milli_pixels @
            ia.milli$ map =>" FWHM"
        ia.summary IA_FRAME_SUMMARY.mad_fwhm_milli_pixels @
            ia.milli$ map =>" FWHMMAD"
        ia.summary IA_FRAME_SUMMARY.median_ellipticity_milli @
            ia.milli$ map =>" ELLIP"
        ia.summary IA_FRAME_SUMMARY.median_snr_milli @
            ia.milli$ map =>" STARSNR"
    then
    ia.summary IA_FRAME_SUMMARY.background_milli_adu @
        ia.milli$ map =>" BKG"
    ia.summary IA_FRAME_SUMMARY.noise_milli_adu @
        ia.milli$ map =>" BKGNOIS"
;

: compute-starStats { img | status -- status }
\ Analyze the current in-memory 16-bit frame; Forth owns rejection policy.
    <IA_FRAME_SUMMARY> ia.summary IA_FRAME_SUMMARY.struct_size !
    img FRAME_BITMAP
    img FRAME_WIDTH @
    img FRAME_HEIGHT @
    img FRAME_WIDTH @
    ia.config
    ia.workspace
    ia.workspace.bytes
    ia.stars
    IA.MAX_STARS
    ia.summary
    IA_AnalyzeFrame -> status
    status IA.OK = if
        ia.reject-flags ia.summarize drop
    then
    img ia.add-FITS
    status
;

BEGIN-STRUCTURE <FRAME_STATISTICS>
    0x40000 +FIELD HISTOGRAM                        \ one 32 bit cell for each 16 bit brightness value
    0x40000 +FIELD HISTOGRAM_ABSOLUTE_DEVIATION     \ histogram of the absolute deviations of the pixels from the median
          4 +FIELD TOTAL_PIXELS
          4 +FIELD MEAN
          4 +FIELD MEDIAN
          4 +FIELD MEDIAN_ABSOLUTE_DEVIATION
          4 +FIELD df.B                             \ display function parameters
          4 +FIELD df.C 
          4 +FIELD df.MADN
          4 +FIELD df.a
          4 +FIELD df.t
          4 +FIELD df.s
          4 +FIELD df.h
          4 +FIELD df.m          
END-STRUCTURE

: allocate-frameStats ( -- imageStats)
    <FRAME_STATISTICS> allocate abort" unable to allocate image statistics"
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

: compute-histogram { image | imageStats -- }
\ prepare a full-resolution histogram for an image 
    image FRAME_STATISTICS @ -> imageStats
    image FRAME_BITMAP
	imageStats HISTOGRAM dup 0x40000 erase
	image FRAME_SIZE_BYTES @ 2/ dup imageStats TOTAL_PIXELS !       \ each pixel is 2 bytes
    ( bitmap histogram n ) <histogram> 
;

: compute-ASBDhistogram { image | imageStats -- }
\ prepare a full-resolution histogram of the absolute deviation values of image 
    image FRAME_STATISTICS @ -> imageStats
    imageStats MEDIAN @
    image FRAME_BITMAP
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

: histogram.saturated ( imageStats -- count )
\ count the pixels in the saturated 16-bit histogram bin
	HISTOGRAM 0x3fffc + @
;

: add-ImageAnalysisFITS { image | imageStats map -- }
\ add key value pairs for FITS observation parameters
    image FRAME_STATISTICS @ -> imageStats
    image FRAME_METADATA @ -> map
	s"  "                                       map =>" #STATS"         \ a header to indicate the source of these FITS values	
    imageStats MEAN @ (.)                       map =>" MEAN"
    imageStats MEDIAN @ (.)                     map =>" MEDIAN"
    imageStats MEDIAN_ABSOLUTE_DEVIATION @ (.)  map =>" MEDIANAD"
    imageStats histogram.saturated (.)           map =>" SATPIX"
;

: compute-imageStats { image | imageStats }
    image FRAME_STATISTICS @ ?dup if free throw then
    0 image FRAME_STATISTICS !
    allocate-frameStats dup -> imageStats image FRAME_STATISTICS !
    image compute-histogram 
    imageStats compute-mean 
    imageStats compute-median
    image compute-ASBDhistogram              \ must compute the median first
    imagestats compute-median_absolute_deviation
    image add-ImageAnalysisFITS
    image compute-starStats drop
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
			( mean) dest i + w!
	2 +LOOP
;

\ utility functions

: .imageStats ( image --)
    FRAME_STATISTICS @
    cr ." Mean      " dup mean ?
    cr ." Median    " dup median ? 
    cr ." MedianAD  " dup median_absolute_deviation ?
    cr ." DF B      " dup df.B ?                             
    cr ." DF C      " dup df.C ? 
    cr ." DF MADN   " dup df.MADN ?
    cr ." DF a      " dup df.a ?
    cr ." DF t      " dup df.t ?
    cr ." DF s      " dup df.s ?
    cr ." DF h      " dup df.h ?
    cr ." DF m      " dup df.m ?
    drop
;
