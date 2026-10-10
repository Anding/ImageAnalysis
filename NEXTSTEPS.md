## Next steps for ImageAnalysis

1. Calibrate the saved L-filter piecewise-V and hyperbolic models from clear
   multi-star brackets and retain both until repeated nights distinguish them.

2. Build the Forth autofocus planning loop around the native calculations:
   discover stars, choose and lock 1-7, revise exposure, quantize the ideal
   subframe, acquire a configurable fixed-exposure bracket, fit each star,
   and accept or retry according to Forth policy.

3. Configure predetermined R, G, and B focuser offsets from the fitted L
   position. The native library remains filter-agnostic.

4. Accumulate quality statistics across ordinary science frames and establish
   Forth acceptance/reporting policy from HFD, Gaussian-equivalent FWHM,
   ellipticity, SNR, background, saturation, and model residuals.