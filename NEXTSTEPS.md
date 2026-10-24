## Next steps for ImageAnalysis

The native package provides stateless star records and an independent
piecewise-V fitter consuming those records. `StarAnalysis.f` binds both;
AstroImagingInForth supplies a short pipeline-driven Forth autofocus sequence.
The items below
are future work, not functionality supplied by the current DLL; the old
combined analysis/autofocus API has been removed.

1. Calibrate the saved L-filter piecewise-V and hyperbolic models from clear
   multi-star brackets and retain both until repeated nights distinguish them.

2. Field-test the implemented Forth bracket, positional tracking, exposure
   retries and uncertainty threshold. Automatic ROI selection is not part of
   the current routine.

3. Configure predetermined R, G, and B focuser offsets from the fitted L
   position. The native library remains filter-agnostic.

4. Accumulate quality statistics across ordinary science frames and establish
   Forth acceptance/reporting policy from HFD, Gaussian-equivalent FWHM,
   ellipticity, SNR, background, saturation, and model residuals.