FUNCTION ADAPT_LIBRARY, library, data, mask

;# Checking that the wavelength limits of the library are within the desired mask limits
IF exp(min(data.lambda)) LT min(library.lambda) OR exp(max(data.lambda)) GT max(library.lambda) THEN BEGIN
   message, " ERROR: library's wavelength range shorter than data", /INFORMATIONAL, /NONAME
   stop
ENDIF


;# Adapting the library
adapted_library = {spec:dblarr(data.npix,library.nspec), lambda:dblarr(data.npix,library.nspec), $
                   nspec:library.nspec, npix:data.npix, velscale:data.velscale}

FOR i=0, library.nspec-1 DO BEGIN

   ;1.- Convolving to match spectral resolution
   IF data.fwhm LT library.fwhm THEN BEGIN
      message, ' ERROR: library has worse spectral resolution than the data.', /INFORMATIONAL, /NONAME
      stop
   ENDIF

   sigma     = sqrt(data.fwhm^2 - library.fwhm^2)/2.355  ; sigma to convolve with in Angstromgs
   sigma_pix = sigma / (library.lambda[1]-library.lambda[0])

   IF sigma_pix GE 0.1 THEN BEGIN
      lsf      = psf_gaussian(NPIXEL=2*ceil(4*sigma_pix)+1, ST_DEV=sigma_pix, /NORM, NDIM=1)
      tmp_spec = convol(library.spec[*,i],lsf,/EDGE_TRUNCATE) ; Degrading the input spectrum to match the data's spectral resolution
   ENDIF ELSE BEGIN
      tmp_spec = library.spec[*,i]
   ENDELSE

   ;2.- Log-rebinning the template
   log_rebin, minmax(library.lambda), tmp_spec, tmp_spec_ln, tmp_lam_ln, VELSCALE=data.velscale

   ;3.- Interpolating the library to match the data wavelength vector
   ; NOTE: I should really rebin the data, but for the moment let's just do it this way ;-)
   newlam = interpol(findgen(n_elements(tmp_spec_ln)), tmp_lam_ln, data.lambda)
   newtmp = interpolate(tmp_spec_ln,newlam,CUBIC=-0.5, MISSING=0.0)

   ;4.- Normalizing the templates to match the mean flux of the galaxy data in the goodpixels regions
   w = where(mask EQ 1)
   newtmp = newtmp * mean(data.spec[w]) / mean(newtmp[w])

   adapted_library.spec[*,i]   = newtmp
   adapted_library.lambda[*,i] = data.lambda

ENDFOR
print,''

return, adapted_library
END