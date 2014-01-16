FUNCTION PREPARE_LIBRARY, library, data, mask, priors

;# Checking that the wavelength limits of the library are within the desired mask limits
IF exp(min(data.lambda)) LT min(library.lambda) OR exp(max(data.lambda)) GT max(library.lambda) THEN BEGIN
   message, " ERROR: library's wavelength range shorter than data", /INFORMATIONAL, /NONAME
   stop
ENDIF


;# Defining the kinematic parameters to be used to prepare the library
step = 50.0d
w         = where(priors.parameter EQ 'Velocity')
nvel      = ceil((priors.max_val[w[0]]-priors.min_val[w[0]])/step)
vel_range = range(priors.min_val[w[0]],priors.max_val[w[0]],nvel)

w         = where(priors.parameter EQ 'Dispersion')
nsig      = ceil((priors.max_val[w[0]]-priors.min_val[w[0]])/step)
sig_range = range(priors.min_val[w[0]],priors.max_val[w[0]],nsig)

meshgrid, vel_range, sig_range, vellist, siglist
vellist = reform(vellist, n_elements(vellist))
siglist = reform(siglist, n_elements(siglist))
ncases  = n_elements(vellist)
ncases = 1

;# Adapting the library
adapted_library = {spec:dblarr(data.npix,library.nspec*ncases), lambda:dblarr(data.npix,library.nspec*ncases), $
                   nspec:library.nspec*ncases, npix:data.npix, velscale:data.velscale, $
                   age:dblarr(library.nspec*ncases), met:dblarr(library.nspec*ncases), $
                   imf_slope:dblarr(library.nspec*ncases), mgfe:dblarr(library.nspec*ncases), $
                   velocity:dblarr(library.nspec*ncases), dispersion:dblarr(library.nspec*ncases)}

k = 0
FOR i=0L, library.nspec-1 DO BEGIN

   ;1.- Convolving to match data's spectral resolution
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

   ;3.- Convolving with all kinematic cases, resampling, normalizing and storing the relevant information
   FOR j=0L, ncases-1 DO BEGIN
       adapted_library.age[k]        = library.age[i]
       adapted_library.met[k]        = library.met[i]
       adapted_library.imf_slope[k]  = library.imf_slope[i]
       adapted_library.mgfe[k]       = library.mgfe[i]
       adapted_library.velocity[k]   = vellist[j]
       adapted_library.dispersion[k] = siglist[j]

       pars    = [vellist[j],siglist[j],0.0,0.0]
       conv_tmp_spec_ln = broaden_template(tmp_spec_ln,data.velscale,pars)


       ;4.- Interpolating the library to match the data wavelength vector
       ; NOTE: I should really rebin the data, but for the moment let's just do it this way ;-)
       newlam = interpol(findgen(n_elements(conv_tmp_spec_ln)), tmp_lam_ln, data.lambda)
       newtmp = interpolate(conv_tmp_spec_ln,newlam,CUBIC=-0.5, MISSING=0.0)

       ;5.- Normalizing the templates to match the mean flux of the galaxy data in the goodpixels regions
       w = where(mask EQ 1)
       newtmp = newtmp * mean(data.spec[w]) / mean(newtmp[w])

       adapted_library.spec[*,k]   = newtmp
       adapted_library.lambda[*,k] = data.lambda

       k = k + 1
   ENDFOR

ENDFOR
print,''

return, adapted_library
END