;==============================================================================
;
; TUNE_SPECTRA
;
; Creates a spectrum from a set of spectra given in the structure 'library'

; Parameters of 'library' used here are: nspec, velscale, spec
;
; Jesus Falcon-Barroso, IAC, January 2013
;==============================================================================
FUNCTION TUNE_SPECTRA, library, PARS=pars

;## Some basic error checking
IF n_elements(library) EQ 0  THEN message, 'ERROR (tune_spectra): Library file is not correct'
IF NOT KEYWORD_SET(PARS)     THEN message, 'ERROR (tune_spectra): no PARS defined'

npars = size(pars,/dimensions)
IF library.nspec*3 NE npars[0] THEN message, 'ERROR (tune_spectra): number of PARS not equal to number of spectra in library'

;## Reforming the PARS vector into an array
; Column 0: weights
; Column 1: velocity (in km/s)
; Column 2: velocity dispersion (in km/s)
new_pars = transpose(reform(pars,3,library.nspec))

;## Generate the LOSVDs using new_pars (vectorizing as much as possible)
mx        = max(new_pars[*,2]/library.velscale,nmx)
sigma_pix = new_pars[nmx,2]/library.velscale
npix_lsf  = 2*ceil(4*sigma_pix)+1
lsf       = fltarr(npix_lsf,library.nspec)
midpix    = npix_lsf/2 ; it's important to divide up by an integer!!
x_pix     = findgen(npix_lsf)     # replicate(1,library.nspec)
vel_pix   = replicate(1,npix_lsf) # (new_pars[*,1]/library.velscale)
sig_pix   = replicate(1,npix_lsf) # (new_pars[*,2]/library.velscale)
z         = (x_pix - (midpix+vel_pix))/sig_pix
lsf       = exp(-(z*z)*0.5)
lsf       = lsf / (replicate(1,npix_lsf) # total(lsf,1))  ; normalizing so that total of each LOSVD EQ 1

;## Convolving the spectra with the corresponding LOSVD
tmp_outspec = library.spec*0.0
FOR i=0, library.nspec-1 DO tmp_outspec[*,i] = convol(library.spec[*,i],lsf[*,i])

;## Applying the weights to each spectra to prodice the final spectrum
outspec = (tmp_outspec # new_pars[*,0])/total(new_pars[*,0])

return, outspec
END