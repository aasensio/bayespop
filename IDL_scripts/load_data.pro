;==============================================================================
;
; LOAD_DATA
; This function loads all the spectra and required info from a configuration file
; into a structure. Note that from now on we will work in ln-space
;
; Jesus Falcon-Barroso, IAC, February 2013
;==============================================================================
FUNCTION LOAD_DATA, datafile, idx

cvel = 299792.458d

;# Reading the required spectrum (taking the input format into account)
hdr = headfits('../data/'+datafile.filename[idx],/SILENT)
e3d = strcompress(string(sxpar(hdr,'EURO3D')),/re)
IF e3d EQ 'T' THEN dataformat = 'E3D' ELSE dataformat = 'FITS'

CASE dataformat OF
  'FITS': BEGIN
             spectrum = mrdfits('../data/'+datafile.filename[idx],0,hdr,/SILENT,/DSCALE)
             s_spec = size(spectrum,/DIMENSIONS)
             ndim   = n_elements(s_spec)
             npix   = s_spec[ndim-1]
             IF ndim EQ 1 THEN BEGIN
                nspec = 1
                npix  = s_spec[0]
                dlam  = sxpar(hdr,'CDELT1',/SILENT)
                lam   = sxpar(hdr,'CRVAL1',/SILENT) + dindgen(npix)*dlam
             ENDIF
             IF ndim EQ 2 THEN BEGIN
                nspec = s_spec[1]
                npix  = s_spec[0]
                dlam  = sxpar(hdr,'CDELT1',/SILENT)
                lam   = sxpar(hdr,'CRVAL1',/SILENT) + dindgen(npix)*dlam
             ENDIF
             IF ndim EQ 3 THEN BEGIN
                nspec = s_spec[0]*s_spec[1]
                npix  = s_spec[2]
                dlam  = sxpar(hdr,'CDELT3',/SILENT)
                lam   = sxpar(hdr,'CRVAL3',/SILENT) + dindgen(npix)*dlam
                spectrum = transpose(reform(spectrum,nspec,npix))
             ENDIF
          END
  'E3D':  BEGIN
             tmp = mrdfits('../data/'+datafile.filename[idx],1,hdr,/SILENT,/DSCALE)
             nspec = n_elements(tmp.xpos)
             npix  = n_elements(tmp[0].data_spe)
             dlam  = sxpar(hdr,'CDELTS',/SILENT)
             lam   = sxpar(hdr,'CRVALS',/SILENT) + dindgen(npix)*dlam
             spectrum = tmp.data_spe
          END
  ELSE: message,'ERROR: cannot read '+datafile.dataformat[idx]+' formatted files. Only FITS(1D,2D,3D) or E3D allowed.'
ENDCASE
velscale = cvel*double(dlam)/mean(double(lam))


;# Shifting the wavelength vector to restframe
z   = sqrt((1.0d + datafile.redshift[idx]/cvel)/(1.0d - datafile.redshift[idx]/cvel)) - 1.0d ; Relativistic formula (because the Universe is large!)
lam = lam / (1.0d + z)

;# Log-rebinning the spectra
log_rebin, minmax(lam), spectrum[*,0], spectrum_ln, lam_ln, VELSCALE=velscale
npix_ln = n_elements(lam_ln)
outspec = dblarr(npix_ln,nspec)
FOR i=0, nspec-1 DO BEGIN
    log_rebin, minmax(lam), spectrum[*,i], spectrum_ln, lam_ln, VELSCALE=velscale
    outspec[*,i] = spectrum_ln
ENDFOR

;# Cutting the data outside the LMIN, LMAX (with some margin to account for border effects in convolution)
goodpix = where(lam_ln GE alog(datafile.lmin[0]-50.) AND lam_ln LE alog(datafile.lmax[0]+50.)) ; we leave a buffer of 50AA
IF goodpix[0] EQ -1 THEN BEGIN
   message, ' ERROR: LMIN,LMAX are both outside the spectral range of the data', /INFORMATIONAL, /NONAME
   stop
ENDIF
lam_ln  = lam_ln[goodpix]
outspec = outspec[goodpix,*]
npix_ln = n_elements(goodpix)

;# Creating the structure that will contain the data spectra and information
;  Note that we rebin the spectra to ln-scale
data = {name:datafile.filename[idx], lambda:lam_ln, spec:outspec, npix:npix_ln, nspec:nspec, $
        velscale:velscale, redshift:datafile.redshift[idx], snr:datafile.snr[idx], $
        fwhm:datafile.fwhm[idx], z:datafile.redshift[idx], fix_losvd:datafile.fix_losvd[idx], $
        lmin:alog(datafile.lmin[idx]), lmax:alog(datafile.lmax[idx]), $
        libfile:datafile.libfile[idx], maskfile:datafile.maskfile[idx], priorfile:datafile.priorfile[idx]}

IF data.lmin LT min(lam_ln) THEN BEGIN
   print, 'WARNING: LMIN LT min(lambda). Re-setting LMIN to min(lambda)'
   data.lmin = min(lam_ln)
ENDIF

IF data.lmax GT max(lam_ln) THEN BEGIN
   print, 'WARNING: LMAX GT max(lambda). Re-setting LMAX to max(lambda)'
   data.lmax = max(lam_ln)
ENDIF


;# Printing some info about the spectrum
print, ''
print, ' - Data format:                ', dataformat, FORMAT='(A31,A)'
print, ' - Wave. range (restframe):    ', exp(minmax(data.lambda)), FORMAT='(A31,F8.2,2x,F8.2)'
print, ' - Fitting range (restframe):  ', exp(data.lmin), exp(data.lmax), FORMAT='(A31,F8.2,2x,F8.2)'
print, ' - Velscale (km/s):            ', data.velscale, FORMAT='(A31,F6.2)'
print, ' - Npix (log-rebinned), Nspec: ', data.npix, data.nspec, FORMAT='(A31,I6,2x,I5)'
print, ' - Redshift (km/s):            ', data.redshift, FORMAT='(A31,F8.2)'
print, ' - FWHM (\AA):                 ', data.fwhm, FORMAT='(A31,F5.2)'
print, ' - SNR (per pixel):            ', data.snr, FORMAT='(A31,F6.2)'
print, ''

return, data
END
;-----------------------------------------------
