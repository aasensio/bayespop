;==============================================================================
;
; LOAD_LIBRARY
; This function loads all the spectra and required info from a configuration file
; into a structure.
;
; libfile: a configuration file with the following format:
;         IMF_type (char) IMF_slope (float) Met(float) Age (float) [Mg/Fe] (float)  FILENAME (char) FWHM in angstroms (float)
;
; NOTE: The code assumes all the spectra have the same wavelength, number of pixels and FWHM
; Jesus Falcon-Barroso, IAC, December 2012
;==============================================================================
FUNCTION LOAD_LIBRARY, datastruct

cvel = 299792.458d

libfile = '../config_files/'+datastruct.libfile

;# Reading the library file
readcol, libfile, dummy, libname, FORMAT='A,A', NUMLINE=1, /SILENT
readcol, libfile, dummy, fwhm,    FORMAT='A,F', NUMLINE=2, /SILENT
readcol, libfile, imf_type, imf_slope, met, age,  mgfe, specfile, FORMAT='A,F,F,F,F,A', COMMENT='#', /SILENT
nspec = n_elements(age)

;# Reading the first spectrum of the library to get the dimensions of the models
tmp  = mrdfits(specfile[0],0,hdr,/SILENT,/DSCALE)
npix = n_elements(tmp)
dlam = sxpar(hdr,'CDELT1')
lam  = sxpar(hdr,'CRVAL1') + dindgen(npix)*dlam
velscale = cvel*dlam/mean(lam)

;# Creating the structure that will contain the library spectra and information
library = {age:age,met:met,imf_type:imf_type,imf_slope:imf_slope,mgfe:mgfe, $
           lambda:lam, spec:dblarr(npix,nspec), $
           fwhm:fwhm[0], npix:npix, velscale:velscale, nspec:nspec, libname:libname[0]}

;# Storing all the library spectra into the structure
FOR i=0, nspec-1 DO BEGIN
    spec = mrdfits(specfile[i],0,hdr,/SILENT,/DSCALE)
    library.spec[*,i] = spec
ENDFOR

;# Printing some info about the spectrum
print, ''
print, ' - Wave. range (restframe): ', minmax(library.lambda), FORMAT='(A28,F8.2,2x,F8.2)'
print, ' - Velscale (km/s):         ', library.velscale, FORMAT='(A28,F6.2)'
print, ' - Npix, Nspec:             ', library.npix, library.nspec, FORMAT='(A28,I6,2x,I5)'
print, ' - FWHM (\AA):              ', library.fwhm, FORMAT='(A28,F5.2)'
print, ''

return, library
END
;-----------------------------------------------
