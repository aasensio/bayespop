FUNCTION CREATE_MASK, datastruct

cvel = 299792.458d

;# Reading the mask file
IF datastruct.maskfile EQ 'none' THEN BEGIN
   mask = lonarr(datastruct.npix)*0 + 1L
   goto, fin
ENDIF

readcol, '../config_files/'+datastruct.maskfile, eml_name,eml_lambda,eml_s, FORMAT='(A,F,F)',COMMENT='#',/SILENT
nlines = n_elements(eml_name)

;# Setting up the basic arrays
mask   = lonarr(datastruct.npix)*0 + 1L

;# Masking all the emission lines defined in the maskfile
FOR i = 0, nlines-1 DO BEGIN

    IF eml_name[i] NE 'sky' THEN BEGIN
       line = eml_lambda[i]
    ENDIF ELSE BEGIN
       z    = sqrt((1.0d + datastruct.redshift/cvel)/(1.0d - datastruct.redshift/cvel)) - 1.0d ; Relativistic formula (because the Universe is large!)
       line = eml_lambda[i]/(1.0d + z) ; Skylines have to be de-redshifted
    ENDELSE

    w = where(datastruct.lambda GE alog(line-0.5*eml_s[i]) AND datastruct.lambda LE alog(line+0.5*eml_s[i]))
    IF w[0] NE -1 THEN mask[w] = 0

ENDFOR

fin:
;# Setting the mask to zero outside the LMIN, LMAX limits for the fit
w = where(datastruct.lambda LT datastruct.lmin OR datastruct.lambda GT datastruct.lmax, COMPLEMENT=goodpix)
mask[w] = 0L

print, ''
print, ' - '+strcompress(string(n_elements(goodpix),FORMAT='(I)'),/re)+' out of '+strcompress(string(datastruct.npix,FORMAT='(I)'),/re)+' will be fitted'
print, ''



return, mask
END
