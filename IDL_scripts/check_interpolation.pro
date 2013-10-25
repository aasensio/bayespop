PRO CREATE_FAKE_DATA, modellist, vellist, sigmalist, OUTSPEC=outspec


; Reading in the input SSP models to use as galaxy spectra and convolving to the required resolution -------------------------
spec     = mrdfits(modellist[0],0,hdr_spec, /SILENT)
npix     = n_elements(spec)
dlam     = sxpar(hdr_spec,'CDELT1')
lam      = sxpar(hdr_spec,'CRVAL1') + findgen(npix)*dlam
ssp_spec = spec

; Generating all the cases to study (output in log-lambda) ------------------------------------------------

   ; log-rebinning the "galaxy" spectra
   log_rebin, minmax(lam), ssp_spec, ln_spec, ln_lam, VELSCALE=velscale

   ; Convolving with the desired LOSVD, and normalizing
   outspec = dblarr(n_elements(ln_lam),n_elements(vellist))
   FOR i=0, n_elements(vellist)-1 DO BEGIN
       pars    = [vellist[i],sigmalist[i],0.0,0.0]
       outspec[*,i] = broaden_template(ln_spec,velscale,pars)
   ENDFOR

return
END
;------------------------------------------------------------------------------
PRO CHECK_INTERPOLATION

;# Common parameters
tempdir   = '../templates/'
modellist = tempdir+['PEGASE_HR_Zp0.0T10.0000.fits']
; modellist = tempdir+['Mku1.30Zp0.00T10.0000.fits']

;# Creating the database
step = 50.
vel= -300.0 + findgen(9) * step
sig =  20.0 + findgen(10) * step
meshgrid, vel, sig, vellist, sigmalist
vellist   = reform(vellist, n_elements(vellist))
sigmalist = reform(sigmalist, n_elements(sigmalist))

CREATE_FAKE_DATA, modellist, vellist, sigmalist, OUTSPEC=outspec
vel_database  = vellist
sig_database  = sigmalist
spec_database = outspec

;# Create random individual cases and compare with interpolated from database
ncases = 5000

vellist   = scale_vector(randomu(seed,ncases), min(vel), max(vel))
sigmalist = scale_vector(randomu(seed,ncases), min(sig),max(sig))

CREATE_FAKE_DATA, modellist, vellist, sigmalist, OUTSPEC=outspec
good_spec = outspec


;# Computing the interpolated spectra for each testcase
interp_spec = good_spec * 0.0
ave = dblarr(ncases)
dev = dblarr(ncases)

FOR i=0, ncases-1 DO BEGIN

    d = sqrt((vel_database - vellist[i])^2 + (sig_database - sigmalist[i])^2)
    idx = where(abs(vel_database - vellist[i]) LT step AND abs(sig_database - sigmalist[i]) LT step)

    weights = (1.0d/d[idx])^2.
    interp_spec[*,i] = (1.0d/total(weights)) * (spec_database[*,idx] # weights) ; weighting with the inverse of the distance squared

;     plot, vel_database, sig_database, psym=1, /xs, /ys
;     oplot, vel_database[idx], sig_database[idx], psym=sym(1), color=cgcolor('red')
;     oplot, [vellist[i]], [sigmalist[i]], psym=sym(1)    , color=cgcolor('green')

    ave[i] = mean_biweight(good_spec[*,i]/interp_spec[*,i], STDDEV=tmp)
    dev[i] = tmp

; stop
ENDFOR

sauron_colormap
cgerase
!P.MULTI=[0,2,2]
!P.CHARSIZE=2
cgcontour, ave, vellist, sigmalist, /irregular, /cell_fill, levels=range(0.8,1.2,50), $;/xs, /ys, $
           XTITLE='Velocity (km/s)', YTITLE='Velocity Dispersion (km/s)'
oplot, vel_database, sig_database, psym=1

histogauss, ave, kk, color=0 ,CHARSIZE=1.2

cgcontour, dev, vellist, sigmalist, /irregular, /cell_fill, levels=range(0.0,max(dev),50), $;/xs, /ys, $
           XTITLE='Velocity (km/s)', YTITLE='Velocity Dispersion (km/s)'
oplot, vel_database, sig_database, psym=1

histogauss, dev, kk, color=0 ,CHARSIZE=1.2

print, 'Minmax (%)',minmax(dev)*100.

stop
END
