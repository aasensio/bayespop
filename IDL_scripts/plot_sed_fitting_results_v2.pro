@tune_spectra.pro
@compute_mcmc_stats.pro
;------------------------------------------------------------------------------
PRO PLOT_SED_FITTING_RESULTS_v2, rootfile

savfile  = '../preproc_data/'+rootfile+'.sav'
mcmcfile = '../results/'+rootfile+'.results'

print, '# Plotting the MCMC results for: '+savfile

; Computing the MCMC statistics for the plotting
compute_mcmc_stats, mcmcfile, savfile, SPEC=spec, STATS=stats, LOSVD=losvd, /PLOT

;# Preparing the output file with the summary ---------------------------------------------------
PS_Start, FILENAME='../figures/'+FILE_BASENAME(savfile,'.sav')+'_summary.ps', CHARSIZE=0.8, /METRIC, /BOLD, /truetype, $
          DEFAULT_THICKNESS=3, XSIZE=21.0, YSIZE=29.7, /NOMATCH,  XOFFSET=0.0, YOFFSET=0.0


!P.MULTI=0
; ctload,1,/reverse
; sauron_colormap

; Computing the 1, 2 and 3sigma levels for the contours
tmp = stats.weights[sort(stats.weights)]
cum = total(tmp,/cumulative)
cum = cum / max(cum)
levels = interpol(tmp,cum,1.0d - [0.9973,0.9545,0.6827])
minlev = interpol(tmp,cum,1.0d - 0.999)


; Plotting the spectrum and bestfit ---------------
mn = min(spec.bfit[spec.mask], MAX=mx)
cgplot, exp(spec.lambda), spec.data, POSITION=[0.075,0.80,0.98,0.97], /NORMAL, $
        XTITLE='Restframe Wavelength ($\angstrom$)', YTITLE='Normalised Flux', $
        /xs, YRANGE=[mn*0.4,mx*1.1], ystyle=1, THICK=4, TITLE=FILE_BASENAME(savfile)

cgoplot, exp(spec.lambda), spec.data
cgoplot, exp(spec.lambda[spec.mask]), spec.bfit[spec.mask], COL='red', thick=3

cgoplot, minmax(exp(spec.lambda)),[0,0]+mn*0.5, THICK=2
cgoplot, exp(spec.lambda[spec.mask]), spec.resid[spec.mask]+mn*0.5, psym=1, SYMSIZE=0.5, COL='Dark Green'

al_legend, ['Data','Bestfit'], color=['black','red'], LINESTYLE=0, THICK=2, $
           LINSIZE=0.35, CHARSIZE=0.8;, POSITION=[!X.CRANGE[0]+0.025*(!X.CRANGE[1]-!X.CRANGE[0]), 0.25]


; SFH ---------------
cgplot, [0],[0], POSITION=[0.075,0.65,0.98,0.75], /NORMAL, /NOERASE, /NODATA, $
         XTICKFORMAT='(A1)', YTITLE='$\Sigma$ W!DL!N', $
         XRANGE=[0.05,20.], /xs, /xlog, YRANGE=[-0.01,1.1*max(stats.tot_weights+stats.err_weights)], /ys

w = where(stats.tot_weights GE 1D-3)

cgerrplot, stats.uage[w], stats.tot_weights[w]-stats.err_weights[w], stats.tot_weights[w]+stats.err_weights[w], thick=2;, COLOR='ORG6'
cgoplot, stats.uage[w], stats.tot_weights[w], PSYM=10, COLOR='ORG6'

; Metallcity Evolution ---------------
cgplot, [0],[0], POSITION=[0.075,0.55,0.98,0.65], /NORMAL, /NOERASE, /NODATA, $
         XTITLE='Age (Gyr)', YTITLE='$\Sigma$ W!DL!N*[Z/Z$\sun$]', $
         XRANGE=[0.05,20.], /xs, /xlog, YRANGE=[-2.45,0.3], /ys

cgerrplot, stats.uage[w], stats.tot_met[w]-stats.err_met[w], stats.tot_met[w]+stats.err_met[w], thick=2;, COLOR='ORG6'
cgoplot, stats.uage[w], stats.tot_met[w], PSYM=10, COLOR='ORG6'
cgoplot, stats.age, stats.met, PSYM=sym(1), SYMSIZE=0.2



; Age vs Metallicity [Weights] ---------------
loadct,13,/SILENT
limits = sigrange(stats.weights)
levels = range(min(limits),max(limits),150)

cgcontour, stats.weights, stats.age, stats.met, /irregular, LEVELS=levels,  /FILL, $
         POSITION=[0.075,0.35,0.5,0.5], /NORMAL, /NOERASE, $
         YTITLE='[Z/Z$\sun$]', XTICKFORMAT='(A1)', $
         XRANGE=[0.075,20.], /xs, /xlog, YRANGE=[-2.45,0.3], /ys, MISSINGVALUE=-999
cgoplot, stats.age, stats.met, PSYM=sym(1), SYMSIZE=0.2

cgcontour, stats.weights_std, stats.age, stats.met, /irregular, LEVELS=levels,  /FILL, $
         POSITION=[0.5,0.35,0.925,0.5], /NORMAL, /NOERASE, $
         XTICKFORMAT='(A1)', YTICKFORMAT='(A1)',$
         XRANGE=[0.075,20.], /xs, /xlog, YRANGE=[-2.45,0.3], /ys
cgoplot, stats.age, stats.met, PSYM=sym(1), SYMSIZE=0.2

cgcolorbar, /vertical, /right, POSITION=[0.93,0.36,0.94,0.49], RANGE=minmax(levels), TITLE='Weights', FORMAT='(F4.2)'


; Age vs Metallicity [Velocity] ---------------
sauron_colormap
limits = sigrange(stats.velocity)
levels = range(-1.0*max(abs(limits)),max(abs(limits)),151)
levels = range(-150,150,151)

cgcontour, stats.velocity, stats.age, stats.met, /irregular, LEVELS=levels,  /FILL, $
         POSITION=[0.075,0.20,0.5,0.35], /NORMAL, /NOERASE, $
         YTITLE='[Z/Z$\sun$]', XTICKFORMAT='(A1)', $
         XRANGE=[0.075,20.], /xs, /xlog, YRANGE=[-2.45,0.3], /ys
cgoplot, stats.age, stats.met, PSYM=sym(1), SYMSIZE=0.2

cgcontour, stats.velocity_std, stats.age, stats.met, /irregular, LEVELS=levels,  /FILL, $
         POSITION=[0.5,0.20,0.925,0.35], /NORMAL, /NOERASE, $
         XTICKFORMAT='(A1)', YTICKFORMAT='(A1)',$
         XRANGE=[0.075,20.], /xs, /xlog, YRANGE=[-2.45,0.3], /ys
cgoplot, stats.age, stats.met, PSYM=sym(1), SYMSIZE=0.2

cgcolorbar, /vertical, /right, POSITION=[0.93,0.21,0.94,0.34], RANGE=minmax(levels), TITLE='Velocity (km/s)', FORMAT='(I4)'


; Age vs Metallicity [Velocity Dispersion] ---------------
loadct,13,/SILENT
limits = sigrange(stats.sigma)
levels = range(min(limits)-max(stats.sigma_std),max(limits)+max(stats.sigma_std),151)
levels = range(25,300,151)

cgcontour, stats.sigma, stats.age, stats.met, /irregular, LEVELS=levels,  /FILL, $
         POSITION=[0.075,0.05,0.5,0.20], /NORMAL, /NOERASE, $
         YTITLE='[Z/Z$\sun$]', XTITLE='Age (Gyr)', $
         XRANGE=[0.075,20.], /xs, /xlog, YRANGE=[-2.45,0.3], /ys
cgoplot, stats.age, stats.met, PSYM=sym(1), SYMSIZE=0.2

cgcontour, stats.sigma_std, stats.age, stats.met, /irregular, LEVELS=levels,  /FILL, $
         POSITION=[0.5,0.05,0.925,0.20], /NORMAL, /NOERASE, $
         YTICKFORMAT='(A1)', XTITLE='Age (Gyr)', $
         XRANGE=[0.075,20.], /xs, /xlog, YRANGE=[-2.45,0.3], /ys
cgoplot, stats.age, stats.met, PSYM=sym(1), SYMSIZE=0.2

cgcolorbar, /vertical, /right, POSITION=[0.93,0.06,0.94,0.19], RANGE=minmax(levels), TITLE='Velocity dispersion (km/s)'

PS_END

stop
END
;---------------------------------------------------------------------------------------------------------
PRO RUN_PLOT_SED_FITTING_RESULTS

; PLOT_SED_FITTING_RESULTS_v2, 'testdata_LR_PEGASE_HR_none_basic_priors'
; PLOT_SED_FITTING_RESULTS_v2, 'testdata_LR_PEGASE_HR_full_none_basic_priors'
PLOT_SED_FITTING_RESULTS_v2, 'testdata_LR_MILESssp_none_basic_priors'

; PLOT_SED_FITTING_RESULTS_v2, 'miles_twopop_singlekin_SN40_MILESssp_none_basic_priors'

END