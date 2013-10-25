@tune_spectra.pro
@compute_mcmc_stats.pro
;------------------------------------------------------------------------------
PRO PLOT_SED_FITTING_RESULTS, savfile

;savfile  = '../preproc_data/testdata_MILESssp_none_basic_priors.sav'
;savfile  = '../preproc_data/testdata_HR_PEGASE_HR_none_basic_priors.sav'
savfile  = '../preproc_data//miles_twopop_singlekin_MILESssp_full_none_basic_priors.sav'
mcmcfile = '../test.extract.txt'

print, '# Plotting the MCMC results for: '+savfile

; Computing the MCMC statistics for the plotting
compute_mcmc_stats, mcmcfile, savfile, SPEC=spec, STATS=stats, LOSVD=losvd, /PLOT

;# Preparing the output file with the summary ---------------------------------------------------
PS_Start, FILENAME='../figures/'+FILE_BASENAME(savfile,'.sav')+'_summary.ps', CHARSIZE=0.8, /METRIC, /BOLD, /truetype, $
          DEFAULT_THICKNESS=3, XSIZE=29, YSIZE=21, /LANDSCAPE, XOFFSET=-0.75, YOFFSET=27.5


!P.MULTI=0
; ctload,1,/reverse
; sauron_colormap

; Computing the 1, 2 and 3sigma levels for the contours
tmp = stats.weights[sort(stats.weights)]
cum = total(tmp,/cumulative)
cum = cum / max(cum)
levels = interpol(tmp,cum,1.0d - [0.9973,0.9545,0.6827])
minlev = interpol(tmp,cum,1.0d - 0.999)


c_colors = [cgcolor('BLU1'),cgcolor('BLU4'),cgcolor('BLU6')]
c_colors = [cgcolor('GRN1'),cgcolor('GRN3'),cgcolor('GRN6')]
c_colors = [cgcolor('ORG1'),cgcolor('ORG3'),cgcolor('ORG6')]
c_annotation = ['3'+cgsymbol('sigma'),'2'+cgsymbol('sigma'),'1'+cgsymbol('sigma')]
c_linestyle  = [1,3,0]

; Plotting the spectrum and bestfit ---------------
mn = min(spec.bfit[spec.mask], MAX=mx)
cgplot, exp(spec.lambda), spec.data, POSITION=[0.0,0.670,0.53,0.975], /NORMAL, $
        XTITLE='Restframe Wavelength ($\angstrom$)', YTITLE='Normalised Flux', $
        /xs, YRANGE=[0.0,mx*1.1], ystyle=1, THICK=4

cgoplot, exp(spec.lambda), spec.data
cgoplot, exp(spec.lambda[spec.mask]), spec.bfit[spec.mask], COL='red', thick=3

cgoplot, minmax(exp(spec.lambda)),[0,0]+mn*0.5, THICK=2
cgoplot, exp(spec.lambda[spec.mask]), spec.resid[spec.mask]+mn*0.5, psym=1, SYMSIZE=0.5, COL='Dark Green'

al_legend, ['Data','Bestfit'], color=['black','red'], LINESTYLE=0, THICK=2, $
           LINSIZE=0.35, CHARSIZE=0.8, POSITION=[!X.CRANGE[0]+0.025*(!X.CRANGE[1]-!X.CRANGE[0]), 0.25]

; SFH ---------------
cgplot, [0],[0], POSITION=[0.0,0.35,0.53,0.6], /NORMAL, /NOERASE, /NODATA, $
         XTICKFORMAT='(A1)', YTITLE='$\Sigma$ W!DL!N', $
         XRANGE=[0.05,20.], /xs, /xlog, YRANGE=[-0.01,1.1*max(stats.tot_weights+stats.err_weights)], /ys

w = where(stats.tot_weights GE 1D-3)

cgerrplot, stats.uage[w], stats.tot_weights[w]-stats.err_weights[w], stats.tot_weights[w]+stats.err_weights[w], thick=2;, COLOR='ORG6'
cgoplot, stats.uage[w], stats.tot_weights[w], PSYM=10, COLOR='ORG6'

; Metallcity Evolution ---------------
cgplot, [0],[0], POSITION=[0.0,0.1,0.53,0.35], /NORMAL, /NOERASE, /NODATA, $
         XTITLE='Age (Gyr)', YTITLE='$\Sigma$ W!DL!N*[Z/Z$\sun$]', $
         XRANGE=[0.05,20.], /xs, /xlog, YRANGE=[-2.45,0.3], /ys

cgerrplot, stats.uage[w], stats.tot_met[w]-stats.err_met[w], stats.tot_met[w]+stats.err_met[w], thick=2;, COLOR='ORG6'
cgoplot, stats.uage[w], stats.tot_met[w], PSYM=10, COLOR='ORG6'
cgoplot, stats.age, stats.met, PSYM=sym(1), SYMSIZE=0.2


; Age vs Metallicity [Weights] ---------------
cgcontour, stats.weights, stats.age, stats.met, /irregular, LEVELS=levels,  /FILL, /OUTLINE, LABEL=1, $
         C_ANNOTATION=c_annotation, C_COLORS=c_colors, C_LINESTYLE=c_linestyle, $
         POSITION=[0.60,0.682,0.88,0.975], /NORMAL, /NOERASE, C_CHARSIZE=1., $
         YTITLE='[Z/Z$\sun$]', XTICKFORMAT='(A1)', $
         XRANGE=[0.075,20.], /xs, /xlog, YRANGE=[-2.45,0.3], /ys
cgoplot, stats.age, stats.met, PSYM=sym(1), SYMSIZE=0.2
cgoplot, stats.age_max, stats.met_max, PSYM=sym(1), color='Charcoal'
cgoplot, stats.age_max, stats.met_max, PSYM=sym(6)

; Age vs Metallicity [Velocity] ---------------
sauron_colormap

w = where(stats.weights LE minlev)
; w = where(stats.weights GE max(levels))
; stats.velocity[w] = !Values.F_NAN
; stats.sigma[w] = !Values.F_NAN


cgcontour, stats.velocity, stats.age, stats.met, /irregular, LEVELS=range(-300,300,151),  /FILL, $
         POSITION=[0.6,0.391,0.88,0.682], /NORMAL, /NOERASE, C_CHARSIZE=1., $
         YTITLE='[Z/Z$\sun$]', XTICKFORMAT='(A1)', $
         XRANGE=[0.075,20.], /xs, /xlog, YRANGE=[-2.45,0.3], /ys
cgoplot, stats.age, stats.met, PSYM=sym(1), SYMSIZE=0.2
cgoplot, stats.age_max, stats.met_max, PSYM=sym(1), color='Charcoal'
cgoplot, stats.age_max, stats.met_max, PSYM=sym(6)
cgcontour, stats.weights, stats.age, stats.met, /irregular, LEVELS=levels,  /OUTLINE, LABEL=1, $
         C_ANNOTATION=c_annotation, C_COLORS='Charcoal', C_LINESTYLE=c_linestyle, C_CHARSIZE=1., /overplot

cgcolorbar, /vertical, /right, POSITION=[0.89,0.397,0.9,0.676], RANGE=[-300,300], TITLE='Velocity (km/s)'

; Age vs Metallicity [Velocity Dispersion] ---------------
loadct,13,/SILENT
cgcontour, stats.sigma, stats.age, stats.met, /irregular, LEVELS=range(50,300,151),  /FILL, $
         POSITION=[0.6,0.1,0.88,0.391], /NORMAL, /NOERASE, C_CHARSIZE=1., $
         YTITLE='[Z/Z$\sun$]', XTITLE='Age (Gyr)', $
         XRANGE=[0.075,20.], /xs, /xlog, YRANGE=[-2.45,0.3], /ys
cgoplot, stats.age, stats.met, PSYM=sym(1), SYMSIZE=0.2
cgoplot, stats.age_max, stats.met_max, PSYM=sym(1), color='Charcoal'
cgoplot, stats.age_max, stats.met_max, PSYM=sym(6)
cgcontour, stats.weights, stats.age, stats.met, /irregular, LEVELS=levels,  /OUTLINE, LABEL=1, $
         C_ANNOTATION=c_annotation, C_COLORS='Charcoal', C_LINESTYLE=c_linestyle, C_CHARSIZE=1., /overplot

cgcolorbar, /vertical, /right, POSITION=[0.89,0.106,0.9,0.385], RANGE=[50,300], TITLE='Velocity dispersion (km/s)'


PS_END
stop

END
