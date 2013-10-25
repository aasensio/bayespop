@tune_spectra.pro
@compute_mcmc_stats.pro
;------------------------------------------------------------------------------
PRO PLOT_SED_FITTING_RESULTS_v4, rootfile

savfile  = '../preproc_data/'+rootfile+'.sav'
mcmcfile = '../results/'+rootfile+'.results'

print, '# Plotting the MCMC results for: '+savfile

;# Computing the MCMC statistics for the plotting
compute_mcmc_stats, mcmcfile, savfile, SPEC=spec, STATS=stats, LOSVD=losvd, /PLOT
stats.uage = stats.uage*1D9

;# Gridding the data for visualization purposes
cleanplot

method = "Kriging"
weights      = griddata(alog10(stats.age),stats.met,stats.weights,      DIMENSION=[stats.nuage,stats.nuage], METHOD=method)
weights_std  = griddata(alog10(stats.age),stats.met,stats.weights_std,  DIMENSION=[stats.nuage,stats.nuage], METHOD=method)
velocity     = griddata(alog10(stats.age),stats.met,stats.velocity,     DIMENSION=[stats.nuage,stats.nuage], METHOD=method)
velocity_std = griddata(alog10(stats.age),stats.met,stats.velocity_std, DIMENSION=[stats.nuage,stats.nuage], METHOD=method)
sigma        = griddata(alog10(stats.age),stats.met,stats.sigma,        DIMENSION=[stats.nuage,stats.nuage], METHOD=method)
sigma_std    = griddata(alog10(stats.age),stats.met,stats.sigma_std,    DIMENSION=[stats.nuage,stats.nuage], METHOD=method)
new_age      = range(min(alog10(stats.uage)),max(alog10(stats.uage)),stats.nuage)
new_met      = range(min(stats.umet),max(stats.umet),stats.nuage)


;# Preparing the output file with the summary ---------------------------------------------------
PS_Start, FILENAME='../figures/'+FILE_BASENAME(savfile,'.sav')+'_summary.ps', CHARSIZE=0.8, /METRIC, /BOLD, /truetype, $
          DEFAULT_THICKNESS=3, XSIZE=21.0, YSIZE=29.7, /NOMATCH,  XOFFSET=0.0, YOFFSET=0.0

!P.MULTI=0

age_range = [0.09,22.]
met_range = [-1.9,0.6]

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
        /xs, YRANGE=[mn*0.4,mx*1.1], ystyle=1, THICK=4, TITLE=FILE_BASENAME(savfile,'.sav')

cgoplot, exp(spec.lambda), spec.data
cgoplot, exp(spec.lambda[spec.mask]), spec.bfit[spec.mask], COL='red', thick=3

cgoplot, minmax(exp(spec.lambda)),[0,0]+mn*0.5, THICK=2
cgoplot, exp(spec.lambda[spec.mask]), spec.resid[spec.mask]+mn*0.5, psym=1, SYMSIZE=0.5, COL='Dark Green'


; SFH ---------------
cgplot, [0],[0], POSITION=[0.075,0.60,0.5,0.75], /NORMAL, /NOERASE, /NODATA, $
         XTITLE='Age (Gyr)', YTITLE='$\Sigma$ W!DL!N', $
         XRANGE=age_range, /xs, xlog=1, YRANGE=[0.0,1.1*max(stats.tot_weights+stats.err_weights)], /ys


cgerrplot, stats.uage*1D-9, stats.tot_weights-stats.err_weights, stats.tot_weights+stats.err_weights, thick=2;, COLOR='ORG6'
cgoplot, stats.uage*1D-9, stats.tot_weights, PSYM=10, COLOR='ORG6'


; Metallicity evolution ---------------
cgplot, [0],[0], POSITION=[0.58,0.60,0.98,0.75], /NORMAL, /NOERASE, /NODATA, $
         XTITLE='Age (Gyr)', YTITLE='Z!DL!N', $
         XRANGE=age_range, /xs, xlog=1, YRANGE=met_range, /ys

cgoplot, stats.uage*1D-9, stats.ave_met, PSYM=10, COLOR='ORG6'
cgoplot, stats.age, stats.met, PSYM=sym(1), SYMSIZE=0.2


;----------------------------------------------
; Age vs Metallicity [Weights] ---------------
;----------------------------------------------
ctload,3,/reverse

ref_lage = [8.0,8.5,9.0,9.5,10.0]
loc_age  = interpol(findgen(n_elements(new_age)),new_age,ref_lage)+0.5
nx       = n_elements(ref_lage)-1
xminor   = 5

ref_met  = [-1.6,-1.2,-0.8,-0.4,0.0,0.4]
loc_met  = interpol(findgen(n_elements(new_met)),new_met,ref_met)+0.5
ny       = n_elements(ref_met)-1
yminor   = 4

struc  = {xticks:nx, xtickv:loc_age, xminor:xminor, xtickformat:'(A1)', $
          yticks:ny, ytickv:loc_met, yminor:yminor, ytickname:string(ref_met,FORMAT='(F5.2)')}

cgimage, weights, /axes, AXKEYWORDS=struc, $
         POSITION=[0.075,0.38,0.46,0.55], /NORMAL, /NOERASE, $
         YTITLE='[Z/Z$\sun$]', STRETCH=1

cgcolorbar, /vertical, /right, POSITION=[0.46,0.39,0.47,0.54], RANGE=[0.0,max(weights)], TITLE='W!DL!N', FORMAT='(F4.2)'

;
; -------------
;

struc    = {xticks:nx, xtickv:loc_age, xminor:xminor, xtickformat:'(A1)', $
            yticks:ny, ytickv:loc_met, yminor:yminor, ytickformat:'(A1)'}

cgimage, weights_std, /axes, AXKEYWORDS=struc, $
         POSITION=[0.55,0.38,0.935,0.55], /NORMAL, /NOERASE, $
         STRETCH=1

cgcolorbar, /vertical, /right, POSITION=[0.935,0.39,0.945,0.54], RANGE=[0.0,max(weights_std)], TITLE='$\Delta$W!DL!N', FORMAT='(F5.3)', DIVISIONS=2

;----------------------------------------------
; Age vs Metallicity [Velocity] ---------------
;----------------------------------------------
sauron_colormap
levels = range(-150,150,151)

struc  = {xticks:nx, xtickv:loc_age, xminor:xminor, xtickformat:'(A1)', $
          yticks:ny, ytickv:loc_met, yminor:yminor, ytickname:string(ref_met,FORMAT='(F5.2)')}

cgimage, velocity, /axes, AXKEYWORDS=struc, $
         POSITION=[0.075,0.21,0.46,0.38], /NORMAL, /NOERASE, $
         YTITLE='[Z/Z$\sun$]', STRETCH=1, MINVALUE=min(levels), MAXVALUE=max(levels)

cgcolorbar, /vertical, /right, POSITION=[0.46,0.22,0.47,0.37], RANGE=minmax(levels), TITLE='V (km/s)', FORMAT='(I4)'

;
; -------------
;

loadct,13,/SILENT
dvlevels = range(0.,10.,11)
struc    = {xticks:nx, xtickv:loc_age, xminor:xminor, xtickformat:'(A1)',  $
            yticks:ny, ytickv:loc_met, yminor:yminor, ytickformat:'(A1)'}

cgimage, velocity_std, /axes, AXKEYWORDS=struc, $
         POSITION=[0.55,0.21,0.935,0.38], /NORMAL, /NOERASE, $
         STRETCH=1, MINVALUE=min(dvlevels), MAXVALUE=max(dvlevels)

cgcolorbar, /vertical, /right, POSITION=[0.935,0.22,0.945,0.37], RANGE=minmax(dvlevels), TITLE='$\Delta$V (km/s)', FORMAT='(I4)', DIVISIONS=2



;---------------------------------------------------------
; Age vs Metallicity [Velocity Dispersion] ---------------
;---------------------------------------------------------
loadct,13,/SILENT
limits = sigrange(stats.sigma)
levels = range(min(limits)-max(stats.sigma_std),max(limits)+max(stats.sigma_std),151)
levels = range(25,300,151)

struc  = {xticks:nx, xtickv:loc_age, xminor:xminor, xtickname:string(ref_lage,FORMAT='(F4.1)'), $
          yticks:ny, ytickv:loc_met, yminor:yminor, ytickname:string(ref_met,FORMAT='(F5.2)')}


cgimage, sigma, /axes, AXKEYWORDS=struc, $
         POSITION=[0.075,0.04,0.46,0.21], /NORMAL, /NOERASE, $
         YTITLE='[Z/Z$\sun$]', XTITLE='Log(Age)  [Gyr]', STRETCH=1, MINVALUE=min(levels), MAXVALUE=max(levels)

cgcolorbar, /vertical, /right, POSITION=[0.46,0.05,0.47,0.20], RANGE=minmax(levels), TITLE='$\sigma$ (km/s)'


;
; -------------
;

dslevels = range(0.,10.,11)
struc    = {xticks:nx, xtickv:loc_age, xminor:xminor, xtickname:string(ref_lage,FORMAT='(F4.1)'),  $
          yticks:ny, ytickv:loc_met, yminor:yminor, ytickformat:'(A1)'}



cgimage, sigma_std, /axes, AXKEYWORDS=struc, $
         POSITION=[0.55,0.04,0.935,0.21], /NORMAL, /NOERASE, $
         XTITLE='Log(Age)  [Gyr]', STRETCH=1, MINVALUE=min(dslevels), MAXVALUE=max(dslevels)

cgcolorbar, /vertical, /right, POSITION=[0.935,0.05,0.945,0.20], RANGE=[minmax(dslevels)], TITLE='$\Delta$$\sigma$ (km/s)', DIVISIONS=2


print, ' - Mean velocity (std):', mean(stats.velocity), mean(stats.velocity_std)
print, ' - Mean sigma (std):', mean(stats.sigma),mean(stats.sigma_std)



;------------------------------------------------
; Reconstructed LOSVD per age bin ---------------
;------------------------------------------------

!P.MULTI=0

ctload,3,/reverse

ref_vel  = range(min(losvd.xvel),max(losvd.xvel),11)
loc_vel  = interpol(findgen(n_elements(losvd.xvel)),losvd.xvel,ref_vel)+0.5
nx       = n_elements(ref_vel)-1
xminor   = 4

ref_lage = alog10(stats.uage)
loc_age  = findgen(n_elements(new_age))+0.5
ny       = n_elements(ref_lage)-1
yminor   = 1


struc    = {xticks:nx, xtickv:loc_vel, xminor:xminor, xtickname:strcompress(string(ref_vel,FORMAT='(I)'),/re), $
            yticks:ny, ytickv:loc_age, yminor:yminor, ytickname:string(ref_lage,FORMAT='(F5.2)')}

cgimage, losvd.losvd, /axes, AXKEYWORDS=struc, POSITION=[0.15,0.5,0.9,0.95], /NORMAL, $
         XTITLE='Velocity (km/s)', YTITLE='Log(Age)  [Gyr]', STRETCH=1, MINVALUE=min(losvd.losvd), MAXVALUE=max(losvd.losvd)

s = size(losvd.losvd,/dimensions)

FOR i=0, ny DO BEGIN

    horizontal, val=i+0.5, color=cgcolor('black'), LINESTYLE=1
    cgtext, s[0]+10, i+0.35, strcompress(string(stats.tot_weights[i]*100., FORMAT='(F5.1)'),/re)+'%', ALIGNMENT=0


ENDFOR

PS_END

END
;---------------------------------------------------------------------------------------------------------
PRO RUN_PLOT_SED_FITTING_RESULTS

list = FILE_SEARCH('../results/testdata2*.results',COUNT=nfiles)

FOR i=0, nfiles-1 DO PLOT_SED_FITTING_RESULTS_v4, FILE_BASENAME(list[i],'.results')

END
