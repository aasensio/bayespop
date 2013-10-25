;-------------------------------------------------------------------------------
PRO COMPUTE_MCMC_STATS, savfile, SPEC=spec, STATS=stats, LOSVD=losvd, PLOT=plot

print, '# Computing the statistics for MCMC run: '+savfile
;# Reading the results from the MCMC chain
restore, savfile

;# Setting up the structures that will contain the results
npix  = n_elements(galaxy.spec)
nspec = library.nspec
uage  = library.age[uniq(library.age,sort(library.age))]
nuage = n_elements(uage)
xvel  = range(-1000.,1000.)

spec  = {lambda:galaxy.lambda,data:fltarr(npix), bfit:fltarr(npix), bfit_lo:fltarr(npix), bfit_hi:fltarr(npix), resid:fltarr(npix), mask:mask}
stats = {age:library.age, met:library.met, weights:fltarr(nspec), velocity:fltarr(nspec), sigma:fltarr(nspec), $
         weights_std:fltarr(nspec), velocity_std:fltarr(nspec), sigma_std:fltarr(nspec), $
         age_max:0.0d, met_max:0.0d, weights_max:0.0d, velocity_max:0.0d, sigma_max:0.0d, $
         tot_weights:fltarr(nuage), err_weights:fltarr(nuage), uage:uage, nuage:nuage}
losvd = {losvd:fltarr(n_elements(xvel),nuage), xvel:xvel, age:uage}

;# Computing the statistics for each parameter (biweight mean, robust stddev and maximum value)
ave_pars = fltarr(2,n_elements(outpars[0,*]))
FOR i=0, n_elements(outpars[0,*])-1 DO BEGIN
      tmp = mean_biweight(outpars[*,i], STDDEV=stddev)
      ave_pars[*,i] = [tmp,stddev]
ENDFOR


tmp = reform(ave_pars[0,*],3,library.nspec)
stats.weights  = reform(tmp[0,*])
stats.velocity = reform(tmp[1,*])
stats.sigma    = reform(tmp[2,*])

tmp = reform(ave_pars[1,*],3,library.nspec)
stats.weights_std  = reform(tmp[0,*])
stats.velocity_std = reform(tmp[1,*])
stats.sigma_std    = reform(tmp[2,*])

tmp                = max(stats.weights,mx)
stats.age_max      = stats.age[mx]
stats.met_max      = stats.met[mx]
stats.velocity_max = stats.velocity[mx]
stats.sigma_max    = stats.sigma[mx]

;# Computing the total of weights and LOSVD per age bin
tmp_losvd = fltarr(n_elements(xvel),stats.nuage)
FOR i=0, stats.nuage-1 DO BEGIN

   w   = where(library.age EQ stats.uage[i],nw)
   tmp = fltarr(n_elements(xvel),nw)

   ; Weights and errors
   stats.tot_weights[i] = total(stats.weights[w])
   stats.err_weights[i] = total(stats.weights_std[w])

   ; LOSVDs (using all the results for each age bin)
   FOR j=0,nw-1 DO tmp[*,j] = gaussian(xvel,[stats.weights[w[j]],stats.velocity[w[j]],stats.sigma[w[j]]])
   tmp_losvd[*,i] = total(tmp,2)

ENDFOR
stats.err_weights = stats.err_weights / total(stats.tot_weights)
stats.tot_weights = stats.tot_weights / total(stats.tot_weights)
losvd.losvd = tmp_losvd
losvd.xvel  = xvel
losvd.age   = uage

;# Preparing the different data and bestfit spectra
norm_galaxy  = galaxy.spec/mean(galaxy.spec)
tmp          = tune_spectra(models, PARS=reform(ave_pars[0,*]))
tmp_lo       = tune_spectra(models, PARS=reform(ave_pars[0,*]-ave_pars[1,*]))
tmp_hi       = tune_spectra(models, PARS=reform(ave_pars[0,*]+ave_pars[1,*]))
spec.data    = galaxy.spec/mean(galaxy.spec[mask])
spec.bfit    = tmp/mean(tmp[mask])
spec.bfit_lo = tmp_lo/mean(tmp_lo[mask])
spec.bfit_hi = tmp_hi/mean(tmp_hi[mask])
spec.resid   = galaxy.spec-tmp


;# IF requested preparing the output file for histograms ---------------------------------------------------
IF KEYWORD_SET(PLOT) THEN BEGIN
   PS_Start, FILENAME='test_histograms.ps', CHARSIZE=1.5, /METRIC, /BOLD, /truetype, $
          DEFAULT_THICKNESS=3, XSIZE=21, YSIZE=29, /NOMATCH, XOFFSET=-0.4, YOFFSET=0

   !P.MULTI=[0,3,6]
   k = 0
   FOR i=0, n_elements(outpars[0,*])-1 DO BEGIN

       IF i MOD 3 EQ 0 THEN BEGIN
          title ='Age: '+strcompress(string(library.age[k], FORMAT='(F7.4)'),/re)+', [Z/Z$\sun$]: '+strcompress(string(library.met[k], FORMAT='(F5.2)'),/re)
          k = k + 1
       ENDIF ELSE BEGIN
          title=''
       ENDELSE

       IF i MOD 3 EQ 0 THEN xtitle = 'W!Di!N'
       IF i MOD 3 EQ 1 THEN xtitle = 'Velocity'
       IF i MOD 3 EQ 2 THEN xtitle = '$\sigma$'

       plothist, outpars[*,i], /AUTOBIN, /xs, /fill, /FLINE, FORIENT=45, TITLE=title, xtitle=xtitle
       vertical, val=ave_pars[0,i], COLOR=cgcolor('green'), thick=8
       vertical, val=ave_pars[0,i]-ave_pars[1,i], COLOR=cgcolor('red'), thick=4
       vertical, val=ave_pars[0,i]+ave_pars[1,i], COLOR=cgcolor('red'), thick=4

   ENDFOR

   PS_END

ENDIF

return
END
;-------------------------------------------------------------------------------
