PRO CREATE_FAKE_DATA, modellist, fwhm_in, weights, vellist, sigmalist, fwhm_out, snr_out,  OUTFILE=outfile, AGES=ages, METS=mets


; Reading in the input SSP models to use as galaxy spectra and convolving to SAURON resolution -------------------------
spec   = mrdfits(modellist[0],0,hdr_spec, /SILENT)
npix   = n_elements(spec)
nspec  = n_elements(modellist)
models = dblarr(npix,nspec)
FOR i  = 0, nspec-1 DO models[*,i] = mrdfits(modellist[i],0,hdr_spec, /SILENT)
dlam   = sxpar(hdr_spec,'CDELT1')
lam    = sxpar(hdr_spec,'CRVAL1') + findgen(npix)*dlam


;# Convolving the models if needed
IF abs(fwhm_in-fwhm_out) GT 0.1 THEN BEGIN

  ssp_spec = dblarr(npix,nspec)
  FOR i=0,nspec-1 DO BEGIN
      fwhm = sqrt(fwhm_out^2 - fwhm_in^2)/dlam ; in pixels
      sigma = fwhm / 2.355d
      lsf  = psf_Gaussian(NPIXEL=2*ceil(4*sigma)+1, ST_DEV=sigma, /NORM, NDIM=1)
      ssp_spec[*,i] = convol(models[*,i],lsf) ; Degrade template to data resolution
  ENDFOR

ENDIF ELSE BEGIN

  ssp_spec = models

ENDELSE

; Generating all the cases to study (output in log-lambda) ------------------------------------------------
FOR i=0, nspec-1 DO BEGIN

   ; log-rebinning the "galaxy" spectra
   log_rebin, minmax(lam), ssp_spec[*,i], ln_spec, ln_lam, VELSCALE=velscale

   IF i EQ 0 THEN newspec = dblarr(n_elements(ln_lam),nspec)

   ; Convolving with the desired LOSVD, and normalizing
   pars         = [vellist[i],sigmalist[i],0.0,0.0]
   newspec[*,i] = broaden_template(ln_spec,velscale,pars)
   newspec[*,i] = newspec[*,i] / mean_biweight(newspec[*,i])

ENDFOR

; Combine spectra and rebin back to linear space
outspec_ln = (newspec # weights) / total(weights)
log_rebin_invert, minmax(ln_lam), outspec_ln, outspec, linlam
outspec = outspec/mean_biweight(outspec)

; Adding random noise for a given S/N
noise   = randomn(seed,n_elements(outspec))/snr_out
outspec = outspec + noise

; Saving the output in a E3D cube
lmin = min(linlam)
dlam = linlam[1] - linlam[0]

mkhdr,    header, outspec
sxaddpar, header, 'CRVAL1', lmin
sxaddpar, header, 'CDELT1', dlam
sxaddpar, header, 'PARAMS', '==========================='
sxaddpar, header, 'AGE', strcompress(strjoin(ages,','),/re)
sxaddpar, header, 'Z', strcompress(strjoin(mets,','),/re)
sxaddpar, header, 'WEIGHTS', strcompress(strjoin(string(weights),','),/re)
sxaddpar, header, 'VEL.', strcompress(strjoin(string(vellist),','),/re)
sxaddpar, header, 'SIGMA', strcompress(strjoin(string(sigmalist),','),/re)
sxaddpar, header, 'SNR', snr_out
sxaddpar, header, 'FWHM', fwhm_out


mwrfits, outspec, outfile, header, /CREATE

END
;------------------------------------------------------------------------------
PRO RUN_CREATE_FAKE_DATA

tempdir = '../templates/'
outdir  = '../data/'
weights   = [0.5,0.5]
vellist   = [0.0d,0.0d]
sigmalist = [100.0d,100.d]


modellist = tempdir+['PEGASE_HR_Zp0.0T10.0000.fits','PEGASE_HR_Zm0.4T01.0000.fits']
fwhm_in   = 0.55d
fwhm_out  = 1.00d
snr_out   = 1000.0d
outfile   = outdir+'testdata_HR.fits'

CREATE_FAKE_DATA, modellist, fwhm_in, weights, vellist, sigmalist, fwhm_out, snr_out, OUTFILE=outfile


modellist = tempdir+['Mku1.30Zp0.00T10.0000.fits','Mku1.30Zm0.40T01.0000.fits']
fwhm_in   = 2.51d
fwhm_out  = 2.51d
snr_out   = 1000.0d
outfile   = outdir+'testdata_LR.fits'

CREATE_FAKE_DATA, modellist, fwhm_in, weights, vellist, sigmalist, fwhm_out, snr_out, OUTFILE=outfile

END
;------------------------------------------------------------------------------
PRO RUN_CREATE_FAKE_DATA_v2

;# Common parameters
tempdir  = '../templates/'
outdir   = '../data/'
confdir  = '../config_files/'

;# Individual cases to test
w1   = [   0.5,    0.85,    0.15,    0.5,    0.85,    0.15,    0.5,    0.5]
w2   = [   0.5,    0.15,    0.85,    0.5,    0.15,    0.85,    0.5,    0.5]
age1 = [  '10',    '10',    '10',   '10',    '10',    '10',   '05',   '10']
age2 = [  '01',    '01',    '01',   '01',    '01',    '01',   '05',   '01']
met1 = ['p0.0',  'p0.0',  'p0.0', 'm0.4',  'm0.4',  'm0.4', 'm0.4', 'm0.4']
met2 = ['m0.4',  'm0.4',  'm0.4', 'p0.0',  'p0.0',  'p0.0', 'p0.0', 'm0.4']
vel1 = [   0.0,     0.0,     0.0,    0.0,     0.0,     0.0,    0.0,    0.0]
vel2 = [-100.0,  -100.0,  -100.0, -100.0,  -100.0,  -100.0, -100.0, -100.0]
sig1 = [ 100.0,   100.0,   100.0,  100.0,   100.0,   100.0,  100.0,  100.0]
sig2 = [  50.0,    50.0,    50.0,   50.0,    50.0,    50.0,   50.0,   50.0]

; ----------------------------------------
fwhm_in   = 0.55d
fwhm_out  = [2.51d,0.55d]
snr_out   = [1000.0d,100.0d,50.0d,25.0d]
root      = 'testdata2'
libfile   = 'PEGASE_HR.conf'
mask      = 'none'
priors    = 'basic_priors.conf'
waverange = [4200.,5500.]
redshift  = 0.0
fix_losvd = 0
conf_file = confdir+'testcases2.conf'
; ----------------------------------------

OPENW, unit, conf_file, /GET_LUN
printf, unit,'#============================================================================================================================================================'
printf, unit,'# FITS file                       FWHM    Redshift (km/s)  S/N   LIBFILE            MASKFILE               PRIORS              FIX_LOSVD   FIT_RANGE MIN MAX'
printf, unit,'#============================================================================================================================================================'


q = 1
FOR k = 0, n_elements(fwhm_out)-1 DO BEGIN
FOR j = 0, n_elements(snr_out)-1 DO BEGIN
FOR i = 0,n_elements(w1)-1 DO BEGIN

    outfile   = outdir+root+'_case'+strcompress(string(q,FORMAT='(I)'),/re)+'.fits'
    print, '- '+outfile

    modellist = tempdir+['PEGASE_HR_Z'+met1[i]+'T'+age1[i]+'.0000.fits','PEGASE_HR_Z'+met2[i]+'T'+age2[i]+'.0000.fits']
    ages      = [age1[i],age2[i]]
    mets      = [met1[i],met2[i]]
    weights   = [w1[i],w2[i]]
    vellist   = [vel1[i],vel2[i]]
    sigmalist = [sig1[i],sig2[i]]

    CREATE_FAKE_DATA, modellist, fwhm_in, weights, vellist, sigmalist, fwhm_out[k], snr_out[j], OUTFILE=outfile, AGES=ages, METS=mets

    printf, unit, outfile, fwhm_out[k], redshift, snr_out[j], libfile, mask, priors, fix_losvd, waverange[0],waverange[1], $
            FORMAT='(A,3x,F5.2,3x,F5.2,3x,F6.1,3x,A,3x,A,3x,A,3x,I,3x,F7.2,3x,F7.2)'


    q = q + 1

ENDFOR
ENDFOR
ENDFOR
CLOSE, unit
print, '# Configuration file: '+conf_file+' created.'

END