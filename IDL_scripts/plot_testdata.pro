@compute_mcmc_stats.pro
;-----------------------------------------------------------------------------
PRO PLOT_TESTDATA, test, OUTFILE=outfile, LABEL=label

ntests    = n_elements(test)
rootfiles = 'testdata_'+test+'_PEGASE_HR_none_basic_priors'
savfiles  = '../preproc_data/'+rootfiles+'.sav'
mcmcfiles = '../results/'+rootfiles+'.results'


IF ntests NE 8 THEN return

;# Computing the MCMC statistics for each case
FOR i=0, ntests-1 DO BEGIN
    idx = strcompress(string(i,FORMAT='(I)'),/re)
    tmp = execute('compute_mcmc_stats, [mcmcfiles[i]], [savfiles[i]], STATS=stats'+idx+', LOSVD=losvd'+idx)
    tmp = execute('stats'+idx+'.uage = stats'+idx+'.uage*1D9')
ENDFOR


;# Preparing the output file with the summary ---------------------------------------------------
PS_Start, FILENAME='../figures/'+outfile+'.ps', CHARSIZE=1.5, /METRIC, /BOLD, /truetype, PAGETYPE="A4", $
          DEFAULT_THICKNESS=3, YSIZE=15.05, XSIZE=26.79, /NOMATCH,  XOFFSET=2., YOFFSET=27.7, /LANDSCAPE

!P.MULTI=[0,4,2]


params = ['Weights','Velocity','Sigma']


;------------------------------------
;# Age-Met [Weights,Velocity, Sigma]
;------------------------------------
ctload,3,/reverse

FOR j=0, n_elements(params)-1 DO BEGIN
FOR i=0, ntests-1 DO BEGIN

    idx = strcompress(string(i,FORMAT='(I)'),/re)
    tmp = execute('data    = griddata(alog10(stats'+idx+'.age),stats'+idx+'.met,stats'+idx+'.'+params[j]+', DIMENSION=[stats'+idx+'.nuage,stats'+idx+'.nuage], METHOD="Kriging")')
    tmp = execute('new_age = range(min(alog10(stats'+idx+'.uage)),max(alog10(stats'+idx+'.uage)),stats'+idx+'.nuage)')
    tmp = execute('new_met = range(min(stats'+idx+'.umet),max(stats'+idx+'.umet),stats'+idx+'.nuage)')

    ref_lage = [8.0,8.5,9.0,9.5,10.0]
    loc_age  = interpol(findgen(n_elements(new_age)),new_age,ref_lage)+0.5
    nx       = n_elements(ref_lage)-1
    xminor   = 5

    ref_met  = [-1.6,-1.2,-0.8,-0.4,0.0,0.4]
    loc_met  = interpol(findgen(n_elements(new_met)),new_met,ref_met)+0.5
    ny       = n_elements(ref_met)-1
    yminor   = 4


    IF i EQ 0 THEN BEGIN
       struc = {xticks:nx, xtickv:loc_age, xminor:xminor, xtickformat:'(A1)', $
                yticks:ny, ytickv:loc_met, yminor:yminor, ytickname:string(ref_met,FORMAT='(F5.2)')}
       xtitle = ''
       ytitle = '[Z/Z$\sun$]'
    ENDIF
    IF i GE 1 AND i LE 3 THEN BEGIN
       struc = {xticks:nx, xtickv:loc_age, xminor:xminor, xtickformat:'(A1)', $
                yticks:ny, ytickv:loc_met, yminor:yminor, ytickformat:'(A1)'}
       xtitle = ''
       ytitle = ''
    ENDIF
    IF i EQ 4 THEN BEGIN
       struc = {xticks:nx, xtickv:loc_age, xminor:xminor, xtickname:string(ref_lage,FORMAT='(F4.1)'), $
                yticks:ny, ytickv:loc_met, yminor:yminor, ytickname:string(ref_met,FORMAT='(F5.2)')}
       xtitle = 'Log(Age)  [Gyr]'
       ytitle = '[Z/Z$\sun$]'
    ENDIF
    IF i GE 5 THEN BEGIN
       struc = {xticks:nx, xtickv:loc_age, xminor:xminor, xtickname:string(ref_lage,FORMAT='(F4.1)'), $
                yticks:ny, ytickv:loc_met, yminor:yminor, ytickformat:'(A1)'}
       xtitle = 'Log(Age)  [Gyr]'
       ytitle = ''
    ENDIF

    IF j EQ 0 THEN BEGIN
       limits = minmax(data)
       limits = [0.0,0.2]
       ctload, 3, /reverse
    ENDIF
    IF j EQ 1 THEN  BEGIN
       limits = [-150.,150.]
       sauron_colormap
    ENDIF
    IF j EQ 2 THEN  BEGIN
       limits = [50.,300.]
       ctload, 13
    ENDIF

    cgimage, data, /axes, AXKEYWORDS=struc, XTITLE=xtitle, YTITLE=ytitle, STRETCH=1, MAXVALUE=max(limits), MINVALUE=min(limits), /SAVE

ENDFOR

;  cgcolorbar, POSITION=[0.3,0.95,0.7,1.1],  RANGE=minmax(limits), TITLE=params[j], DIVISIONS=2, /FIT
xyouts, 0.5, 1.1, params[j], /NORMAL, ALIGNMENT=0.5, CHARSIZE=2, COLOR=cgcolor('black')
xyouts, 0.5, 1.05, label, /NORMAL, ALIGNMENT=0.5, CHARSIZE=1., COLOR=cgcolor('black')

ENDFOR

;--------------
;# Age-LOSVD
;--------------
!P.MULTI=[0,4,2]

ctload,3,/reverse
FOR i=0, ntests-1 DO BEGIN

    idx = strcompress(string(i,FORMAT='(I)'),/re)
    tmp = execute('data    = losvd'+idx+'.losvd')
    tmp = execute('new_age = range(min(alog10(stats'+idx+'.uage)),max(alog10(stats'+idx+'.uage)),stats'+idx+'.nuage)')
    tmp = execute('ref_vel = [-900,-600,-300,0,300,600,900]')
    tmp = execute('loc_vel = interpol(findgen(n_elements(losvd'+idx+'.xvel)),losvd'+idx+'.xvel,ref_vel)+0.5')
    nx       = n_elements(ref_vel)-1
    xminor   = 3

    tmp = execute('ref_lage = alog10(stats'+idx+'.uage)')
    loc_age  = findgen(n_elements(new_age))+0.5
    ny       = n_elements(ref_lage)-1
    yminor   = 1

    IF i EQ 0 THEN BEGIN
       struc    = {xticks:nx, xtickv:loc_vel, xminor:xminor, xtickformat:'(A1)', $
                   yticks:ny, ytickv:loc_age, yminor:yminor, ytickname:string(ref_lage,FORMAT='(F5.2)')}
       xtitle = ''
       ytitle = 'Log(Age)  [Gyr]'
    ENDIF
    IF i GE 1 AND i LE 3 THEN BEGIN
       struc = {xticks:nx, xtickv:loc_vel, xminor:xminor, xtickformat:'(A1)', $
                yticks:ny, ytickv:loc_age, yminor:yminor, ytickformat:'(A1)'}
       xtitle = ''
       ytitle = ''
    ENDIF
    IF i EQ 4 THEN BEGIN
       struc    = {xticks:nx, xtickv:loc_vel, xminor:xminor, xtickname:strcompress(string(ref_vel,FORMAT='(I)'),/re), $
                   yticks:ny, ytickv:loc_age, yminor:yminor, ytickname:string(ref_lage,FORMAT='(F5.2)')}
       xtitle = 'Velocity (km/s)'
       ytitle = 'Log(Age)  [Gyr]'
    ENDIF
    IF i GE 5 THEN BEGIN
       struc    = {xticks:nx, xtickv:loc_vel, xminor:xminor, xtickname:strcompress(string(ref_vel,FORMAT='(I)'),/re), $
                   yticks:ny, ytickv:loc_age, yminor:yminor, ytickformat:'(A1)'}
       xtitle = 'Velocity (km/s)'
       ytitle = ''
    ENDIF


    cgimage, data, /axes, AXKEYWORDS=struc, XTITLE=xtitle, YTITLE=ytitle, STRETCH=1, MINVALUE=min(data), MAXVALUE=max(data)

    FOR j=0, ny DO horizontal, val=j+0.5, color=cgcolor('black'), LINESTYLE=1


ENDFOR
xyouts, 0.5, 1.1, 'LOSVD', /NORMAL, ALIGNMENT=0.5, CHARSIZE=2, COLOR=cgcolor('black')
xyouts, 0.5, 1.05, label, /NORMAL, ALIGNMENT=0.5, CHARSIZE=1., COLOR=cgcolor('black')

PS_END,/PDF

END
;------------------------------------------------------------------------------
PRO RUN_PLOT_TESTDATA


;# Reading all the available info about each case
print, '# Loading info of testcases'

list = FILE_SEARCH('../data/testdata_case*.fits',COUNT=ncases)

test = strarr(ncases)
age1 = fltarr(ncases)
age2 = fltarr(ncases)
met1 = strarr(ncases)
met2 = strarr(ncases)
vel1 = fltarr(ncases)
vel2 = fltarr(ncases)
sig1 = fltarr(ncases)
sig2 = fltarr(ncases)
wei1 = fltarr(ncases)
wei2 = fltarr(ncases)
snr  = fltarr(ncases)
fwhm = fltarr(ncases)
FOR i=0, ncases-1 DO BEGIN

    hdr = headfits(list[i])

    tmp  = strsplit(FILE_BASENAME(list[i],'.fits'),'_',/EXTRACT)
    test[i] = strcompress(tmp[1],/re)

    tmp_age = strsplit(sxpar(hdr,'AGE'),',',/extract)
    age1[i] = strcompress(tmp_age[0],/re)
    age2[i] = strcompress(tmp_age[1],/re)

    tmp_met = strsplit(sxpar(hdr,'Z'),',',/extract)
    met1[i] = strcompress(tmp_met[0],/re)
    met2[i] = strcompress(tmp_met[1],/re)

    tmp_wei = strsplit(sxpar(hdr,'WEIGHTS'),',',/extract)
    wei1[i] = strcompress(tmp_wei[0],/re)
    wei2[i] = strcompress(tmp_wei[1],/re)

    tmp_vel = strsplit(sxpar(hdr,'VEL.'),',',/extract)
    vel1[i] = strcompress(tmp_vel[0],/re)
    vel2[i] = strcompress(tmp_vel[1],/re)

    tmp_sig = strsplit(sxpar(hdr,'SIGMA'),',',/extract)
    sig1[i] = strcompress(tmp_sig[0],/re)
    sig2[i] = strcompress(tmp_sig[1],/re)

    fwhm[i] = sxpar(hdr,'FWHM')
    snr[i]  = sxpar(hdr,'SNR')

ENDFOR


;# Selecting cases and producing plots
w1 = [   0.5,    0.85,    0.15,    0.5,    0.85,    0.15,    0.5,    0.5]
w2 = [   0.5,    0.15,    0.85,    0.5,    0.15,    0.85,    0.5,    0.5]
a1 = [   10.,     10.,     10.,    10.,     10.,     10.,    5.0,    10.]
a2 = [   1.0,     1.0,     1.0,    1.0,     1.0,     1.0,    5.0,    1.0]
m1 = ['p0.0',  'p0.0',  'p0.0', 'm0.4',  'm0.4',  'm0.4', 'm0.4', 'm0.4']
m2 = ['m0.4',  'm0.4',  'm0.4', 'p0.0',  'p0.0',  'p0.0', 'p0.0', 'm0.4']
v1 = [   0.0,     0.0,     0.0,    0.0,     0.0,     0.0,    0.0,    0.0]
v2 = [   0.0,     0.0,     0.0,    0.0,     0.0,     0.0,    0.0,    0.0]
s1 = [ 100.0,   100.0,   100.0,  100.0,   100.0,   100.0,  100.0,  100.0]
s2 = [ 100.0,   100.0,   100.0,  100.0,   100.0,   100.0,  100.0,  100.0]

FOR i=0, n_elements(w1)-1 DO BEGIN

    w = where(wei1 EQ w1[i]  AND wei2 EQ w2[i]  AND $
              age1 EQ a1[i]  AND age2 EQ a2[i]  AND $
              met1 EQ m1[i]  AND met2 EQ m2[i]  AND $
              vel1 EQ v1[i]  AND vel2 EQ v2[i]  AND $
              sig1 EQ s1[i]  AND sig2 EQ s2[i])

    label = 'Weights = ['+strcompress(string(wei1[w[0]]),/re)+','+strcompress(string(wei2[w[0]]),/re)+'], '+ $
            'Age = ['+strcompress(string(age1[w[0]]),/re)+','+strcompress(string(age2[w[0]]),/re)+'], '+ $
            'Z = ['+met1[w[0]]+','+met2[w[0]]+'], '+ $
            'Vel. = ['+strcompress(string(vel1[w[0]]),/re)+','+strcompress(string(vel2[w[0]]),/re)+'], '+ $
            'Sig. = ['+strcompress(string(sig1[w[0]]),/re)+','+strcompress(string(sig2[w[0]]),/re)+']'
    seq = multisort(fwhm[w],snr[w],ORDER=[1,-1])
    PLOT_TESTDATA, test[w[seq]], OUTFILE='results'+strcompress(string(i+1,FORMAT='(I)'),/re), LABEL=label

;     FILE_DELETE, '../figures/results'+strcompress(string(i+1,FORMAT='(I)'),/re)+'.ps', /ALLOW_NONEXISTENT

ENDFOR


END