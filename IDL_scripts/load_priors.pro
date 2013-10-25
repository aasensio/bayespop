;==============================================================================
;
; LOAD_PRIORS
; This function loads all the priors required for the fitting into a structure.
;
; conffile: a configuration file with the following format:
;   Parameter (char):  parameter to impose a prior to. 6 options: Velocity, Dispersion, Age, Metallicity, IMF, [Mg/Fe]
;   Type (char):  type of prior. Two options: UNIFORM or NORMAL
;   Min (float):  lower bound for prior
;   Max (float):  upper bound for prior
;   Mean (float): in case of NORMAL, Mean of the distribution. Otherwise ignored
;   Sig (float):  in case of NORMAL, dispersion of the distribution. Otherwise ignored
;
; Jesus Falcon-Barroso, IAC, October 2013
;==============================================================================
FUNCTION LOAD_PRIORS, datastruct, libstruct

cvel = 299792.458d

priorsfile = '../config_files/'+datastruct.priorfile

;# Reading the data file
readcol, priorsfile, parameter, type, min_val, max_val, mean_val, sig_val, FORMAT='A,A,F,F,F,F', COMMENT='#', /SILENT
npriors = n_elements(parameter)

;# Creating the structure that will contain the information stored in the datafile
prior_struc = {parameter:strarr(6),type:strarr(6), min_val:fltarr(6), max_val:fltarr(6), mean_val:fltarr(6), sig_val:fltarr(6)}


;# Checking that all the required priors are set.
;NOTE: If not specified, library parameters are used as default assuming an UNIFORM distribution

w = where(parameter EQ 'Velocity')
IF w[0] EQ -1 THEN BEGIN
   prior_struc.parameter[0] = 'Velocity'
   prior_struc.type[0]      = 'UNIFORM'
   prior_struc.min_val[0]   = -300.0
   prior_struc.max_val[0]   =  300.0
   print, '  [Using default values for priors in Velocity]'
ENDIF ELSE BEGIN
   prior_struc.parameter[0] = parameter[w[0]]
   prior_struc.type[0]      = type[w[0]]
   prior_struc.min_val[0]   = min_val[w[0]]
   prior_struc.max_val[0]   = max_val[w[0]]
   prior_struc.mean_val[0]  = mean_val[w[0]]
   prior_struc.sig_val[0]   = sig_val[w[0]]
ENDELSE


w = where(parameter EQ 'Dispersion')
IF w[0] EQ -1 THEN BEGIN
   prior_struc.parameter[1] = 'Dispersion'
   prior_struc.type[1]      = 'UNIFORM'
   prior_struc.min_val[1]   =   50.0
   prior_struc.max_val[1]   =  400.0
   print, '  [Using default values for priors in Dispersion]'
ENDIF ELSE BEGIN
   prior_struc.parameter[1] = parameter[w[0]]
   prior_struc.type[1]      = type[w[0]]
   prior_struc.min_val[1]   = min_val[w[0]]
   prior_struc.max_val[1]   = max_val[w[0]]
   prior_struc.mean_val[1]  = mean_val[w[0]]
   prior_struc.sig_val[1]   = sig_val[w[0]]
ENDELSE

w = where(parameter EQ 'Age')
IF w[0] EQ -1 THEN BEGIN
   prior_struc.parameter[2] = 'Age'
   prior_struc.type[2]      = 'UNIFORM'
   prior_struc.min_val[2]   = min(libstruct.age)
   prior_struc.max_val[2]   = max(libstruct.age)
   print, '  [Using default values for priors in Age]'
ENDIF ELSE BEGIN
   prior_struc.parameter[2] = parameter[w[0]]
   prior_struc.type[2]      = type[w[0]]
   prior_struc.min_val[2]   = min_val[w[0]]
   prior_struc.max_val[2]   = max_val[w[0]]
   prior_struc.mean_val[2]  = mean_val[w[0]]
   prior_struc.sig_val[2]   = sig_val[w[0]]
ENDELSE

w = where(parameter EQ 'Metallicity')
IF w[0] EQ -1 THEN BEGIN
   prior_struc.parameter[3] = 'Metallicity'
   prior_struc.type[3]      = 'UNIFORM'
   prior_struc.min_val[3]   = min(libstruct.met)
   prior_struc.max_val[3]   = max(libstruct.met)
   print, '  [Using default values for priors in Metallicity]'
ENDIF ELSE BEGIN
   prior_struc.parameter[3] = parameter[w[0]]
   prior_struc.type[3]      = type[w[0]]
   prior_struc.min_val[3]   = min_val[w[0]]
   prior_struc.max_val[3]   = max_val[w[0]]
   prior_struc.mean_val[3]  = mean_val[w[0]]
   prior_struc.sig_val[3]   = sig_val[w[0]]
ENDELSE

w = where(parameter EQ 'IMF')
IF w[0] EQ -1 THEN BEGIN
   prior_struc.parameter[4] = 'IMF'
   prior_struc.type[4]      = 'UNIFORM'
   prior_struc.min_val[4]   = min(libstruct.imf_slope)
   prior_struc.max_val[4]   = max(libstruct.imf_slope)
   print, '  [Using default values for priors in IMF]'
ENDIF ELSE BEGIN
   prior_struc.parameter[4] = parameter[w[0]]
   prior_struc.type[4]      = type[w[0]]
   prior_struc.min_val[4]   = min_val[w[0]]
   prior_struc.max_val[4]   = max_val[w[0]]
   prior_struc.mean_val[4]  = mean_val[w[0]]
   prior_struc.sig_val[4]   = sig_val[w[0]]
ENDELSE

w = where(parameter EQ '[Mg/Fe]')
IF w[0] EQ -1 THEN BEGIN
   prior_struc.parameter[5] = '[Mg/Fe]'
   prior_struc.type[5]      = 'UNIFORM'
   prior_struc.min_val[5]   = min(libstruct.mgfe)
   prior_struc.max_val[5]   = max(libstruct.mgfe)
   print, '  [Using default values for priors in Mg/Fe]'
ENDIF ELSE BEGIN
   prior_struc.parameter[5] = parameter[w[0]]
   prior_struc.type[5]      = type[w[0]]
   prior_struc.min_val[5]   = min_val[w[0]]
   prior_struc.max_val[5]   = max_val[w[0]]
   prior_struc.mean_val[5]  = mean_val[w[0]]
   prior_struc.sig_val[5]   = sig_val[w[0]]
ENDELSE

;# Printing some info about the priors
print, ''
print, ' - Parameters: ', prior_struc.parameter
print, ' - Type: ', prior_struc.type
print, ' - Min_val: ', prior_struc.min_val
print, ' - Max_val: ', prior_struc.max_val
print, ' - Mean_val: ', prior_struc.mean_val
print, ' - Sig_val : ', prior_struc.sig_val
print, ''


return, prior_struc
END
;-----------------------------------------------
