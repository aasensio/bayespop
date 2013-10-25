;==============================================================================
;
; LOAD_CONFFILE
; This function loads all the required info from a configuration file into a structure. 
;
; conffile: a configuration file with the following format:
;   FITS file (char)
;   FWHM (float) 
;   Redshift (float)  
;   S/N (float) 
;   Libfile (char)
;   Maskfile (char)
;   Fix_LOSVD (byte)
;   Lmin (float) 
;   Lmax (float)
;
; Jesus Falcon-Barroso, IAC, March 2013
;==============================================================================
FUNCTION LOAD_CONFFILE, conffile

cvel = 299792.458d

;# Reading the data file
readcol, conffile, filename, fwhm, redshift, snr, libfile, maskfile, priorsfile, fix_losvd, lmin, lmax, FORMAT='A,F,F,F,A,A,A,B,F,F', COMMENT='#', /SILENT
nfiles = n_elements(filename)

;# Creating the structure that will contain the information stored in the datafile
conf_struc = {filename:filename, fwhm:fwhm, redshift:redshift, snr:snr, $
              lmin:lmin, lmax:lmax, maskfile:maskfile, libfile:libfile, priorfile:priorsfile, $
              fix_losvd:fix_losvd,nfiles:nfiles}

return, conf_struc
END
;-----------------------------------------------
