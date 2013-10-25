;------------------------------------------------------------------------------
@load_data.pro
@load_library.pro
@create_mask.pro
@tune_spectra.pro
@adapt_library.pro
;------------------------------------------------------------------------------
PRO SED_FITTING, sedfit_file

datafile = '../config_files/galaxies_list.conf'
libfile  = '../config_files/Library_MILESssp.conf'

cvel = 299792.458d

;## Loading the datafile
print, '# Reading the datafile: '+datafile
readcol, datafile, filename, FORMAT='A', COMMENT='#', /SILENT
nfiles = n_elements(filename)
print, '- '+strcompress(string(nfiles),/re)+' FITS file found'
print, ''

;## Reading the libfile and loading contents to a structure
print, '# Loading library file: '+libfile
library  = load_library(libfile)
iweights = 0+uindgen(library.nspec)*3
ivel     = 1+uindgen(library.nspec)*3
idisp    = 2+uindgen(library.nspec)*3


;## Looping on every spectrum of the datafile ---------------
FOR i=0, nfiles-1 DO BEGIN

    ;## Loading contents of each spectrum to a structure and define mask
    ; NOTE: the data will be put in restframe using the redshift info and log-rebinned
    galaxy = load_data(datafile,i)

   ;# Creating the mask (that includes the limits in wavelength for the fitting)
   print, '# Creating the mask to use on the data'
   mask = create_goodpixels_mask(galaxy,/LOG)

   ;# Adapting the library models to the data
   print, '# Adapting library to match data'
   models  = adapt_library(library, galaxy, mask)
   nmodels = models.nspec

   ;# Setting the priors
   print, '# Setting the priors'
   left  = dblarr(nmodels*3)
   right = dblarr(nmodels*3)
   mu    = dblarr(nmodels*3)
   sigma = dblarr(nmodels*3)
   type  = strarr(nmodels*3)

   ; Priors for model weights (a UNIFORM prior)
   left[iweights]  = 0.0d
   right[iweights] = 1.0d
   type[iweights]  = 'UNIFORM'

   ; Priors for model velocities
   left[ivel]  = -50.0d
   right[ivel] =  50.0d
   mu[ivel]    =  0.0d
   sigma[ivel] =  20.0d
   type[ivel]  = 'NORMAL'

   ; Priors for model velocity dispersions
   left[idisp]  = cvel*(galaxy.fwhm/mean(exp(galaxy.lambda[mask])))/2.355d  ; Instrumental resolution in km/s (sigma)
   right[idisp] = 250.0d          ; Maximum allowed velocity dispersion for a population
   type[idisp]  = 'UNIFORM'


   ;### Starting minimization
   print, '# Starting minimization'
   data  = {galaxy:galaxy, models:models, mask:mask}
   prior = {left:left, right:right, mu:mu, sigma:sigma, type:type}

   ; Creating the MCMC object
   a = obj_new('mcmc',nmodels*3,prior, data=data, /metropolis_within_gibbs)


ENDFOR

END