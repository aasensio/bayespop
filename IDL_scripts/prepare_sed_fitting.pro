;------------------------------------------------------------------------------
; Loading all the required modules
;------------------------------------------------------------------------------
@create_mask.pro
@get_preproc_name.pro
@load_conffile.pro
@load_data.pro
@load_library.pro
@load_priors.pro
@prepare_library.pro
;------------------------------------------------------------------------------
PRO PREPARE_SED_FITTING, configuration_file

IF n_elements(configuration_file) EQ 0 THEN configuration_file = '../config_files/galaxies.conf'

;## Loading the datafile
print, '# Reading the configuration file: '+configuration_file
conffile = load_conffile(configuration_file)
print, '- '+strcompress(string(conffile.nfiles),/re)+' case(s) found'
print, ''

;## Looping on every case of the datafile ---------------
FOR i=0, conffile.nfiles-1 DO BEGIN

    ;## Loading contents of each case in the conffile into a structure
    ; NOTE: the data will be put in restframe using the redshift info and log-rebinned
    len = strlen(conffile.filename[i])
    print, strjoin(replicate('=',len+4))
    print, '  '+conffile.filename[i]
    print, strjoin(replicate('=',len+4))
    galaxy = load_data(conffile,i)

    ;## Creating the mask (1 for the good pixels and 0 for the ones to mask for the fitting)
    print, '# Creating the mask to use on the data: '+galaxy.maskfile
    mask = create_mask(galaxy)

    ;## Reading the libfile and loading contents to a structure
    print, '# Loading library file: '+galaxy.libfile
    library  = load_library(galaxy)

    ;## Reading the priors conffile
    print, '# Loading the priors file: '+galaxy.priorfile    
    priors = load_priors(galaxy,library)

    ;## Adapting the library models to the data
    print, '# Adapting library to match data (Be patient. This can take a few minutes)'
    models  = prepare_library(library, galaxy, mask, priors)

    ;## Defining the name of the preproc file for this case
    ;NOTE: We strip whatever is after the last '.'
    print, '# Defining output name for prepoc_data:'
    outfile = get_preproc_name(galaxy)

    ;## Writing preprocessed data into a binary file
    print, '# Saving the necessary data to the preproc_data file'

    IF NOT FILE_TEST('../preproc_data',/DIRECTORY) THEN FILE_MKDIR, '../preproc_data'
    OPENW, unit, outfile, /F77_UNFORMATTED, /GET_LUN
    root = string(replicate(32B,100))
    strput, root, "results/test"
    WRITEU, unit, root

      ; Prior's related parameters
      FOR j=0, 5 DO BEGIN
          CASE priors.type[j] OF
            'UNIFORM': BEGIN
                         WRITEU, unit, 0L
                         WRITEU, unit, 1.d0*priors.min_val[j], 1.d0*priors.max_val[j]  ; lower and upper bound for parameter (scalar [double])
                         WRITEU, unit, 1.d0*priors.mean_val[j]                         ; initial guess for the uniform distributionr (scalar [double])
                       END
            'NORMAL':  BEGIN
                         WRITEU, unit, 1L
                         WRITEU, unit, 1.d0*priors.min_val[j], 1.d0*priors.max_val[j]  ; lower and upper bound for parameter (scalar [double])
                         WRITEU, unit, 1.d0*priors.mean_val[j], 1.d0*priors.sig_val[j] ; mean and sigma of normal distribution (scalar [double])
                       END
            ELSE: message,'ERROR: prior for Weights not UNIFORM or NORMAL'
          ENDCASE


      ENDFOR

      ; Model's related parameters
      WRITEU, unit, 1L*models.npix
      WRITEU, unit, 1L*models.nspec
      WRITEU, unit, 1.d0*models.velscale
      WRITEU, unit, 1.d0*models.spec
      WRITEU, unit, 1.d0*models.age
      WRITEU, unit, 1.d0*models.met
      WRITEU, unit, 1.d0*models.imf_slope
      WRITEU, unit, 1.d0*models.mgfe

      ; Data's related parameters
      WRITEU, unit, 1.d0*galaxy.snr
      WRITEU, unit, 1L*galaxy.nspec
      WRITEU, unit, 1.d0*galaxy.spec
      WRITEU, unit, 1L*mask
      WRITEU, unit, 1L*galaxy.fix_losvd            

   CLOSE, unit
   stop
   FREE_LUN, unit
   print, ' '
   print, '# Preprocessed file '+outfile+' written on disk'
   print, ' '

ENDFOR

END