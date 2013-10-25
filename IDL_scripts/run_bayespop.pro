PRO RUN_BAYESPOP

list = FILE_SEARCH('../preproc_data/testdata2_case*input', COUNT=nfiles)

IF  FILE_TEST('../tempdir/',/DIRECTORY) EQ 0 THEN FILE_MKDIR, '../tempdir'

nfiles = 1

FOR i=0, nfiles-1 DO BEGIN

    dirname = '../tempdir/tmp'+strcompress(string(i+1),/re)

;     IF  FILE_TEST(dirname,/DIRECTORY) EQ 0 THEN BEGIN

        FILE_MKDIR, dirname

        FILE_COPY, '../mcmc',dirname+'/.', /OVERWRITE
        FILE_COPY, list[i], dirname+'/internalConf.input', /OVERWRITE

        cd, dirname

        print, '# '+list[i]+'  ------------------'

        SPAWN, './mcmc internalConf.input'
        FILE_COPY, 'test.extract.txt','../../results/'+FILE_BASENAME(list[i],'.input')+'.results', /OVERWRITE

        cd, '../../IDL_scripts'

;     ENDIF

ENDFOR


END