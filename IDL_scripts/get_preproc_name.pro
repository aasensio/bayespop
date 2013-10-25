FUNCTION GET_PREPROC_NAME, datastruct

    out1 = FILE_BASENAME(datastruct.name)
    pos1 = strsplit(out1,'.')
    out1 = strmid(out1,0,max(pos1)-1)

    out2 = FILE_BASENAME(datastruct.libfile)
    pos2 = strsplit(out2,'.')
    out2 = strmid(out2,0,max(pos2)-1)

    out3 = FILE_BASENAME(datastruct.maskfile)
    pos3 = strsplit(out3,'.')
    out3 = strmid(out3,0,max(pos3)-1)
    IF datastruct.maskfile EQ 'none' THEN out3 = 'none'

    out4 = FILE_BASENAME(datastruct.priorfile)
    pos4 = strsplit(out4,'.')
    out4 = strmid(out4,0,max(pos4)-1)

    outfile = '../preproc_data/'+out1+'_'+out2+'_'+out3+'_'+out4+'.input'

    print, ''
    print, ' - '+outfile
    print, ''

return, outfile
END