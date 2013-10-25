FUNCTION MEAN_BIWEIGHT_ARRAY, array, dimension, UNC=unc


s = size(array,/DIMENSIONS)

CASE dimension OF 

 1: BEGIN
      newarray = dblarr(s[0])
      newunc   = dblarr(s[0])

      FOR i=0, s[0]-1 DO BEGIN

         tmp = mean_biweight(array[i,*],STDDEV=tmpstd,/NAN)
         newarray[i] = tmp
	 newunc[i]   = tmpstd

      ENDFOR


     END
     
  2: BEGIN
      newarray = dblarr(s[1])
      newunc   = dblarr(s[1])

      FOR i=0, s[1]-1 DO BEGIN

         tmp = mean_biweight(array[*,i],STDDEV=tmpstd,/NAN)
         newarray[i] = tmp
	 newunc[i]   = tmpstd

      ENDFOR

     END  
ENDCASE


unc = newunc

return, newarray
END
