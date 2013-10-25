; Simple code to show the results of the code

; doplot, 'OBSERVATIONS/test.obs', 'RESULTS/test'
pro doplot, length
	
	openr,2,'../bestFit.dat'
	readf,2,npar
	best = dblarr(npar)
	readf,2,best
	readf,2,nl
	spec = dblarr(2,nl)
	readf,2,spec
	close,2

stop
	cgplot, spec[0,*]
	cgoplot, spec[1,*], col='red'	
	
	t = dblarr(npar,length)
	
	openr,2,'../test.extract.txt'
	readf,2,t
	close,2
	
	cgdisplay,xsize=1200,ysize=800

	!p.multi = [0,6,5]
	for i = 0, npar-1 do begin
		if (stddev(t[i,*]) gt 1.d-6) then begin
			cgplot, t[i,*], /ynozero
		endif

	endfor
	
; 	for i = 0, npar-1 do begin
; 		if (stddev(t[i,*]) gt 1.d-6) then begin
; 			cghistoplot, t[i,*]
; 		endif
; 	endfor
	
	nSpec = npar / 3
	mn = dblarr(npar)
	std = dblarr(npar)
	for i = 0, npar-1 do begin
		mn[i] = mean(t[i,*])
		std[i] = stddev(t[i,*])
	endfor
	
	cgplot, mn[0:nSpec-1], psym=4
	cgerrplot, mn[0:nSpec-1]-std[0:nSpec-1], mn[0:nSpec-1]+std[0:nSpec-1]
	
	cgplot, mn[nSpec:2*nSpec-1], psym=4
	cgerrplot, mn[nSpec:2*nSpec-1]-std[nSpec:2*nSpec-1], mn[nSpec:2*nSpec-1]+std[nSpec:2*nSpec-1]
	
	cgplot, mn[2*nSpec:3*nSpec-1], psym=4
	cgerrplot, mn[2*nSpec:3*nSpec-1]-std[2*nSpec:3*nSpec-1], mn[2*nSpec:3*nSpec-1]+std[2*nSpec:3*nSpec-1]
	
	samples = dblarr(nl,length)
	temp = dblarr(nl)
	openr,2,'../test.spectra.txt'
	for i = 0, length-1 do begin
		readf,2,temp
		samples[*,i] = temp
	endfor
	close,2
	
	stop
	
 	!p.multi = 0
 		
	noise=mean(spec[1,*])/100.
	
	cgplot, spec[1,*], /xs; xran=[0,200]
	cgerrplot, spec[1,*]-noise, spec[1,*]+noise, col='red'
	for i = 0, length-1 do begin
		cgoplot, samples[*,i], col='blue'
	endfor

	!p.multi = 0
	
	chi2 = dblarr(100)
	std = dblarr(nl)	
	for i =0,99 do begin
		chi2[i]=total((samples[*,i]-spec[1,*])^2/noise^2)
	endfor
	for i =0,nl-1 do begin
		std[i] = stddev(samples[i,*])
	endfor
		
	print, std[0], stddev(t[0,*])
	

	stop
end
