! Include file for example nested sampler program 'Gaussian Rings'

module init
use params
implicit none

contains

!-----------------------------------------------------------------
! Read the problem configuration
!-----------------------------------------------------------------
	subroutine setProblem
	integer :: i, j
	real(kind=8) :: l0, par(2)

		open(unit=12,file='internalConf.input',action='read',status='old',form='unformatted')		
		read(12) nest_root
		
		allocate(priors(7))

! ! Prior for weights
! 		read(12) priorWeight%typ
! 		select case (priorWeight%typ)
! 			case(0)
! 				read(12) priorWeight%lower, priorWeight%upper
! 				read(12) priorWeight%mu
! 			case(1)
! 				read(12) priorWeight%lower, priorWeight%upper
! 				read(12) priorWeight%mu, priorWeight%sigma
! 		end select		

! Prior for velocity, dispersion, age, metallicity, imf_slope, mgfe
		do i = 1, 6
			read(12) priors(i)%typ
			select case (priors(i)%typ)
				case(0)
					read(12) priors(i)%lower, priors(i)%upper
					read(12) priors(i)%mu
				case(1)
					read(12) priors(i)%lower, priors(i)%upper
					read(12) priors(i)%mu, priors(i)%sigma
			end select
		enddo
						
		print *, 'Priors'
		do i = 1, 6
			print *, variableNames(i), priors(i)%lower, priors(i)%upper
		enddo
				
		read(12) library%nSpec
		read(12) library%nPixAdapted
		read(12) library%velScale
		
		print *, library%nSpec, library%nPixAdapted, library%velScale
		stop

		galaxy%nPix = library%nPixAdapted
		
		allocate(library%specAdapted(library%nPixAdapted,library%nSpec))
		read(12) library%specAdapted

		read(12) galaxy%snr
		read(12) galaxy%nSpec
		
		allocate(galaxy%spec(galaxy%nPix,galaxy%nSpec))
		allocate(galaxy%maskIndex(galaxy%nPix,galaxy%nSpec))
		allocate(galaxy%synth(galaxy%nPix))
		allocate(galaxy%synth2(galaxy%nPix))
		allocate(galaxy%convolvedSpectrum(galaxy%nPix))
		
		read(12) galaxy%spec		
		read(12) galaxy%maskIndex
		read(12) galaxy%fixVLOS
				
		close(12)

! Fix all velocities to a common value
		if (galaxy%fixVLOS == 0) then
			sdim = 3*library%nSpec
		else
			sdim = library%nSpec+2
		endif
		
		nest_nPar = sdim

! Allocate memory for gradient of spectrum used in the calculation of the derivative of the log-posterior
		allocate(galaxy%synthGrad(galaxy%nPix,sdim))

		allocate(map_pars(sdim))		

! For periodic boundary conditions
		allocate(pWrap(sdim))
		pWrap = 0

		allocate(galaxy%trialWeight(library%nSpec))
		allocate(galaxy%trialVelocity(library%nSpec))
		allocate(galaxy%trialDispersion(library%nSpec))
				
		allocate(galaxy%trial(sdim))

		maxslhood = -1.d100

	end subroutine setProblem

! The following routines might be part of the code in the future, when it
! reads the library and observation by itself
! !-------------------------------------------------------------
! ! Split a string into its fields
! !-------------------------------------------------------------
! 	subroutine splitString(string, n, words)
! 	character(len=120) :: string, words(n)
! 	character :: chr
! 	integer :: i, j, n, pos1, pos2
! 
! 		pos1 = 1
! 		pos2 = 2
! 		do i = 1, n
! 			chr = string(pos2:pos2)
! 			do while (chr /= ' ')
! 				pos2 = pos2 + 1
! 				chr = string(pos2:pos2)
! 			enddo
! 			words(i) = string(pos1:pos2)
! 
! 			do while (chr == ' ')
! 				pos2 = pos2 + 1
! 				chr = string(pos2:pos2)
! 			enddo
! 			
! 			pos1 = pos2			
! 		enddo
! 
! 	end subroutine splitString
! 
! !-------------------------------------------------------------
! ! Read spectrum from library
! !-------------------------------------------------------------
! 	subroutine readSpectrumFromLibrary(file, indexLibrary)
! 	character(len=120) :: file, comment
! 	integer :: indexLibrary, status, unit, naxes, blocksize, i, hdutype
! 	real(kind=8) :: dLambda, zero	
! 	logical :: anynull
! 
! 		print *, 'Reading : ', trim(adjustl(file))
! 		
! 		status = 0
!       call ftgiou(unit, status)
!       if (status /= 0)  then
! 	      call ftrprt('STDOUT',status)
! 	      stop
!       endif
! 
!       call ftopen(unit,file,1,blocksize,status) !open FITS file
! 
! 		if (status /= 0)  then
! 			call ftrprt('STDOUT',status)
! 			stop
! 		endif
! 
! ! Read some information from the header
! 		call ftgkyj(unit,'NAXIS1',library%nPix,comment,status)
! 		call ftgkyd(unit,'CDELT1',dLambda,comment,status)
! 		call ftgkyd(unit,'CRVAL1',zero,comment,status)
! 
! ! Generate the wavelength axis
! 		if (.not. associated(library%lambda)) then
! 			allocate(library%lambda(library%nPix))
! 			do i = 1, library%nPix
! 				library%lambda(i) = dLambda*(i-1.d0) + zero
! 			enddo
! 		endif
! 
! ! Allocate memory for library if not done
! 		if (.not. associated(library%spec)) then
! 			allocate(library%spec(library%nPix,library%nSpec))
! 		endif
! 		
! ! Do the actual reading		
!  		call ftgpvd(unit,1,1,library%nPix,0.0,library%spec(:,indexLibrary),anynull,status)
! 
! 		call ftclos(unit, status)
! 		call ftfiou(unit, status)
! 
! 		library%velscale = lightSpeed_c * dLambda / (sum(library%lambda) / library%nPix)
! 		
! 	end subroutine readSpectrumFromLibrary
! 
! !-------------------------------------------------------------
! ! Read spectrum from observations
! !-------------------------------------------------------------
! 	subroutine readSpectrumFromObservation(file, indexObservation)
! 	character(len=120) :: file, comment
! 	integer :: indexObservation, status, unit, naxes, blocksize, i, hdutype
! 	real(kind=8) :: dLambda, zero
! 	logical :: anynull
! 
! 		print *, 'Reading : ', trim(adjustl(file))
! 
! 		status = 0
!       call ftgiou(unit, status)
!       if (status /= 0)  then
! 	      call ftrprt('STDOUT',status)
! 	      stop
!       endif
! 
!       call ftopen(unit,file,1,blocksize,status) !open FITS file
! 
! 		if (status /= 0)  then
! 			call ftrprt('STDOUT',status)
! 			stop
! 		endif
! 
! ! Read some information from the header
! 		call ftgkyj(unit,'NAXIS1',galaxy%nPix,comment,status)
! 		call ftgkyd(unit,'CDELT1',dLambda,comment,status)
! 		call ftgkyd(unit,'CRVAL1',zero,comment,status)
! 
! ! Generate the wavelength axis
! 		if (.not. associated(galaxy%lambda)) then
! 			allocate(library%lambda(library%nPix))
! 			do i = 1, library%nPix
! 				library%lambda(i) = dLambda*(i-1.d0) + zero
! 			enddo
! 		endif
! 
! ! Allocate memory for library if not done
! 		if (.not. associated(galaxy%spec)) then
! 			allocate(galaxy%spec(library%nPix,library%nSpec))
! 		endif
! 
! ! Do the actual reading
!  		call ftgpvd(unit,1,1,galaxy%nPix,0.0,galaxy%spec(:,indexObservation),anynull,status)
! 
! 		call ftclos(unit, status)
! 		call ftfiou(unit, status)
! 
! 		galaxy%velscale = lightSpeed_c * dLambda / (sum(galaxy%lambda) / galaxy%nPix)
! 
! 	end subroutine readSpectrumFromObservation
! 
! !-------------------------------------------------------------
! ! Read spectrum from library
! !-------------------------------------------------------------
! 	subroutine adaptSpectrumFromLibrary(indexLibrary)
! 	integer :: indexLibrary
! 
! 		
! 
! 	end subroutine adaptSpectrumFromLibrary
! 	
! !-------------------------------------------------------------
! ! Read all data, including model spectra and observation
! !-------------------------------------------------------------
! 	subroutine initData
! 	integer :: i, pos, eof, nTotal
! 	character(len=120) :: dummy, a, f
! 	character(len=120), allocatable :: res(:)
! 	real(kind=8) :: b, c, d, e
! 
! ! Read the observations
! 		print *, 'Reading observations : ', trim(adjustl(galaxy%filename))
! 		open(unit=12,file=galaxy%filename,action='read',status='old')
! 		read(12,*)
! ! Count the number of uncommented models
! 		eof = 0
! 		nTotal = 0
! 		do while(eof == 0)
! 			read(12,*,iostat=eof) dummy
! 			if (eof == 0) then
! 				if (dummy(1:1) /= '#') then
! 					nTotal = nTotal + 1
! 				endif
! 			endif
! 		enddo
! 		close(12)
! 
! 		print *, 'Number of selected observations : ', nTotal
! 
! 		galaxy%nSpec = nTotal
! 
! 		allocate(res(7))
! 		open(unit=12,file=galaxy%filename,action='read',status='old')
! 		read(12,*)
! 		nTotal = 0
! 		eof = 0
! 		do while(eof == 0)
! 			read(12,FMT='(A)',iostat=eof) dummy
! 			if (eof == 0) then
! 				if (dummy(1:1) /= '#') then
! 					nTotal = nTotal + 1
! 					call splitString(dummy, 7, res)
! 					call readSpectrumFromObservation(res(1), nTotal)
! 				endif
! 			endif
! 		enddo
! 		close(12)
! 		deallocate(res)
! 		
! 		
! ! Read the library
! 		
! 		print *, 'Reading library : ', trim(adjustl(library%filename))
! 		open(unit=12,file=library%filename,action='read',status='old')
! 		read(12,FMT='(A)') dummy
! 		pos = index(dummy, ':')
! 		library%name = dummy(pos+1:len_trim(dummy))
! 		print *, 'Library name : ', trim(adjustl(library%name))
! 
! 		read(12,FMT='(A)') dummy
! 		pos = index(dummy, ':')
! 		read(dummy(pos+1:len_trim(dummy)),FMT='(F)') library%fwhm
! 		print *, 'Library FWHM : ', library%fwhm
! 
! ! Count the number of uncommented models
! 		read(12,*)
! 		eof = 0
! 		nTotal = 0
! 		do while(eof == 0)
! 			read(12,*,iostat=eof) dummy
! 			if (eof == 0) then
! 				if (dummy(1:1) /= '#') then
! 					nTotal = nTotal + 1
! 				endif
! 			endif
! 		enddo
! 		close(12)
! 
! 		print *, 'Number of selected models : ', nTotal
! 
! 		library%nSpec = nTotal
! 		
! 		open(unit=12,file=library%filename,action='read',status='old')
! 		read(12,*)
! 		read(12,*)
! 		read(12,*)
! 		nTotal = 0
! 		eof = 0
! 		allocate(res(6))
! 		do while(eof == 0)
! 			read(12,FMT='(A)',iostat=eof) dummy
! 			if (eof == 0) then
! 				if (dummy(1:1) /= '#') then
! 					nTotal = nTotal + 1
! 					call splitString(dummy, 6, res)
! 					call readSpectrumFromLibrary(res(6), nTotal)
! 					call adaptSpectrumFromLibrary(nTotal)
! 				endif
! 			endif
! 		enddo
! 		close(12)
! 		deallocate(res)
! 
! 		stop
! 
! 
! 	end subroutine initData


end module init