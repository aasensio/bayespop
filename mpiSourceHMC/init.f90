! Include file for example nested sampler program 'Gaussian Rings'

module init
use params
implicit none

contains

!-----------------------------------------------------------------
! Read the problem configuration
!-----------------------------------------------------------------
	subroutine setProblem(fileIn, whichGalaxy)
	integer :: i, j, whichGalaxy
	character(len=120) :: fileIn
	real(kind=8) :: l0, par(2)
				
		open(unit=12,file='preproc_data/'//trim(adjustl(fileIn))//".input",action='read',status='old',form='unformatted')
				
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
												
		read(12) library%nPixAdapted
		read(12) library%nSpec
		read(12) library%velScale
		
		galaxy%nPix = library%nPixAdapted
		
		if (associated(library%specAdapted)) deallocate(library%specAdapted)
		if (associated(library%age)) deallocate(library%age)
		if (associated(library%metallicity)) deallocate(library%metallicity)
		if (associated(library%imfSlope)) deallocate(library%imfSlope)
		if (associated(library%mgfe)) deallocate(library%mgfe)
		
		allocate(library%specAdapted(library%nPixAdapted,library%nSpec))
		allocate(library%age(library%nSpec))
		allocate(library%metallicity(library%nSpec))
		allocate(library%imfSlope(library%nSpec))
		allocate(library%mgfe(library%nSpec))
		
		read(12) library%specAdapted
		read(12) library%age
		read(12) library%metallicity
		read(12) library%imfSlope
		read(12) library%mgfe
				
		read(12) galaxy%snr
		read(12) galaxy%nSpec
		
		
		if (associated(galaxy%spec)) deallocate(galaxy%spec)
		if (associated(galaxy%maskIndex)) deallocate(galaxy%maskIndex)
		if (associated(galaxy%synth)) deallocate(galaxy%synth)
		if (associated(galaxy%synth2)) deallocate(galaxy%synth2)
		if (associated(galaxy%convolvedSpectrum)) deallocate(galaxy%convolvedSpectrum)
		
		allocate(galaxy%specRead(galaxy%nPix,galaxy%nSpec))
		allocate(galaxy%spec(galaxy%nPix))
		allocate(galaxy%maskIndex(galaxy%nPix))
		allocate(galaxy%synth(galaxy%nPix))
		allocate(galaxy%synth2(galaxy%nPix))
		allocate(galaxy%convolvedSpectrum(galaxy%nPix))
					
! Read all spectra and choose only the one we want
		read(12) galaxy%specRead
		galaxy%spec = galaxy%specRead(:,whichGalaxy)
		deallocate(galaxy%specRead)
		read(12) galaxy%maskIndex		
		read(12) galaxy%fixVLOS				
				
		close(12)		

! ! Fix all velocities to a common value
! 		if (galaxy%fixVLOS == 0) then
! 			sdim = 3*library%nSpec
! 		else
! 			sdim = library%nSpec+2
! 		endif
! 		
! 		nest_nPar = sdim
! 
! ! Allocate memory for gradient of spectrum used in the calculation of the derivative of the log-posterior
! 		if (associated(galaxy%synthGrad)) deallocate(galaxy%synthGrad)
! 		allocate(galaxy%synthGrad(galaxy%nPix,sdim))
! 
! 		if (allocated(map_pars)) deallocate(map_pars)
! 		allocate(map_pars(sdim))		
! 
! ! For periodic boundary conditions
! 		if (allocated(pWrap)) deallocate(pWrap)
! 		allocate(pWrap(sdim))
! 		pWrap = 0
! 
! 		if (associated(galaxy%trialWeight)) deallocate(galaxy%trialWeight)
! 		if (associated(galaxy%trialVelocity)) deallocate(galaxy%trialVelocity)
! 		if (associated(galaxy%trialDispersion)) deallocate(galaxy%trialDispersion)
! 		
! 		allocate(galaxy%trialWeight(library%nSpec))
! 		allocate(galaxy%trialVelocity(library%nSpec))
! 		allocate(galaxy%trialDispersion(library%nSpec))
! 				
! 		if (associated(galaxy%trial)) deallocate(galaxy%trial)
! 		allocate(galaxy%trial(sdim))
! 
! 		maxslhood = -1.d100
				
	end subroutine setProblem

end module init