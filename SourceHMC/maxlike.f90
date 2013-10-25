module maximumlikelihood_m
use params
use like_m, only : slikelihoodNoNormalization
use geodesiclm_m, only : maximumLikelihoodGeodesicLM
use maths, only : sigmoid
implicit none
	
contains

!-----------------------------------------------------------------------
! Do the actual nested sampling to estimate the evidence and the posterior distribution
!-----------------------------------------------------------------------
	subroutine getMaximumLikelihood
   integer :: i, loop
   integer :: nbin, mmax, ndata, iprint, k
	character(len=60) :: task, csave
	logical :: lsave(4)
	integer, allocatable :: nbd(:), iwa(:), isave(:)
	real(kind=8) :: factr, pgtol, f, f2
	real(kind=8), allocatable :: x(:), l(:), u(:), g(:), dsave(:), wa(:), x2(:), g2(:)

! Set the mask for this spectrum
		galaxy%nMask = count(galaxy%maskIndex(:,galaxy%whichComputing) == 1)
		allocate(galaxy%mask(galaxy%nMask))
		loop = 1
		do i = 1, galaxy%nPix
			if (galaxy%maskIndex(i,galaxy%whichComputing) == 1) then
				galaxy%mask(loop) = i
				loop = loop + 1
			endif
		enddo

! Estimate noise level : sigma = <spec> / SNR
		galaxy%noise = (sum(galaxy%spec(galaxy%mask,galaxy%whichComputing)) / galaxy%nMask) / galaxy%snr

		allocate(x(sdim),l(sdim),u(sdim),g(sdim))

! Boundaries
		l = -20.d0
		u = 20.d0
		do i = 1, sdim
			call random_number(f)
			x(i) = 1.d0 * f - 0.5d0
		enddo
		
! 		x = log((priorWeight%lower-1.d0/45.d0)/(1.d0/45.d0-priorWeight%upper))
! 		x(sdim-1) = log((priorVelocity%lower-0.d0)/(0.d0-priorVelocity%upper))
! 		x(sdim) = log((priorDispersion%lower-380.d0)/(380.d0-priorDispersion%upper))
				
! Test derivatives
  		do i = 1, sdim
  			x(i) = -2.d0 + 4.d0/sdim * i
  		enddo
!  				  		
!   		allocate(x2(sdim),g2(sdim))
!   		x2 = x  		
!   		call slikelihoodNoNormalization(x,f,g,.TRUE.)
!   		g2 = galaxy%synthGrad(galaxy%mask(100),:)
!   		galaxy%synth2 = galaxy%synth
!   		do i = 1, sdim
!   			x2 = x
!   			x2(i) = x2(i) + 1.d-5
!   			call slikelihoodNoNormalization(x2,f2,g,.TRUE.)
! !   			write(*,FMT='(I3,3(2X,E12.5))') i, (f2-f) / 1.d-5, g(i)
!   			
!    			print *, g2(i), (galaxy%synth(galaxy%mask(100)) - galaxy%synth2(galaxy%mask(100))) / 1.d-5
!   			
!   		enddo
!   		deallocate(x2,g2) 		
!    	stop

! SCALCG
! 		call maximumLikelihoodScalCG(sdim, x)

! DESCON
! 		call maximumLikelihoodDescon(sdim, x)

! LBFG
!  		call maximumLikelihoodLBFGSB(sdim, x, u, l)

! Geodesic Levenberg-Marquardt
  		call maximumLikelihoodGeodesicLM(sdim, x)

! Compute the spectrum at the final value of the parameters
 		call slikelihoodNoNormalization(x,f,g,.TRUE.)
		
		galaxy%trial = x
				
		deallocate(x,l,u,g)
								
!  		do i = 1, library%nSpec
! 			if (galaxy%trial(i) < 1.d-10) then
! 				galaxy%trial(i) = sigmoid(galaxy%trial(i) + 1.d-5)
! 			endif
! 			if (galaxy%trial(i) == 1.5d0) then
! 				galaxy%trial(i) = sigmoid(galaxy%trial(i) - 1.d-5)
! 			endif
! 		enddo

! Use the sigmoid transformation
!  		do i = 1, library%nSpec
! 			galaxy%trial(i) = sigmoid(galaxy%trial(i) + 1.d-5)
! 		enddo

		open(unit=15,file='bestFit.dat',action='write',status='replace')
  		if (galaxy%fixVLOS == 0) then
			write(15,*) library%nSpec, library%nSpec, library%nSpec
		else
			write(15,*) library%nSpec, 1, 1
		endif
		do i = 1, library%nSpec
			write(15,*) galaxy%trialWeight(i), galaxy%trialVelocity(i) * library%velScale, galaxy%trialDispersion(i) * library%velScale
		enddo
		write(15,*) size(galaxy%mask)
		do i = 1, size(galaxy%mask)
			write(15,*) galaxy%synth(galaxy%mask(i)), galaxy%spec(galaxy%mask(i),galaxy%whichComputing)
		enddo
		close(15)
				
		deallocate(galaxy%mask)
		
	end subroutine getMaximumLikelihood
	
end module maximumlikelihood_m