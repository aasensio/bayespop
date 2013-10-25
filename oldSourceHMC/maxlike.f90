module maximumlikelihood_m
use params
use like_m, only : slikelihoodNoNormalization
use l_bfgs_b, only : setulb
use scalcg_m, only : maximumLikelihoodScalCG
use descon_m, only : descon
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
		loop = 1		
		do i = 1, library%nSpec
			l(loop) = 0.d0
			u(loop) = 1.5d0
			x(loop) = 0.001d0
			loop = loop + 1
		enddo
		if (galaxy%fixVLOS == 0) then
			do i = 1, library%nSpec
				l(loop) = priorVelocity%lower
				u(loop) = priorVelocity%upper
				x(loop) = priorVelocity%mu
				loop = loop + 1
			enddo
			do i = 1, library%nSpec
				l(loop) = priorDispersion%lower
				u(loop) = priorDispersion%upper
				x(loop) = priorDispersion%mu
				loop = loop + 1
			enddo
		else
			l(loop) = priorVelocity%lower
			u(loop) = priorVelocity%upper
			x(loop) = priorVelocity%mu
			loop = loop + 1
			
			l(loop) = priorDispersion%lower
			u(loop) = priorDispersion%upper
			x(loop) = priorDispersion%mu
			loop = loop + 1
		endif
				
!   		x2 = x
!  		
!  		call slikelihoodNoNormalization(x,f,g,.TRUE.)
!  		do i = 1, sdim
!  			x2 = x
!  			x2(i) = x2(i) + 1.d-5
!  			call slikelihoodNoNormalization(x2,f2,g2,.TRUE.)
!  			write(*,FMT='(I3,3(2X,E12.5))') i, (f2-f) / 1.d-5, g(i), x(i)
!  		enddo
! 		
!  		stop

! SCALCG
! 		call maximumLikelihoodScalCG(sdim, x)

! DESCON
! 		call maximumLikelihoodDescon(sdim, x)

! LBFG
! 		call maximumLikelihoodLBFGSB(sdim, x, u, l)

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
 		do i = 1, library%nSpec
			galaxy%trial(i) = sigmoid(galaxy%trial(i) + 1.d-5)
		enddo

		open(unit=15,file='bestFit.dat',action='write',status='replace')
  		if (galaxy%fixVLOS == 0) then
			write(15,*) library%nSpec, library%nSpec, library%nSpec
		else
			write(15,*) library%nSpec, 1, 1
		endif
		write(15,*) galaxy%trial
		write(15,*) size(galaxy%mask)
		do i = 1, size(galaxy%mask)
			write(15,*) galaxy%synth(galaxy%mask(i)), galaxy%spec(galaxy%mask(i),galaxy%whichComputing)
		enddo
		close(15)
				
		deallocate(galaxy%mask)
		
	end subroutine getMaximumLikelihood
	
end module maximumlikelihood_m