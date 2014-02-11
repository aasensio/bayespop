module sampling_m
use params
use like_m, only : slikelihoodNoNormalization
use maths, only : invSigmoid, sigmoid
implicit none

contains

!-----------------------------------------------------------------------
! Transform samples back to the sigmoid parameters
!-----------------------------------------------------------------------
	function getBackToSigmoid(x)
	real(kind=8) :: x(:)
	real(kind=8) :: getBackToSigmoid(size(x))
	
		getBackToSigmoid(1:library%nSpec) = invsigmoid(x(1:library%nSpec), 0.d0, 1.d0)
		
		if (galaxy%fixVLOS == 0) then
			getBackToSigmoid(library%nSpec+1:2*library%nSpec) = invsigmoid(x(library%nSpec+1:2*library%nSpec), priors(1)%lower, priors(1)%upper)
			getBackToSigmoid(2*library%nSpec+1:3*library%nSpec) = invsigmoid(x(2*library%nSpec+1:3*library%nSpec), priors(2)%lower, priors(2)%upper)
		else
			getBackToSigmoid(library%nSpec+1) = invsigmoid(x(library%nSpec+1), priors(1)%lower, priors(1)%upper)
			getBackToSigmoid(library%nSpec+2) = invsigmoid(x(library%nSpec+2), priors(2)%lower, priors(2)%upper)
		endif
	end function getBackToSigmoid

!-----------------------------------------------------------------------
! Do the actual nested sampling to estimate the evidence and the posterior distribution
!-----------------------------------------------------------------------
	subroutine hmcSample
   integer :: nclusters, context !total number of clusters found
   integer :: maxNode !variables used by the posterior routine
   integer :: i, j, loop, nSamplesEstimationVariance
   real(kind=8) :: f, mean
   real(kind=8), allocatable :: samples(:,:), temp(:)


! Set the mask for this spectrum
		galaxy%nMask = count(galaxy%maskIndex == 1)
		
		if (associated(galaxy%mask)) deallocate(galaxy%mask)				
		allocate(galaxy%mask(galaxy%nMask))
		
		loop = 1
		do i = 1, galaxy%nPix
			if (galaxy%maskIndex(i) == 1) then
				galaxy%mask(loop) = i
				loop = loop + 1
			endif
		enddo
				
! Estimate noise level : sigma = <spec> / SNR		
		galaxy%noise = (sum(galaxy%spec(galaxy%mask)) / galaxy%nMask) / galaxy%snr		
		
! Initial parameters
		if (allocated(st)) deallocate(st)
		if (allocated(stepSize)) deallocate(stepSize)
		if (allocated(savePars)) deallocate(savePars)
		
		allocate(st(sdim))
		allocate(stepSize(sdim))
		allocate(savePars(sdim))

! Initial values obtained with a optimization routine
! These values are obtained 
		st = galaxy%trial
		
! 		call testDerivatives(st)
		
! Size of the steps and initial conditions
! Transform from the initial conditions to the sigmoid functions
		stepSize(1:library%nSpec) = 0.01d0
		
		if (galaxy%fixVLOS == 0) then
			stepSize(library%nSpec+1:2*library%nSpec) = 0.01d0
			stepSize(2*library%nSpec+1:3*library%nSpec) = 0.01d0			
		else
		
			stepSize(library%nSpec+1) = 0.01d0
			stepSize(library%nSpec+2) = 0.01d0
		endif
									
		scaleFactor = 1.d0
 		seed = 1234
		fbInt = 10
		maxStep = 100
		resume = 0
		flPfx = 'temp/test'//trim(adjustl(myrankStr))
		
! 		stepSize = stepSize / maxStep
			
! First 100 iterations to estimate the stepSize
		nSamplesEstimationVariance = 200
		open(unit=20,file= (trim(flPfx)//".samples"),form='unformatted',action='write',access='stream')
		open(unit=21,file= (trim(flPfx)//".spectra"),form='unformatted',action='write',access='stream')
		write(20) nSamplesEstimationVariance
		write(21) galaxy%nPix
		
! 		do i = 1, nSamplesEstimationVariance
! 			call writeHMCProcess(sdim,st,1.d0,st)
! 		enddo
		
   		call run_guided_hmc(sdim, st, scaleFactor, maxStep, stepSize, flPfx(1:len_trim(flPfx)), seed, 0, &
   			fbInt, negLogPosterior, writeHMCProcess, 0, nSamplesEstimationVariance)
 						
		close(20)
		close(21)
				
! 		if (allocated(samples)) deallocate(samples)
! 		if (allocated(temp)) deallocate(temp)				
! 		allocate(samples(sdim,nSamplesEstimationVariance))
! 		allocate(temp(sdim))
! 		
! ! Estimate the variance of the variables		
! 		open(unit=20,file= (trim(flPfx)//".samples"),form='unformatted',action='read',access='stream')
! 		read(20) nSamplesEstimationVariance
! 		do i = 1, nSamplesEstimationVariance
! 			read(20) temp		
! 			samples(:,i) = getBackToSigmoid(temp)
! 		enddo
! 		close(20)
! 		do i = 1, sdim
! 			mean = sum(samples(i,:)) / (1.d0*nSamplesEstimationVariance)
! 			st(i) = mean			
! 			stepSize(i) = sqrt(sum((samples(i,:)-mean)**2) / (nSamplesEstimationVariance-1.d0))			
! 		enddo
! ! 		stepSize = stepSize / maxStep
! 						
! 		nSamplesEstimationVariance = 500
! 		open(unit=20,file= (trim(flPfx)//".samples"),form='unformatted',action='write',access='stream')
! 		open(unit=21,file= (trim(flPfx)//".spectra"),form='unformatted',action='write',access='stream')
! 		write(20) nSamplesEstimationVariance
! 		write(21) galaxy%nPix
! 		
!   		call run_guided_hmc(sdim, st, scaleFactor, maxStep, stepSize, flPfx(1:len_trim(flPfx)), seed, 0, &
!   			fbInt, negLogPosterior, writeHMCProcess, 0, nSamplesEstimationVariance)
!  						
! 		close(20)
! 		close(21)

	end subroutine hmcSample

	subroutine negLogPosterior(ndim,x,v,g)
	integer ndim
	real(kind=8), dimension(ndim) :: x
   real(kind=8), dimension(ndim) :: g
   real(kind=8) :: v
	integer i
				
 		call slikelihoodNoNormalization(x,v,g,.TRUE.)
		v = -v		
		g = -g
  
	end subroutine negLogPosterior

!
! A subroutie to write the extract file
! I have assumed that the unit=20 is opened for
! writing (append) earlier. In general only write
! those parametes which are estimated. The files
! can be really big depending on the dimensionality
!
	subroutine writeHMCProcess(ndim,x,v,g)
	integer ndim
	real(kind=8), dimension(ndim) :: x, x2
   real(kind=8), dimension(ndim) :: g
   real(kind=8) :: v
	integer i	, loop
	
	loop = 1
	do i = 1, library%nSpec
		savePars(loop) = galaxy%trialWeight(i)		
		loop = loop + 1
	enddo
	if (galaxy%fixVLOS == 0) then
		do i = 1, library%nSpec
			savePars(loop) = galaxy%trialVelocity(i)
			loop = loop + 1
		enddo
	else
		savePars(loop) = galaxy%trialVelocity(1)
		loop = loop + 1
	endif
	if (galaxy%fixVLOS == 0) then
		do i = 1, library%nSpec
			savePars(loop) = galaxy%trialDispersion(i)
			loop = loop + 1
		enddo
	else
		savePars(loop) = galaxy%trialDispersion(1)
		loop = loop + 1
	endif
		
	write(20) savePars	
	write(21) galaxy%synth	

end subroutine writeHMCProcess

!------------------------------------------------------------------
! Test derivatives
!------------------------------------------------------------------
	subroutine testDerivatives(pars)
	real(kind=8) :: pars(:), parsNew(size(pars)), logP, logPNew, logPGradient(size(pars)), logPGradientNew(size(pars))
	integer :: n, i
	
		n = size(pars)
		
		call negLogPosterior(n,pars,logP,logPGradient)
		
		do i = 1, n
			parsNew = pars
			parsNew(i) = parsNew(i) + 1.d-5
			call negLogPosterior(n,parsNew,logPNew,logPGradientNew)
			print *, logPGradient(i), (logPNew-logP) / 1.d-5			
		enddo
	end subroutine testDerivatives


end module sampling_m
