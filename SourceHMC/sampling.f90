module sampling_m
use params
use like_m, only : slikelihoodNoNormalization
use maths, only : invSigmoid, sigmoid
implicit none

contains

!-----------------------------------------------------------------------
! Do the actual nested sampling to estimate the evidence and the posterior distribution
!-----------------------------------------------------------------------
	subroutine hmcSample
   integer :: nclusters, context !total number of clusters found
   integer :: maxNode !variables used by the posterior routine
   integer :: i, j, loop, nSamplesEstimationVariance
   real(kind=8) :: f, mean
   real(kind=8), allocatable :: samples(:,:)


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
		
! Initial parameters
		allocate(st(sdim))
		allocate(stepSize(sdim))

! Initial values obtained with a optimization routine
! These values are obtained 
		st = galaxy%trial
		
		call testDerivatives(st)
		
! Size of the steps and initial conditions
! Transform from the initial conditions to the sigmoid functions

		stepSize(1:library%nSpec) = 1.d0						
		
		if (galaxy%fixVLOS == 0) then
			stepSize(library%nSpec+1:2*library%nSpec) = 1.d0
			stepSize(2*library%nSpec+1:3*library%nSpec) = 1.d0			
		else
		
			stepSize(library%nSpec+1) = 1.d0			
			stepSize(library%nSpec+2) = 1.d0			
		endif
				
			
		scaleFactor = 1.d0
 		seed = 1234
		fbInt = 10
		maxStep = 10
		resume = 0
		flPfx = 'test'
		
		stepSize = stepSize / maxStep
			
! First 100 iterations to estimate the stepSize
		nSamplesEstimationVariance = 500
		open(unit=20,file= (trim(flPfx)//".extract.txt"),recl=1500000)
		open(unit=21,file= (trim(flPfx)//".spectra.txt"),recl=1500000)
				
		call run_guided_hmc(sdim, st, scaleFactor, maxStep, stepSize, flPfx(1:len_trim(flPfx)), seed, 0, &
			fbInt, negLogPosterior, writeHMCProcess, 0, nSamplesEstimationVariance)
			
		close(20)
		close(21)
		
		allocate(samples(sdim,nSamplesEstimationVariance))
		
! Estimate the variance of the variables		
		open(unit=20,file= (trim(flPfx)//".extract.txt"),recl=1500000,status='old',action='read')
		do i = 1, nSamplesEstimationVariance
			read(20,*) (samples(j,i),j=1,sdim)			
		enddo
		close(20)
		do i = 1, sdim
			mean = sum(samples(i,:)) / (1.d0*nSamplesEstimationVariance)
			stepSize(i) = sqrt(sum((samples(i,:)-mean)**2) / (nSamplesEstimationVariance-1.d0))			
		enddo
		
		deallocate(samples)
		
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
	integer i	
		
	do i=1,ndim-1		
		if (i <= library%nSpec) then
			write(20,'(E18.10,a)',advance='no') sigmoid(x(i), 0.d0, 1.d0) !1.d0/(1.d0+exp(-x(i))), ' '
		else
			write(20,'(E18.10,a)',advance='no') x(i), ' '
		endif
	enddo
		
	write(20,'(E18.10)') x(ndim)
	
	write(21,*) galaxy%synth(galaxy%mask)

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
