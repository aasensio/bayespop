module like_m
use params
use maths, only : convolution, sigmoid, diffSigmoid
implicit none
      
contains
	
!-----------------------------------------------------------------------
! Likelihood function
!-----------------------------------------------------------------------
	subroutine slikelihoodNoNormalization(trial,slhood,gradient,maxLikelihood)
	logical :: maxLikelihood
	real(kind=8) :: trial(:), slhood, logPrior, logLikelihood, gradient(sdim)
	real(kind=8) :: temp(sdim),dist,loclik, sigma, xi, L, maxDispersion, sumKernel, covar(sdim,sdim), sumKernel2
	integer :: i,j, nPixKernel, midPix, activeWeight(library%nSpec)

! Take also into account the priors
		logPrior = 0.d0
		gradient = 0.d0
						
! ********* Prior for the weights
		galaxy%trialWeight = sigmoid(trial(1:library%nSpec), priorWeight%lower, priorWeight%upper)
		
! Gaussian prior
		if (priorWeight%typ == 1) then
			logPrior = logPrior - sum((galaxy%trialWeight - priorWeight%mu)**2 / (2.d0*priorWeight%sigma**2))
			gradient(1:library%nSpec) = gradient(1:library%nSpec) - (galaxy%trialWeight - priorWeight%mu) / priorWeight%sigma**2
		endif

! Dirac prior
		if (priorWeight%typ == 2) then
			galaxy%trialWeight = priorWeight%mu
		endif
		

!******************
! If we do not fix the LOS velocity of each component
!******************
		if (galaxy%fixVLOS == 0) then

! ********* Prior for the velocity
			galaxy%trialVelocity = sigmoid(trial(library%nSpec+1:2*library%nSpec), priorVelocity%lower, priorVelocity%upper)

! Gaussian prior
			if (priorVelocity%typ == 1) then
				logPrior = logPrior - sum((galaxy%trialVelocity - priorVelocity%mu)**2 / (2.d0*priorVelocity%sigma**2))
				gradient(library%nSpec+1:2*library%nSpec) = gradient(library%nSpec+1:2*library%nSpec) - (galaxy%trialVelocity - priorVelocity%mu) / priorVelocity%sigma**2
			endif

! Dirac prior
			if (priorVelocity%typ == 2) then
				galaxy%trialVelocity = priorVelocity%mu
			endif
			
! ********* Prior for the dispersion
			galaxy%trialDispersion = sigmoid(trial(2*library%nSpec+1:3*library%nSpec), priorDispersion%lower, priorDispersion%upper)

! Gaussian prior
			if (priorDispersion%typ == 1) then
				logPrior = logPrior - sum((galaxy%trialDispersion - priorDispersion%mu)**2 / (2.d0*priorDispersion%sigma**2))
				gradient(2*library%nSpec+1:3*library%nSpec) = gradient(2*library%nSpec+1:3*library%nSpec) - (galaxy%trialDispersion - priorDispersion%mu) / priorDispersion%sigma**2
			endif

! Dirac prior
			if (priorDispersion%typ == 2) then
				galaxy%trialDispersion = priorDispersion%mu
			endif
			
		else

!******************
! If we fix the LOS velocity of each component, fill all the trialVelocity array with the same value
!******************			
! ********* Prior for the velocity
			galaxy%trialVelocity = sigmoid(trial(library%nSpec+1), priorVelocity%lower, priorVelocity%upper)

! Gaussian prior
			if (priorVelocity%typ == 1) then
				logPrior = logPrior - (galaxy%trialVelocity(1) - priorVelocity%mu)**2 / (2.d0*priorVelocity%sigma**2)
				gradient(library%nSpec+1) = gradient(library%nSpec+1) - (galaxy%trialVelocity(1) - priorVelocity%mu) / priorVelocity%sigma**2
			endif

! Dirac prior
			if (priorVelocity%typ == 2) then
				galaxy%trialVelocity = priorVelocity%mu
			endif
			
! ********* Prior for the dispersion
			galaxy%trialDispersion = sigmoid(trial(library%nSpec+2), priorDispersion%lower, priorDispersion%upper)
					
! Gaussian prior
			if (priorDispersion%typ == 1) then
				logPrior = logPrior - (galaxy%trialDispersion(1) - priorDispersion%mu)**2 / (2.d0*priorDispersion%sigma**2)
				gradient(library%nSpec+2) = gradient(library%nSpec+2) - (galaxy%trialDispersion(1) - priorDispersion%mu) / priorDispersion%sigma**2
			endif

! Dirac prior
			if (priorDispersion%typ == 2) then
				galaxy%trialDispersion = priorDispersion%mu
			endif
					
		endif									
		
! If the weight of the component is smaller than 0.001, set its weight to zero
		activeWeight = 1
! 		if (sigmoidTransform) then
! 			do i = 1, library%nSpec
! 				if (galaxy%trialWeight(i) < 0.001d0) then
! 					activeWeight(i) = 0
! 					trial(i) = -50.d0
! 					trial(i+library%nSpec) = 0.d0
! 					trial(i+2*library%nSpec) = 200.d0
! 				endif
! 			enddo
! 		endif

! **************************
! ********** DATA LIKELIHOOD
! **************************		
				
! Normalize the velocity and dispersion to a global scale
		galaxy%trialDispersion = galaxy%trialDispersion / library%velScale
		galaxy%trialVelocity = galaxy%trialVelocity / library%velScale

		maxDispersion = maxval(galaxy%trialDispersion)

! Find the number of pixels of the kernel
		nPixKernel = 2*ceiling(4*maxDispersion)+1

! Do something more efficient that allocating memory for each iteration
		allocate(galaxy%xKernel(nPixKernel))
		allocate(galaxy%yKernel(nPixKernel))		
		allocate(galaxy%yKernelDerivative(nPixKernel))
		
		galaxy%synth = 0.d0		
		galaxy%synthGrad = 0.d0
		midPix = nPixKernel / 2
		
! Compute the kernel, do the convolution, add the component and compute the derivatives
		do i = 1, library%nSpec

			if (activeWeight(i) == 1) then
				do j = 1, nPixKernel
					galaxy%xKernel(j) = (j-1 - (midPix+galaxy%trialVelocity(i))) / galaxy%trialDispersion(i)
				enddo
				galaxy%yKernel = exp(-0.5d0*galaxy%xKernel**2)
				sumKernel = sum(galaxy%yKernel)				

				galaxy%convolvedSpectrum = convolution(library%specAdapted(:,i), galaxy%yKernel / sumKernel)
				galaxy%synth = galaxy%synth + activeWeight(i) * galaxy%trialWeight(i) * galaxy%convolvedSpectrum
												

! Term used to compute the derivatives with respect to the weights
				galaxy%synthGrad(:,i) = galaxy%convolvedSpectrum

				if (galaxy%fixVLOS == 0) then

! Term used to compute the derivatives with respect to the velocities
					galaxy%yKernelDerivative = -(-galaxy%xKernel / galaxy%trialDispersion(i) * galaxy%yKernel * sumKernel - galaxy%yKernel * &
						sum(-galaxy%xKernel / galaxy%trialDispersion(i) * galaxy%yKernel)) / sumKernel**2

					galaxy%synthGrad(:,i+library%nSpec) = activeWeight(i) * galaxy%trialWeight(i) * &
						convolution(library%specAdapted(:,i), galaxy%yKernelDerivative ) / library%velScale

! Term used to compute the derivatives with respect to the dispersion
					galaxy%yKernelDerivative = (galaxy%xKernel**2 / galaxy%trialDispersion(i) * galaxy%yKernel * sumKernel - galaxy%yKernel * &
						sum(galaxy%xKernel**2 / galaxy%trialDispersion(i) * galaxy%yKernel)) / sumKernel**2

					galaxy%synthGrad(:,i+2*library%nSpec) = activeWeight(i) * galaxy%trialWeight(i) * &
						convolution(library%specAdapted(:,i), galaxy%yKernelDerivative ) / library%velScale
																
				else
! Compute the derivative of the kernel 
					galaxy%yKernelDerivative = -(-galaxy%xKernel / galaxy%trialDispersion(i) * galaxy%yKernel * sumKernel - galaxy%yKernel * &
						sum(-galaxy%xKernel / galaxy%trialDispersion(i) * galaxy%yKernel)) / sumKernel**2
							
					galaxy%synthGrad(:,1+library%nSpec) = galaxy%synthGrad(:,1+library%nSpec) + activeWeight(i) * galaxy%trialWeight(i) * &
						convolution(library%specAdapted(:,i), galaxy%yKernelDerivative ) / library%velScale

! Term used to compute the derivatives with respect to the dispersion
					galaxy%yKernelDerivative = (galaxy%xKernel**2 / galaxy%trialDispersion(i) * galaxy%yKernel * sumKernel - galaxy%yKernel * &
						sum(galaxy%xKernel**2 / galaxy%trialDispersion(i) * galaxy%yKernel)) / sumKernel**2
					
					galaxy%synthGrad(:,2+library%nSpec) = galaxy%synthGrad(:,2+library%nSpec) + activeWeight(i) * galaxy%trialWeight(i) * &
						convolution(library%specAdapted(:,i), galaxy%yKernelDerivative ) / library%velScale
				endif
												
			endif
			
		enddo
							
		deallocate(galaxy%xKernel)
		deallocate(galaxy%yKernel)
		deallocate(galaxy%yKernelDerivative)
		
! Compute the log-likelihood
		logLikelihood = -0.5d0 * sum( (galaxy%synth(galaxy%mask)-galaxy%spec(galaxy%mask,galaxy%whichComputing))**2 / galaxy%noise**2 )
		
! Derivative of the log-likelihood, obtained applying the chain rule
		do i = 1, sdim			
			gradient(i) = gradient(i) - sum( (galaxy%synth(galaxy%mask)-galaxy%spec(galaxy%mask,galaxy%whichComputing)) / galaxy%noise**2 &
				* galaxy%synthGrad(galaxy%mask,i) )
		enddo
				
! Multiply by the derivative of the sigmoid
		if (galaxy%fixVLOS == 0) then			
			gradient(1:library%nSpec) = gradient(1:library%nSpec) * diffSigmoid(trial(1:library%nSpec), priorWeight%lower, priorWeight%upper)
			gradient(library%nSpec+1:2*library%nSpec) = gradient(library%nSpec+1:2*library%nSpec) * diffSigmoid(trial(library%nSpec+1:2*library%nSpec), priorVelocity%lower, priorVelocity%upper)
			gradient(2*library%nSpec+1:3*library%nSpec) = gradient(2*library%nSpec+1:3*library%nSpec) * diffSigmoid(trial(2*library%nSpec+1:3*library%nSpec), priorDispersion%lower, priorDispersion%upper)
			
			if (maxLikelihood) then
				do i = 1, library%nSpec
					galaxy%synthGrad(:,i) = galaxy%synthGrad(:,i) * diffSigmoid(trial(i), priorWeight%lower, priorWeight%upper)
					galaxy%synthGrad(:,i+library%nSpec) = galaxy%synthGrad(:,i+library%nSpec) * diffSigmoid(trial(i+library%nSpec), priorVelocity%lower, priorVelocity%upper)
					galaxy%synthGrad(:,i+2*library%nSpec) = galaxy%synthGrad(:,i+2*library%nSpec) * diffSigmoid(trial(i+2*library%nSpec), priorDispersion%lower, priorDispersion%upper)
				enddo
			endif
		else
			gradient(1:library%nSpec) = gradient(1:library%nSpec) * diffSigmoid(trial(1:library%nSpec), priorWeight%lower, priorWeight%upper)
			gradient(library%nSpec+1) = gradient(library%nSpec+1) * diffSigmoid(trial(library%nSpec+1), priorVelocity%lower, priorVelocity%upper)
			gradient(library%nSpec+2) = gradient(library%nSpec+2) * diffSigmoid(trial(library%nSpec+2), priorDispersion%lower, priorDispersion%upper)
			
			if (maxLikelihood) then
				do i = 1, library%nSpec
					galaxy%synthGrad(:,i) = galaxy%synthGrad(:,i) * diffSigmoid(trial(i), priorWeight%lower, priorWeight%upper)
				enddo
				galaxy%synthGrad(:,library%nSpec+1) = galaxy%synthGrad(:,library%nSpec+1) * diffSigmoid(trial(library%nSpec+1), priorVelocity%lower, priorVelocity%upper)
				galaxy%synthGrad(:,library%nSpec+2) = galaxy%synthGrad(:,library%nSpec+2) * diffSigmoid(trial(library%nSpec+2), priorDispersion%lower, priorDispersion%upper)
			endif
		endif
					
! Total posterior
		slhood = logLikelihood + logPrior
		
! Save the parameters that give the best likelihood
		if (slhood > maxslhood) then
			map_pars = trial
			maxslhood = slhood
 		endif
 		 		 		
		return
			
	end subroutine slikelihoodNoNormalization

end module like_m