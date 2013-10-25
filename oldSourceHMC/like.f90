module like_m
use params
use maths, only : convolution, sigmoid, diffSigmoid
implicit none
      
contains
	
!-----------------------------------------------------------------------
! Likelihood function
!-----------------------------------------------------------------------
	subroutine slikelihoodNoNormalization(trial,slhood,gradient,sigmoidTransform)
	logical :: sigmoidTransform
	real(kind=8) :: trial(:), slhood, logPrior, logLikelihood, gradient(sdim)
	real(kind=8) :: temp(sdim),dist,loclik, sigma, xi, L, maxDispersion, sumKernel, covar(sdim,sdim), sumKernel2
	integer :: i,j, loop, nPixKernel, midPix, activeWeight(library%nSpec)

! Take also into account the priors
		logPrior = 0.d0
		gradient = 0.d0

		loop = 1
				
! ********* Prior for the weights
		do i = 1, library%nSpec
			if (sigmoidTransform) then
				galaxy%trialWeight(i) = sigmoid(trial(loop)) !1.d0 / (1.d0+exp(-trial(loop)))
			else
				galaxy%trialWeight(i) = trial(loop)
			endif

			if (trial(loop) < priorWeight%lower .or. trial(loop) > priorWeight%upper) then
				slhood = -1e15
				gradient = -1e15
 				print *, 'out w', galaxy%trialWeight(i), trial(loop)
				return
			endif
			
! Gaussian prior
			if (priorWeight%typ == 1) then
				logPrior = logPrior - (galaxy%trialWeight(i) - priorWeight%mu)**2 / (2.d0*priorWeight%sigma**2)
				gradient(loop) = gradient(loop) - (galaxy%trialWeight(i) - priorWeight%mu) / priorWeight%sigma**2
			endif

! Dirac prior
			if (priorWeight%typ == 2) then
				galaxy%trialWeight(i) = priorWeight%mu
			endif
			
			loop = loop + 1
		enddo
		

!******************
! If we do not fix the LOS velocity of each component
!******************

! ********* Prior for the velocity
		if (galaxy%fixVLOS == 0) then
			do i = 1, library%nSpec
				galaxy%trialVelocity(i) = trial(loop)

				if (trial(loop) > priorVelocity%upper .or. trial(loop) < priorVelocity%lower) then
					slhood = -1e15
					gradient = -1e15
 					print *, 'out v', galaxy%trialVelocity(i)
					return
				endif

! Gaussian prior
				if (priorVelocity%typ == 1) then
					logPrior = logPrior - (galaxy%trialVelocity(i) - priorVelocity%mu)**2 / (2.d0*priorVelocity%sigma**2)
					gradient(loop) = gradient(loop) - (galaxy%trialVelocity(i) - priorVelocity%mu) / priorVelocity%sigma**2
				endif

! Dirac prior
				if (priorVelocity%typ == 2) then
					galaxy%trialVelocity(i) = priorVelocity%mu
				endif

				loop = loop + 1
			enddo

! ********* Prior for the dispersion
			do i = 1, library%nSpec
				galaxy%trialDispersion(i) = trial(loop)

				if (trial(loop) > priorDispersion%upper .or. trial(loop) < priorDispersion%lower) then
					slhood = -1e15
					gradient = -1e15
					print *, 'out d', galaxy%trialDispersion(i)
					return
				endif

	! Gaussian prior
				if (priorDispersion%typ == 1) then
					logPrior = logPrior - (galaxy%trialDispersion(i) - priorDispersion%mu)**2 / (2.d0*priorDispersion%sigma**2)
					gradient(loop) = gradient(loop) - (galaxy%trialDispersion(i) - priorDispersion%mu) / priorDispersion%sigma**2
				endif

	! Dirac prior
				if (priorDispersion%typ == 2) then
					galaxy%trialDispersion(i) = priorDispersion%mu
				endif
				
				loop = loop + 1
			enddo
		else
		
!******************
! If we fix the LOS velocity of each component, fill all the trialVelocity array with the same value
!******************
			galaxy%trialVelocity = trial(loop)

			if (trial(loop) > priorVelocity%upper .or. trial(loop) < priorVelocity%lower) then
				slhood = -1e15
				gradient = -1e15
				return
			endif

! Dirac prior
			if (priorVelocity%typ == 2) then
				galaxy%trialVelocity = priorVelocity%mu
			endif
			
! Gaussian prior
			if (priorVelocity%typ == 1) then
				logPrior = logPrior - (galaxy%trialVelocity(1) - priorVelocity%mu)**2 / (2.d0*priorVelocity%sigma**2)
				gradient(loop) = gradient(loop) - (galaxy%trialVelocity(i) - priorVelocity%mu) / priorVelocity%sigma**2
			endif
			loop = loop + 1
			
			
			galaxy%trialDispersion = trial(loop)

			if (trial(loop) > priorDispersion%upper .or. trial(loop) < priorDispersion%lower) then
				slhood = -1e15
				gradient = -1e15
				return
			endif

! Dirac prior
			if (priorDispersion%typ == 2) then
				galaxy%trialDispersion = priorDispersion%mu
			endif
			
! Gaussian prior
			if (priorDispersion%typ == 1) then
				logPrior = logPrior - (galaxy%trialDispersion(1) - priorDispersion%mu)**2 / (2.d0*priorDispersion%sigma**2)
				gradient(loop) = gradient(loop) - (galaxy%trialDispersion(i) - priorDispersion%mu) / priorDispersion%sigma**2
			endif
			loop = loop + 1
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
				
! Normalize the velocity and dispersion to the global scale
		galaxy%trialDispersion = galaxy%trialDispersion / library%velScale
		galaxy%trialVelocity = galaxy%trialVelocity / library%velScale

		maxDispersion = maxval(galaxy%trialDispersion)

! Find the number of pixels of the kernel
		nPixKernel = 2*ceiling(4*maxDispersion)+1

! Do something more efficient that allocating memory for each iteration
		allocate(galaxy%xKernel(nPixKernel))
		allocate(galaxy%yKernel(nPixKernel))		
		allocate(galaxy%yKernelDerivative(nPixKernel))
		
		allocate(galaxy%synth2(galaxy%nPix))

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
					
		do i = 1, library%nSpec
			if (sigmoidTransform) then
				galaxy%synthGrad(:,i) = galaxy%synthGrad(:,i) * diffSigmoid(trial(i)) !exp(-trial(i)) * galaxy%trialWeight(i)**2 * activeWeight(i)			
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