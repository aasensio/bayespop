module like_m
use params
use maths, only : convolution, sigmoid, diffSigmoid, hunt
implicit none
      
contains

!-----------------------------------------------------------------------
! Linear interpolation
!-----------------------------------------------------------------------
	subroutine linearInterpolation(pars, modelIndex)
	real(kind=8) :: pars(2), delta, deltaJ
	integer :: modelIndex, i, j, ind, near(2), indices(2)
			
 		call hunt(library%velocity, library%nVelocities, pars(1), near(1))
 		call hunt(library%dispersion, library%nDispersions, pars(2), near(2))
 								
! Extract the values of the function that will be used
		do i = 1, 2**library%nDim
			indices = near + library%ee(:,i)			
			library%wrk(:,i) = library%spec(:,modelIndex,indices(1),indices(2))
  			do j = 1, library%nDim
  				library%wrkJ(:,j,i) = library%spec(:,modelIndex,indices(1),indices(2))
  			enddo
		enddo
		
 		library%jacobian = 0.d0
! Do the actual linear interpolation
		do i = 1, library%nDim
			ind = library%indi(i)			
			delta = -(pars(ind) - library%params(ind,near(ind))) / &
				(library%params(ind,near(ind)) - library%params(ind,near(ind)+1))
			deltaJ = -1.d0 / (library%params(ind,near(ind)) - library%params(ind,near(ind)+1))
			
			do j = 1, 2**(library%nDim-i)				
				library%wrk(:,j) = (library%wrk(:,2*j) - library%wrk(:,2*j-1)) * delta + library%wrk(:,2*j-1)
  				library%wrkJ(:,ind,j) = (library%wrkJ(:,ind,2*j) - library%wrkJ(:,ind,2*j-1)) * deltaJ
			enddo
			
		enddo
		
 		library%spectrum = library%wrk(:,1)
  		library%jacobian = library%wrkJ(:,:,1)
		
	end subroutine linearInterpolation

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
		galaxy%trialWeight = sigmoid(trial(1:library%nSpec), 0.d0, 1.d0)			
		
! 		print *, trial
! 		stop

!******************
! If we do not fix the LOS velocity of each component
!******************
		if (galaxy%fixVLOS == 0) then

! ********* Prior for the velocity
			galaxy%trialVelocity = sigmoid(trial(library%nSpec+1:2*library%nSpec), priors(1)%lower, priors(1)%upper)

! Gaussian prior
			if (priors(1)%typ == 1) then
				logPrior = logPrior - sum((galaxy%trialVelocity - priors(1)%mu)**2 / (2.d0*priors(1)%sigma**2))
				gradient(library%nSpec+1:2*library%nSpec) = gradient(library%nSpec+1:2*library%nSpec) - (galaxy%trialVelocity - priors(1)%mu) / priors(1)%sigma**2
			endif

! Dirac prior
			if (priorVelocity%typ == 2) then
				galaxy%trialVelocity = priors(1)%mu
			endif
			
! ********* Prior for the dispersion
			galaxy%trialDispersion = sigmoid(trial(2*library%nSpec+1:3*library%nSpec), priors(2)%lower, priors(2)%upper)

! Gaussian prior
			if (priors(2)%typ == 1) then
				logPrior = logPrior - sum((galaxy%trialDispersion - priors(2)%mu)**2 / (2.d0*priors(2)%sigma**2))
				gradient(2*library%nSpec+1:3*library%nSpec) = gradient(2*library%nSpec+1:3*library%nSpec) - (galaxy%trialDispersion - priors(2)%mu) / priors(2)%sigma**2
			endif

! Dirac prior
			if (priors(2)%typ == 2) then
				galaxy%trialDispersion = priors(2)%mu
			endif
			
		else

!******************
! If we fix the LOS velocity of each component, fill all the trialVelocity array with the same value
!******************			
! ********* Prior for the velocity
			galaxy%trialVelocity = sigmoid(trial(library%nSpec+1), priors(1)%lower, priors(1)%upper)

! Gaussian prior
			if (priors(1)%typ == 1) then
				logPrior = logPrior - (galaxy%trialVelocity(1) - priors(1)%mu)**2 / (2.d0*priors(1)%sigma**2)
				gradient(library%nSpec+1) = gradient(library%nSpec+1) - (galaxy%trialVelocity(1) - priors(1)%mu) / priors(1)%sigma**2
			endif

! Dirac prior
			if (priors(1)%typ == 2) then
				galaxy%trialVelocity = priors(1)%mu
			endif
			
! ********* Prior for the dispersion
			galaxy%trialDispersion = sigmoid(trial(library%nSpec+2), priors(2)%lower, priors(2)%upper)
					
! Gaussian prior
			if (priors(2)%typ == 1) then
				logPrior = logPrior - (galaxy%trialDispersion(1) - priors(2)%mu)**2 / (2.d0*priors(2)%sigma**2)
				gradient(library%nSpec+2) = gradient(library%nSpec+2) - (galaxy%trialDispersion(1) - priors(2)%mu) / priors(2)%sigma**2
			endif

! Dirac prior
			if (priors(2)%typ == 2) then
				galaxy%trialDispersion = priors(2)%mu
			endif
					
		endif									
		
! If the weight of the component is smaller than 0.001, set its weight to zero
		activeWeight = 1

! **************************
! ********** DATA LIKELIHOOD
! **************************		
						
		galaxy%synth = 0.d0		
		galaxy%synthGrad = 0.d0
		
! Compute the kernel, do the convolution, add the component and compute the derivatives
		do i = 1, library%nSpec

			if (activeWeight(i) == 1) then
								
  				call linearInterpolation((/galaxy%trialVelocity(i),galaxy%trialDispersion(i)/), i)
				
				galaxy%synth = galaxy%synth + activeWeight(i) * galaxy%trialWeight(i) * library%spectrum
				
				galaxy%synthGrad(:,i) = library%spectrum

! Every component has different kinematics
				if (galaxy%fixVLOS == 0) then

					galaxy%synthGrad(:,i+library%nSpec) = activeWeight(i) * galaxy%trialWeight(i) * library%jacobian(:,1)
					galaxy%synthGrad(:,i+2*library%nSpec) = activeWeight(i) * galaxy%trialWeight(i) * library%jacobian(:,2)
																
				else
							
					galaxy%synthGrad(:,1+library%nSpec) = galaxy%synthGrad(:,1+library%nSpec) + activeWeight(i) * galaxy%trialWeight(i) * library%jacobian(:,1)
					galaxy%synthGrad(:,2+library%nSpec) = galaxy%synthGrad(:,2+library%nSpec) + activeWeight(i) * galaxy%trialWeight(i) * library%jacobian(:,2)
					
				endif
												
			endif
			
		enddo
									
! Compute the log-likelihood
		logLikelihood = -0.5d0 * sum( (galaxy%synth(galaxy%mask)-galaxy%spec(galaxy%mask,galaxy%whichComputing))**2 / galaxy%noise**2 )
		
! Derivative of the log-likelihood, obtained applying the chain rule
		do i = 1, sdim			
			gradient(i) = gradient(i) - sum( (galaxy%synth(galaxy%mask)-galaxy%spec(galaxy%mask,galaxy%whichComputing)) / galaxy%noise**2 &
				* galaxy%synthGrad(galaxy%mask,i) )
		enddo
				
! Multiply by the derivative of the sigmoid
		if (galaxy%fixVLOS == 0) then			
			gradient(1:library%nSpec) = gradient(1:library%nSpec) * diffSigmoid(trial(1:library%nSpec), 0.d0, 1.d0)
			gradient(library%nSpec+1:2*library%nSpec) = gradient(library%nSpec+1:2*library%nSpec) * diffSigmoid(trial(library%nSpec+1:2*library%nSpec), priors(1)%lower, priors(1)%upper)
			gradient(2*library%nSpec+1:3*library%nSpec) = gradient(2*library%nSpec+1:3*library%nSpec) * diffSigmoid(trial(2*library%nSpec+1:3*library%nSpec), priors(2)%lower, priors(2)%upper)
			
			if (maxLikelihood) then
				do i = 1, library%nSpec
					galaxy%synthGrad(:,i) = galaxy%synthGrad(:,i) * diffSigmoid(trial(i), 0.d0, 1.d0)
					galaxy%synthGrad(:,i+library%nSpec) = galaxy%synthGrad(:,i+library%nSpec) * diffSigmoid(trial(i+library%nSpec), priors(1)%lower, priors(1)%upper)
					galaxy%synthGrad(:,i+2*library%nSpec) = galaxy%synthGrad(:,i+2*library%nSpec) * diffSigmoid(trial(i+2*library%nSpec), priors(2)%lower, priors(2)%upper)
				enddo
			endif
		else
			gradient(1:library%nSpec) = gradient(1:library%nSpec) * diffSigmoid(trial(1:library%nSpec), 0.d0, 1.d0)
			gradient(library%nSpec+1) = gradient(library%nSpec+1) * diffSigmoid(trial(library%nSpec+1), priors(1)%lower, priors(1)%upper)
			gradient(library%nSpec+2) = gradient(library%nSpec+2) * diffSigmoid(trial(library%nSpec+2), priors(2)%lower, priors(2)%upper)
			
			if (maxLikelihood) then
				do i = 1, library%nSpec
					galaxy%synthGrad(:,i) = galaxy%synthGrad(:,i) * diffSigmoid(trial(i), 0.d0, 1.d0)
				enddo
				galaxy%synthGrad(:,library%nSpec+1) = galaxy%synthGrad(:,library%nSpec+1) * diffSigmoid(trial(library%nSpec+1), priors(1)%lower, priors(1)%upper)
				galaxy%synthGrad(:,library%nSpec+2) = galaxy%synthGrad(:,library%nSpec+2) * diffSigmoid(trial(library%nSpec+2), priors(2)%lower, priors(2)%upper)
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