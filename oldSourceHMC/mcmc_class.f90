!*****************************************************************
! This module defines a class to do MCMC sampling using several algorithms:
!   - Adaptive Metropolis
!   - Metropolis-within-Gibbs
!
! Everything is general except for the function eval_posterior, which evaluates
! the log-posterior distribution of your problem for a given set of parameters
! Since you might need additional data to evaluate the posterior (priors, data, etc.)
! this function is usually overriden by a child class defined outwards
!*****************************************************************
module mcmc_class
use maths, only : randomu, randomn, init_random_seed, mrandomn, computeCholesky
use like_m, only : computeLikelihoodAdaptive, initValuesAdaptive
	implicit none
! 	private
	real(kind=8) :: pi = 3.1415926535897931d0 ! Class-wide private constant

	type, public :: mcmc
		integer :: npar, step, step_accepted, logp, logptrial, steps_accepted_in_batch, steps_in_batch, algorithmID
		integer :: actual_chain_length, nCov
		real(kind=8), pointer :: pars(:), chain(:,:), mean(:), covariance(:,:), identity(:,:), chunk(:,:), chunk2(:,:)
		real(kind=8), pointer :: trial(:), ls(:), lsaccepted(:), cholesky(:,:)
		real(kind=8) :: alpha
		character(len=40) :: algorithm
	contains
		procedure :: initChain
		procedure :: stepChain
		procedure :: stepMGChain
		procedure :: thinChain
		procedure :: burninChain
		procedure :: updateStats
		procedure :: selectTrial
		procedure :: evalPosterior => logPosterior
		procedure :: initialValues => initialValues
	end type mcmc

	contains

!------------------------------------------------------------------
! Markov chain initialization
! npar: number of parameters of the chain
! algorithm: 'ADAPTIVE' / 'METROPOLIS_GIBBS'
!------------------------------------------------------------------
	subroutine initChain(this, npar, algorithm)
	class(mcmc), intent(inout) :: this
	integer, intent(in) :: npar
	character(len=*), intent(in) :: algorithm
	integer :: i, outBounds

! Set number of parameters
		this%npar = npar

! Set sampling algorithm
		this%algorithm = algorithm

! Allocate memory for all variables
		allocate(this%pars(npar))
		allocate(this%chain(10000,npar))
		allocate(this%mean(npar))
		allocate(this%covariance(npar,npar))
		allocate(this%cholesky(npar,npar))
		allocate(this%trial(npar))
		allocate(this%identity(npar,npar))
		
		this%actual_chain_length = 10000

! Fill up the identity matrix
		this%identity = 0.d0
		do i = 1, npar
			this%identity(i,i) = 1.d0
		enddo

! Initialization of a few variables
		this%step = 0
		this%step_accepted = 0
		this%logp = 0.d0
		this%logptrial = 0		

! Adaptive Montecarlo
! algorithmID is set to avoid the comparison with the strings everytime we do a proposal
		if (index(this%algorithm,'ADAPTIVE_METROPOLIS')) then
			this%algorithmID = 1
			this%alpha = 2.4d0 / sqrt(1.d0*npar)
		endif

! If doing Metropolis-within-Gibbs sampling, allocate a few more variables
		if (index(this%algorithm,'METROPOLIS_GIBBS')) then
			allocate(this%ls(npar))
			allocate(this%lsaccepted(npar))
			this%steps_in_batch = 0
			this%algorithmID = 2
		endif
		
! Initialize random seed
		call init_random_seed()
		
! First random sample inside the range of parameters and evaluate posterior
		call this%initialValues()		
		this%trial = this%pars		
		this%logp = this%evalPosterior(outBounds)

! Piece of chain to use to compute the covariance matrix
		this%nCov = 1000
		allocate(this%chunk(this%nCov,npar))
		allocate(this%chunk2(this%nCov,npar))
		
	end subroutine initChain

!------------------------------------------------------------------
! Return the chain
!------------------------------------------------------------------
	function getChain(this) result(out)
	class(mcmc), intent(in) :: this
	real(kind=8) :: out(this%step,this%npar)
	
		out = this%chain(1:this%step,:)

	end function getChain

!------------------------------------------------------------------
! Markov chain step
!------------------------------------------------------------------
	subroutine stepChain(this)
	class(mcmc), intent(inout) :: this
	real(kind=8) :: logp, old_pars(this%npar), trial(this%npar), acceptance, n_batch, delta_n
	real(kind=8) :: logptrial, alpha, r, ran
	real(kind=8), allocatable :: temp(:,:)
	logical :: reset_batch
	integer :: i, outBounds

! Add one step to the chain counters
		this%step = this%step + 1
		this%steps_in_batch = this%steps_in_batch + 1
		reset_batch = .FALSE.

		logp = this%logp

! Save the original set of parameters
		this%trial = this%pars
		old_pars = this%pars

! Propose a new set of parameters
		this%trial = this%selectTrial()

! Evaluate the posterior
		logptrial = this%evalPosterior(outBounds)

! If it is inside the boundaries
		if (outBounds == 0) then

! Metropolis-Hastings step verification
			r = exp(logptrial - logp)
			alpha = min(1.d0,r)
			ran = randomu()

! If the MH step is successful, we keep the modification of the i-th parameter
			if (ran < alpha) then
				this%pars = this%trial
				logp = logptrial
				this%steps_accepted_in_batch = this%steps_accepted_in_batch + 1
			endif

! Modify the logarithm of the variances to adapt the acceptance rate
			if ( modulo(this%step+1, 1000) == 0) then
				acceptance = 1.d0 * this%steps_accepted_in_batch / this%steps_in_batch				

				if (acceptance >= 0.44d0) then
					this%alpha = this%alpha * 1.05d0
				endif
				if (acceptance < 0.23d0) then
					this%alpha = this%alpha / 1.05d0
				endif

				write(*,FMT='(A,F6.3,A,F10.4)') '       Acceptance= ', acceptance*100.d0, ' - alpha=', this%alpha

				reset_batch = .TRUE.

				call this%updateStats()
				
			endif

		endif

		this%logp = logp

! Update the chain and allocate more memory for the new chunk if necessary
		if (this%step > this%actual_chain_length) then
			allocate(temp(this%actual_chain_length,this%npar))
			temp(1:this%actual_chain_length,:) = this%chain
			deallocate(this%chain)
			allocate(this%chain(this%actual_chain_length+10000,this%npar))
			this%chain(1:this%actual_chain_length,:) = temp
			deallocate(temp)
			this%actual_chain_length = this%actual_chain_length + 10000
		endif

		this%chain(this%step, :) = this%pars

! The next iteration starts a new batch
		if (reset_batch) then
			this%steps_in_batch = 0
			this%steps_accepted_in_batch = 0
		endif
		
	end subroutine stepChain

!------------------------------------------------------------------
! Markov chain Metropolis-within-Gibbs step
!------------------------------------------------------------------
	subroutine stepMGChain(this)
	class(mcmc), intent(inout) :: this
	real(kind=8) :: logp, old_pars(this%npar), trial(this%npar), acceptance, n_batch, delta_n
	real(kind=8) :: logptrial, alpha, r, ran
	real(kind=8), allocatable :: temp(:,:)
	logical :: reset_batch
	integer :: i, outBounds
		
! Add one step to the chain counters
		this%step = this%step + 1
		this%steps_in_batch = this%steps_in_batch + 1
		reset_batch = .FALSE.
		
		logp = this%logp

! Save the original set of parameters	
		this%trial = this%pars
		old_pars = this%pars

! Propose a new set of parameters
		trial = this%selectTrial()

! Do the Gibbs updating
		do i = 1, this%npar

! Set the i-th parameter to its new value and do a Metropolis step
			this%trial(i) = trial(i)

! Evaluate the posterior
			logptrial = this%evalPosterior(outBounds)
			
! If it is inside the boundaries
			if (outBounds == 0) then

! Metropolis-Hastings step verification
				r = exp(logptrial - logp)
				alpha = min(1.d0,r)
				ran = randomu()

! If the MH step is successful, we keep the modification of the i-th parameter
				if (ran < alpha) then
					this%pars(i) = this%trial(i)
					this%lsaccepted(i) = this%lsaccepted(i) + 1
					logp = logptrial
				else
! If not, return the original value of the parameter
					this%trial(i) = old_pars(i)
				endif
			else
! If not, return the original value of the parameter
				this%trial(i) = old_pars(i)
			endif

! Modify the logarithm of the variances to adapt the acceptance rate
			if ( modulo(this%step+1, 50) == 0) then
				acceptance = 1.d0 * this%lsaccepted(i) / this%steps_in_batch
				n_batch = this%step / 50.d0

				delta_n = min(0.01d0, 1.d0 / sqrt(n_batch))				
				
				if (acceptance >= 0.44d0) then
					this%ls(i) = this%ls(i) + delta_n
				endif
				if (acceptance < 0.44d0) then
					this%ls(i) = this%ls(i) - delta_n
				endif
				
				reset_batch = .TRUE.
			endif

		enddo

		this%logp = logp
		
! Update the chain and allocate more memory for the new chunk if necessary
		if (this%step > this%actual_chain_length) then
			allocate(temp(this%actual_chain_length,this%npar))
			temp(1:this%actual_chain_length,:) = this%chain
			deallocate(this%chain)
			allocate(this%chain(this%actual_chain_length+10000,this%npar))	
			this%chain(1:this%actual_chain_length,:) = temp
			deallocate(temp)
			this%actual_chain_length = this%actual_chain_length + 10000
		endif
		
		this%chain(this%step, :) = this%pars
		
! The next iteration starts a new batch
		if (reset_batch) then
			this%lsaccepted = 0
			this%steps_in_batch = 0
		endif
		
	end subroutine stepMGChain

!------------------------------------------------------------------
! Markov chain thinning by a factor
!------------------------------------------------------------------
	subroutine thinChain(this, factor)
	class(mcmc), intent(inout) :: this
	integer :: factor
	real(kind=8), allocatable :: temp(:,:)
	integer :: i, loop
	
		allocate(temp(this%step / factor,this%npar))
		loop = 1
		do i = 1, this%step, factor
			temp(loop,:) = this%chain(i,:)
		enddo
		
! Now reallocate memory for the chain
		deallocate(this%chain)
		allocate(this%chain(loop,this%npar))	
		this%chain = temp(1:loop,:)
		deallocate(temp)
		
! Change the length of the chain
		this%step = loop		
		
	end subroutine thinChain

!------------------------------------------------------------------
! Markov chain burn-in
!------------------------------------------------------------------
	subroutine burninChain(this, initial_steps)
	class(mcmc), intent(inout) :: this
	integer :: initial_steps
	real(kind=8), allocatable :: temp(:,:)
	integer :: i, loop
	
		allocate(temp(this%step - initial_steps,this%npar))
		temp = this%chain(initial_steps:this%step,:)
		
! Now reallocate memory for the chain
		deallocate(this%chain)
		allocate(this%chain(this%step - initial_steps,this%npar))	
		this%chain = temp
		deallocate(temp)
		
! Change the length of the chain
		this%step = this%step - initial_steps
	end subroutine burninChain

!------------------------------------------------------------------
! Update mean and covariance matrix
!------------------------------------------------------------------
	subroutine updateStats(this)
	class(mcmc), intent(in) :: this
	integer :: i, j, indexFrom, indexTo

! Extract chunk of chain
		indexFrom = this%step - this%nCov + 1
		indexTo = this%step
		this%chunk = this%chain(indexFrom:indexTo,:)
		
! Update the mean
		do i = 1, this%npar
			this%mean(i) = sum(this%chunk(:,i)) / (indexTo-indexFrom)
		enddo

! Update the covariance matrix
		do i = 1, this%npar
			this%chunk2(:,i) = this%chunk(:,i) - this%mean(i)
		enddo

		this%covariance = matmul(transpose(this%chunk2), this%chunk2) / (1.d0*this%nCov)

		call computeCholesky(this%alpha**2 * this%covariance, this%cholesky)

	end subroutine updateStats

!------------------------------------------------------------------
! Select a new trial point following a certain proposal density
!------------------------------------------------------------------
	function selectTrial(this) result (out)
	class(mcmc), intent(in) :: this
	real(kind=8) :: out(this%npar)
	integer :: i

! Adaptive Metropolis
		if (this%algorithmID == 1) then
			out = mrandomn(this%pars, this%cholesky)
		endif
		
! Adaptive Metropolis-within-Gibbs
! Draw npar random numbers from N(0,sigma^2), with sigma appropriate for each parameter
		if (this%algorithmID == 2) then
			do i = 1, this%npar
				out(i) = this%pars(i) + exp(this%ls(i)) * randomn()
			enddo
		endif
		
		return
		
	end function selectTrial

!------------------------------------------------------------------
! Set initial values for the parameters
!------------------------------------------------------------------
	subroutine initialValues(this)
	class(mcmc), intent(inout) :: this
	real(kind=8), allocatable :: variances(:)
	integer :: i

		allocate(variances(this%npar))

		call initValuesAdaptive(this%pars, variances)

		do i = 1, this%npar
			 this%covariance(i,i) = variances(i)**2
		enddo
		
		deallocate(variances)

		call computeCholesky(this%covariance, this%cholesky)

	end subroutine initialValues

!------------------------------------------------------------------
! Evaluate the posterior
!------------------------------------------------------------------
	function logPosterior(this, outBounds) result(logP)
	class(mcmc), intent(in) :: this
	integer, intent(inout) :: outBounds
	real(kind=8) :: logP
	
		call computeLikelihoodAdaptive(this%trial,outBounds,logP)
		
	end function logPosterior

end module mcmc_class