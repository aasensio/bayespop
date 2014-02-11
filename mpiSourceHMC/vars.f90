! Include file for example nested sampler program 'Gaussian Rings'

module params
implicit none

! Toy Model Parameters

!dimensionality
      	integer :: sdim


      
!tot no. of parameters, should be sdim in most cases but if you need to
!store some additional parameters with the actual parameters then
!you need to pass them through the likelihood routine
			integer :: nest_nPar
      
      	!root for saving posterior files
      	character*100 nest_root
! 			parameter(nest_root='chains/2-')
      	integer maxiter
      	
!=======================================================================
	
	real(kind=8), allocatable :: st(:), stepSize(:), savePars(:), covarianceLM(:,:), hessianLM(:,:)
	real(kind=8) :: scaleFactor
	integer :: seed, fbInt, maxStep, resume
	character(len=128) :: flPfx

	real(kind=8), parameter :: lightSpeed_c = 299792.458d0
	character(len=20), allocatable, dimension(:) :: fvars
	character(len=250) :: feval
	
	real(kind=8) :: sigma, maxslhood, logevidence
	real(kind=8), allocatable, dimension(:) :: map_pars
	real(kind=8), parameter :: PI = 3.14159265359d0, one_sigma = 68.268955e0, &
		two_sigma = 95.449972e0

	integer :: ndata
	character(len=128) :: file_data
	real(kind=8), allocatable :: x(:), y(:), noise(:), model(:), fvals(:)
	integer, allocatable :: pWrap(:)

	type prior_t
		integer :: typ
      real(kind=8) :: mu, sigma
      real(kind=8) :: lower, upper
      integer :: pWrap
	end type prior_t

	type markov_t
		character(len=4) :: typ
		integer :: nparam, nlength		
		real(kind=8), pointer :: chain(:,:), best_parameters(:), most_probable(:)
		real(kind=8) :: evidence, avg_lnL
		character(len=100) :: filename
		character(len=4), pointer :: param_names(:)
	end type markov_t

	type library_t
		character(len=120) :: filename, name
		integer :: nSpec, nPix, nPixAdapted
		real(kind=8) :: velScale, fwhm
		real(kind=8), pointer :: age(:), metallicity(:), imfSlope(:), mgfe(:)
		real(kind=8), pointer :: spec(:,:), specAdapted(:,:), lambda(:)
	end type library_t

	type galaxy_t
		character(len=120) :: filename
		integer :: nPix, nMask, nSpec, whichComputing, fixVLOS
		real(kind=8) :: snr, noise, velScale
		real(kind=8), pointer :: spec(:), synth(:), xKernel(:), yKernel(:), lambda(:), synthGrad(:,:), convolvedSpectrum(:), synth2(:), specRead(:,:)
		real(kind=8), pointer :: trialWeight(:), trialVelocity(:), trialDispersion(:), yKernelDerivative(:), trial(:)
		integer, pointer :: maskIndex(:), mask(:)
	end type galaxy_t

	type(markov_t) :: chain_analysis
	type(prior_t), pointer :: priors(:)
	
	character(len=11) :: variableNames(6)=(/'Velocity   ','Dispersion ','Age        ','Metallicity','IMF Slope  ','MGFE       '/)

	type(library_t) :: library
	type(galaxy_t) :: galaxy
	
	integer :: nCases, myrank
	integer, allocatable :: nGalaxies(:), whichCase(:), whichGalaxy(:)
	character(len=120), allocatable :: fileCases(:)
	character(len=4) :: myrankStr

end module params