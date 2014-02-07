module geodesiclm_m
use like_m, only : slikelihoodNoNormalization
use params, only : galaxy, covarianceLM

	real(kind=8), allocatable :: xvar(:), fvec(:), fjac(:,:), dtd(:,:)
	integer :: m, n, mode, niters, nfev, njev, naev, maxiters, maxfev, maxjev, converged, print_level, info
	integer :: maxaev, print_unit, imethod, iaccel, ibroyden, ibold
	logical :: analytic_jac, analytic_avv, center_diff
	real(kind=8) :: eps, h1, h2, maxlam, artol, cgoal, gtol, xtol, xrtol, ftol, frtol, initial_factor
	real(kind=8) :: factoraccept, factorreject, avmax
		
contains

!-----------------------------------------------------------------------
! Solve the optimization problem using a geodesic Levenberg-Marquardt algorithm
!-----------------------------------------------------------------------
	subroutine maximumLikelihoodGeodesicLM(sdim, x)
	integer :: sdim, i
	real(kind=8) :: x(sdim)

		n = sdim
		m = galaxy%nMask
		
		if (allocated(xvar)) deallocate(xvar)
		if (allocated(fvec)) deallocate(fvec)
		if (allocated(fjac)) deallocate(fjac)
		if (allocated(dtd)) deallocate(dtd)
		if (allocated(covarianceLM)) deallocate(covarianceLM)
		
		allocate(xvar(n), fvec(M), fjac(M,N), dtd(n,n), covarianceLM(n,n))
		
		dtd = 0.d0
		do i = 1, n
			dtd(i,i) = 1.d0
		enddo
		
				
		analytic_jac = .TRUE.
		analytic_Avv = .FALSE.
		center_diff = .FALSE.
		
		xvar = x
		
		eps = 1.d-14
		h1 = 0.1d0
		h2 = 0.1d0
		mode = 1
		
		maxiters = 100
		maxfev = 0
		maxjev = 0
		maxlam = -1
		maxaev = 0
		artol = 1.d-15
		cgoal = 1.d-18
		gtol = 1.d-18
		xtol = 1.d-18
		xrtol = 1.d-18
		frtol = 1.d-18
		ftol = 1.d-15
		print_level = 3
		print_unit = 6
		imethod = 0
		initial_factor = 0.1d0
		factoraccept = 10
		factorreject = 10
		avmax = 2
		info = 0
		ibroyden = 0
		
		iaccel = 1
	 	ibold = 1     ! Type of acceptance
		
		call geolevmar(func, jacobian, Avv, xvar, fvec, fjac, n, m, callback, info,&
					analytic_jac, analytic_Avv, center_diff, eps, h1, h2,&
					dtd, mode, niters, nfev, njev, naev,&
					maxiters, maxfev, maxjev, maxaev, maxlam,&
					artol, Cgoal, gtol, xtol, xrtol, ftol, frtol,&
					converged, print_level, print_unit,&
					imethod, iaccel, ibold, ibroyden,&
					initial_factor, factoraccept, factorreject, avmax)
					
		x = xvar
		
		covarianceLM = matmul( transpose(fjac), fjac)
				
	end subroutine maximumLikelihoodGeodesicLM

!-----------------------------------------------------------------------
! Forward problem
!-----------------------------------------------------------------------
	subroutine func(m, n, x, fvec)
	integer :: m, n
	real(kind=8) :: x(n), fvec(m), f, g(n)
		
 		call slikelihoodNoNormalization(x,f,g,.TRUE.)
		
		fvec = (galaxy%synth(galaxy%mask)-galaxy%spec(galaxy%mask)) / galaxy%noise
									
	end subroutine func

!-----------------------------------------------------------------------
! Jacobian
!-----------------------------------------------------------------------
	subroutine jacobian(m, n, x, fjac)
	integer :: m, n
	real(kind=8) :: x(n), fjac(m,n), f, g(n)
				
 		call slikelihoodNoNormalization(x,f,g,.TRUE.)
		
		fjac = galaxy%synthGrad(galaxy%mask,:) / galaxy%noise		
							
	end subroutine jacobian

!-----------------------------------------------------------------------
! Directional derivative
!-----------------------------------------------------------------------
	subroutine Avv(m, n, x, v, acc)
	integer :: m, n
	real(kind=8) :: x(n), v(n), acc(m)
		
		
	end subroutine Avv

!-----------------------------------------------------------------------
! Callback
!-----------------------------------------------------------------------
	subroutine callback(m,n,x,fvec,info)
	integer :: m, n, info
	real(kind=8) :: x(n), fvec(m)
		return
	end subroutine callback

end module geodesiclm_m