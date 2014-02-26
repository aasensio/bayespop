module maths
use params, only : indexArrayConvolutionCentral
implicit none

contains

!-----------------------------------------------------------------------
! Return a sigmoid function
!-----------------------------------------------------------------------
	elemental function sigmoid(x, lower, upper) result(y)
	real(kind=8), intent(in) :: x, lower, upper
	real(kind=8) :: y
		y = lower + (upper-lower) / (1.d0 + exp(-x))
	end function sigmoid
	
!-----------------------------------------------------------------------
! Return the inverse of the sigmoid function
!-----------------------------------------------------------------------
	elemental function invSigmoid(x, lower, upper) result(y)
	real(kind=8), intent(in) :: x, lower, upper
	real(kind=8) :: y
		y = log( (lower - x) / (x - upper) )
	end function invSigmoid
	
!-----------------------------------------------------------------------
! Return the derivative of the sigmoid function
!-----------------------------------------------------------------------
	elemental function diffSigmoid(x, lower, upper) result(y)
	real(kind=8), intent(in) :: x, lower, upper
	real(kind=8) :: y
		y = (upper-lower) * exp(-x) / (1.d0+exp(-x))**2
	end function diffSigmoid
	
!-----------------------------------------------------------------------
! Return the derivative of the sigmoid function
!-----------------------------------------------------------------------
	elemental function lnJacSigmoid(x, lower, upper) result(y)
	real(kind=8), intent(in) :: x, lower, upper
	real(kind=8) :: y
		y = log(upper-lower) - x - 2.d0 * log(1.d0 + exp(-x))
	end function lnJacSigmoid

!-----------------------------------------------------------------------
! Return the derivative of the sigmoid function
!-----------------------------------------------------------------------
	elemental function difflnJacSigmoid(x, lower, upper) result(y)
	real(kind=8), intent(in) :: x, lower, upper
	real(kind=8) :: y
		y = -1.d0 + 2.d0 * exp(-x) / (1.0+exp(-x))
	end function difflnJacSigmoid

!-----------------------------------------------------------------------
! Convolution of a spectrum with a kernel. Use edge wrapping
!-----------------------------------------------------------------------
	function convolution(spec, kernel)
	real(kind=8) :: spec(:), kernel(:), convolution(size(spec))
	integer :: i, j, nSpec, nKernel, k
		nSpec = size(spec)
		nKernel = size(kernel)

		convolution = 0.d0

		do i = 1, nSpec

			do j = 1, nKernel
				k = mod(i+j-nKernel/2-1, nSpec)
				convolution(i) = convolution(i) + spec(k) * kernel(j)
			enddo
			
		enddo

	end function convolution
	
!-----------------------------------------------------------------------
! Convolution of a spectrum with a kernel. Use edge wrapping
!-----------------------------------------------------------------------
	function convolutionVectorized(spec, kernel)
	real(kind=8) :: spec(:), kernel(:), convolutionVectorized(size(spec))
	integer :: i, j, nSpec, nKernel, k
		nSpec = size(spec)
		nKernel = size(kernel)

		convolutionVectorized = 0.d0

		do i = 1, nSpec									
			convolutionVectorized(i) = convolutionVectorized(i) + sum(spec(mod(i+indexArrayConvolutionCentral, nSpec)) * kernel)
		enddo

	end function convolutionVectorized


!-------------------------------------------------------------
! Initialize the random number generator
!-------------------------------------------------------------
	subroutine init_random_seed()
	integer :: i, n, clock
   integer, dimension(:), allocatable :: seed
          
   	call random_seed(size = n)
      allocate(seed(n))
          
      call system_clock(count=clock)
          
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)
          
      deallocate(seed)
   end subroutine init_random_seed
          
!-------------------------------------------------------------
! Generates a random number following an uniform distribution in the interval [0,1]
! Call it with idum<0 for initialization
!-------------------------------------------------------------
	function randomu()	
	real(kind=8) :: randomu
		
		call random_number(randomu)
							
	end function randomu
	
!-------------------------------------------------------------
! Generates a random number following an normal distribution with zero mean
! and unit variance
!-------------------------------------------------------------
	function randomn()

	real(kind=8) :: randomn

	real(kind=8) :: u, sum
	real(kind=8), save :: v, sln
	logical, save :: second = .false.
	real(kind=8), parameter :: one = 1.0, vsmall = tiny( one )

! if second, use the second random number generated on last call
	if (second) then

		second = .false.
  		randomn = v*sln
	else
! first call; generate a pair of random normals

  		second = .true.
  		do
    		call random_number( u )
    		call random_number( v )
    		u = scale( u, 1 ) - one
    		v = scale( v, 1 ) - one
    		sum = u*u + v*v + vsmall         ! vsmall added to prevent log(zero) / zero
    		if(sum < one) exit
  		end do
  		sln = sqrt(- scale( log(sum), 1 ) / sum)
  		randomn = u*sln
	end if

	return
	end function randomn

!-------------------------------------------------------------
! Generates a random number following a Poisson distribution
!-------------------------------------------------------------
	function random_poisson(lambda)
	real(kind=8) :: lambda, L, u, p
	integer :: k, random_poisson

		L = -lambda
		k = 0
		p = 0.d0
		do while (p >= L)
			k = k + 1
			u = randomu()
			p = p + log(u)
		enddo
		random_poisson = k-1
		
	end function random_poisson

!-------------------------------------------------------------
! Generates a random number following an exponential distribution
!-------------------------------------------------------------
	function random_exp(lambda)
	real(kind=8) :: lambda, random_exp

		random_exp = -lambda * log(1.d0 - randomu())
		
	end function random_exp

!-------------------------------------------------------------
! Return the mean and standard deviation of a vector
!-------------------------------------------------------------
	function mean_stddev(x)
	real(kind=8) :: x(:), mean_stddev(2)

		mean_stddev(1) = sum(x) / size(x)
		mean_stddev(2) = sqrt( sum((x-mean_stddev(1))**2) / size(x) )
		
	end function mean_stddev

!-------------------------------------------------------------
! Generates a random number following a Maxwellian distribution
!-------------------------------------------------------------
	function random_maxw(lambda)
	real(kind=8) :: lambda, random_maxw

		random_maxw = sqrt( (lambda*randomn())**2 + (lambda*randomn())**2 + (lambda*randomn())**2 )
		
	end function random_maxw

!-------------------------------------------------------------
! Carry out the Cholesky decomposition of a symmetric matrix
!-------------------------------------------------------------
	subroutine cholesky(a,n,p)
	integer :: n
	real(kind=8) :: a(n,n), p(n)
	integer :: i, j, k
	real(kind=8) :: sum

		do i = 1, n
			do j = i, n
				sum = a(i,j)
				do k = i-1, 1, -1
					sum = sum - a(i,k)*a(j,k)
				enddo
				if (i == j) then
					if (sum == 0.d0) then
						print *, 'Cholesky decomposition failed...'
					endif
					p(i) = dsqrt(sum)
				else
					a(j,i) = sum / p(i)
				endif
			enddo
		enddo

	end subroutine cholesky

!-------------------------------------------------------------
! Compute the Cholesky decomposition of a covariance matrix
!-------------------------------------------------------------
	subroutine computeCholesky(covar, chol)
	integer :: idum, n, i, j
	real(kind=8) :: covar(:,:), chol(:,:)
	real(kind=8), allocatable :: p(:)

		n = size(chol,dim=1)

		allocate(p(n))

		chol = covar

		do i = 1, n			
			chol(i,i) = chol(i,i) + 1.d-7    ! Some regularization
		enddo

		call cholesky(chol,n,p)

		do j = 1, n
			do i = j, n
				chol(j,i) = 0.d0
			enddo
			chol(j,j) = p(j)
		enddo

		deallocate(p)

	end subroutine computeCholesky

!-------------------------------------------------------------
! Generates a multivariate normal random number with a given
! mean and covariance matrix
!-------------------------------------------------------------
	function mrandomn(rmean,cholesky)
	integer :: idum, n, i, j
	real(kind=8) :: rmean(:), cholesky(:,:), mrandomn(size(rmean))
	real(kind=8) :: eps(size(rmean))

		n = size(rmean)

		do i = 1, n
			eps(i) = randomn()
		enddo

		mrandomn = matmul(cholesky,eps) + rmean


	end function mrandomn

! ---------------------------------------------------------
!	Given an array xx(1:n), and given a value x, returns a value jlo such that x is between
!	xx(jlo) and xx(jlo+1). xx(1:n) must be monotonic, either increasing or decreasing.
!	jlo=0 or jlo=n is returned to indicate that x is out of range. jlo on input is taken as
!	the initial guess for jlo on output.
! ---------------------------------------------------------
	subroutine hunt(xx,n,x,jlo)
	integer :: jlo,n
	real(kind=8) :: x,xx(n)
	integer :: inc,jhi,jm
	logical :: ascnd

		ascnd=xx(n).ge.xx(1)
		if (jlo.le.0.or.jlo.gt.n) then
			jlo=0
			jhi=n+1
			goto 3
		endif
		inc=1
		if (x.ge.xx(jlo).eqv.ascnd) then
1     	jhi=jlo+inc
			if (jhi.gt.n) then
				jhi=n+1
			else if (x.ge.xx(jhi).eqv.ascnd) then
				jlo=jhi
				inc=inc+inc
				goto 1
			endif
		else
			jhi=jlo
2     	jlo=jhi-inc
			if (jlo.lt.1) then
				jlo=0
			else if (x.lt.xx(jlo).eqv.ascnd) then
				jhi=jlo
				inc=inc+inc
				goto 2
			endif
		endif
3 		if (jhi-jlo.eq.1) then
			if(x.eq.xx(n)) jlo=n-1
			if(x.eq.xx(1)) jlo=1
			return
		endif
		jm = (jhi+jlo)/2
		if (x.ge.xx(jm).eqv.ascnd) then
			jlo=jm
		else
			jhi=jm
		endif
		goto 3
	end subroutine hunt

!----------------------------------------------------------------
! This function integrates a tabulated function
!----------------------------------------------------------------
	function int_tabulated(x, f)
	real(kind=8) :: x(:), f(:), int_tabulated, res, error_res
	integer :: n
		n = size(x)
		call cubint (f, x, n, 1, n, res, error_res)
		int_tabulated = res
	end function int_tabulated

!-------------------------------------------------------------
! CUBINT approximates an integral using cubic interpolation of data.
!  Parameters:
!
!    Input, real FTAB(NTAB), contains the tabulated function
!    values, FTAB(I) = F(XTAB(I)).
!
!    Input, real XTAB(NTAB), contains the points at which the
!    function was tabulated.  XTAB should contain distinct
!    values, given in ascending order.
!
!    Input, integer NTAB, the number of tabulated points.
!    NTAB must be at least 4.
!
!    Input, integer IA, the entry of XTAB at which integration
!    is to begin.  IA must be no less than 1 and no greater
!    than NTAB.
!
!    Input, integer IB, the entry of XTAB at which integration
!    is to end.  IB must be no less than 1 and no greater than
!    NTAB.
!
!    Output, real RESULT, the approximate value of the
!    integral from XTAB(IA) to XTAB(IB) of the function.
!
!    Output, real ERROR, an estimate of the error in
!    integration.
!-------------------------------------------------------------
	subroutine cubint ( ftab, xtab, ntab, ia, ib, result, error )

  	integer ntab
!
  	real(kind=8) :: c, d1, d2, d3, error, ftab(ntab), h1, h2, h3, h4
  	integer :: i, ia, ib, ind, it, j, k
  	real(kind=8) r1, r2, r3, r4, result, s, term, xtab(ntab)
!
  	result = 0.0E+00
  	error = 0.0E+00

  	if ( ia == ib ) then
    	return
  	end if

  	if ( ntab < 4 ) then
    	write ( *, '(a)' ) ' '
    	write ( *, '(a)' ) 'CUBINT - Fatal error!'
    	write ( *, '(a,i6)' ) '  NTAB must be at least 4, but input NTAB = ',ntab
    	stop
  	endif

  	if ( ia < 1 ) then
    	write ( *, '(a)' ) ' '
    	write ( *, '(a)' ) 'CUBINT - Fatal error!'
    	write ( *, '(a,i6)' ) '  IA must be at least 1, but input IA = ',ia
    	stop
  	endif

  	if ( ia > ntab ) then
   	write ( *, '(a)' ) ' '
    	write ( *, '(a)' ) 'CUBINT - Fatal error!'
    	write ( *, '(a,i6)' ) '  IA must be <= NTAB, but input IA=',ia
    	stop
  	endif

  	if ( ib < 1 ) then
    	write ( *, '(a)' ) ' '
    	write ( *, '(a)' ) 'CUBINT - Fatal error!'
    	write ( *, '(a,i6)' ) '  IB must be at least 1, but input IB = ',ib
    	stop
  	endif

  	if ( ib > ntab ) then
    	write ( *, '(a)' ) ' '
    	write ( *, '(a)' ) 'CUBINT - Fatal error!'
    	write ( *, '(a,i6)' ) '  IB must be <= NTAB, but input IB=',ib
    	stop
  	endif
!
!  Temporarily switch IA and IB, and store minus sign in IND
!  so that, while integration is carried out from low X's
!  to high ones, the sense of the integral is preserved.
!
  	if ( ia > ib ) then
    	ind = -1
    	it = ib
    	ib = ia
    	ia = it
  	else
    	ind = 1
  	endif

  	s = 0.0E+00
  	c = 0.0E+00
  	r4 = 0.0E+00
  	j = ntab-2
  	if ( ia < ntab-1 .or. ntab == 4 ) then
    	j=max(3,ia)
  	endif

  	k = 4
  	if ( ib > 2 .or. ntab == 4 ) then
    	k=min(ntab,ib+2)-1
  	endif

  	do i = j, k

    	if ( i <= j ) then

      	h2 = xtab(j-1)-xtab(j-2)
      	d3 = (ftab(j-1)-ftab(j-2)) / h2
      	h3 = xtab(j)-xtab(j-1)
      	d1 = (ftab(j)-ftab(j-1)) / h3
      	h1 = h2+h3
      	d2 = (d1-d3)/h1
      	h4 = xtab(j+1)-xtab(j)
      	r1 = (ftab(j+1)-ftab(j)) / h4
      	r2 = (r1-d1) / (h4+h3)
      	h1 = h1+h4
      	r3 = (r2-d2) / h1

      	if ( ia <= 1 ) then
        		result = h2 * (ftab(1)+h2*(0.5*d3-h2*(d2/6.0-(h2+h3+h3)*r3/12.)))
        		s = -h2**3 * (h2*(3.0*h2+5.0*h4)+10.0*h3*h1)/60.0
      	endif

    	else

	      h4 = xtab(i+1)-xtab(i)
      	r1 = (ftab(i+1)-ftab(i))/h4
      	r4 = h4+h3
      	r2 = (r1-d1)/r4
      	r4 = r4+h2
      	r3 = (r2-d2)/r4
      	r4 = (r3-d3)/(r4+h1)

    	endif

    	if ( i > ia .and. i <= ib ) then

      	term = h3*((ftab(i)+ftab(i-1))*0.5-h3*h3*(d2+r2+(h2-h4)*r3) / 12.0 )
      	result = result+term
      	c = h3**3*(2.0E+00 *h3*h3+5.*(h3*(h4+h2) + 2.0 * h2 * h4 ) ) / 120.0E+00
      	error = error+(c+s)*r4

      	if ( i /= j ) then
        		s = c
      	else
        		s = s+c+c
      	endif

    	else

	      error = error+r4*s

    	endif

    	if ( i >= k ) then

      	if ( ib >= ntab ) then
	        	term = h4*(ftab(ntab) - h4*(0.5*r1+h4*(r2/6.0 +(h3+h3+h4)*r3/12.)))
        		result = result + term
        		error = error - h4**3 * r4 * &
          		( h4 * ( 3.0 * h4 + 5.0 * h2 ) &
          		+ 10.0 * h3 * ( h2 + h3 + h4 ) ) / 60.0E+00
      	endif

      	if ( ib >= ntab-1 ) error=error+s*r4
    	else
	      h1 = h2
      	h2 = h3
      	h3 = h4
      	d1 = r1
      	d2 = r2
      	d3 = r3
    	endif

  	enddo
!
!  Restore original values of IA and IB, reverse signs
!  of RESULT and ERROR, to account for integration
!  that proceeded from high X to low X.
!
  	if ( ind /= 1 ) then
    	it = ib
    	ib = ia
    	ia = it
    	result = -result
    	error = -error
  	endif

  	return
	end subroutine

	subroutine qsortd(x,ind,n)

	implicit none
	integer, parameter  :: dp = 4

	real (kind=8), intent(in)  :: x(:)
	integer, intent(out)   :: ind(:)
	integer, intent(in)    :: n

	!***************************************************************************

	!                                                         robert renka
	!                                                 oak ridge natl. lab.

	!   this subroutine uses an order n*log(n) quick sort to sort a real (dp)
	! array x into increasing order.  the algorithm is as follows.  ind is
	! initialized to the ordered sequence of indices 1,...,n, and all interchanges
	! are applied to ind.  x is divided into two portions by picking a central
	! element t.  the first and last elements are compared with t, and
	! interchanges are applied as necessary so that the three values are in
	! ascending order.  interchanges are then applied so that all elements
	! greater than t are in the upper portion of the array and all elements
	! less than t are in the lower portion.  the upper and lower indices of one
	! of the portions are saved in local arrays, and the process is repeated
	! iteratively on the other portion.  when a portion is completely sorted,
	! the process begins again by retrieving the indices bounding another
	! unsorted portion.

	! input parameters -   n - length of the array x.

	!                      x - vector of length n to be sorted.

	!                    ind - vector of length >= n.

	! n and x are not altered by this routine.

	! output parameter - ind - sequence of indices 1,...,n permuted in the same
	!                          fashion as x would be.  thus, the ordering on
	!                          x is defined by y(i) = x(ind(i)).

	!*********************************************************************

	! note -- iu and il must be dimensioned >= log(n) where log has base 2.

	!*********************************************************************

	integer   :: iu(21), il(21)
	integer   :: m, i, j, k, l, ij, it, itt, indx
	real(kind=8)      :: r
	real(kind=8) :: t

	! local parameters -

	! iu,il =  temporary storage for the upper and lower
	!            indices of portions of the array x
	! m =      index for iu and il
	! i,j =    lower and upper indices of a portion of x
	! k,l =    indices in the range i,...,j
	! ij =     randomly chosen index between i and j
	! it,itt = temporary storage for interchanges in ind
	! indx =   temporary index for x
	! r =      pseudo random number for generating ij
	! t =      central element of x

	if (n <= 0) return

	! initialize ind, m, i, j, and r

	do  i = 1, n
	ind(i) = i
	end do

	m = 1
	i = 1
	j = n
	r = .375

	! top of loop

	20 if (i >= j) go to 70
	if (r <= .5898437) then
	r = r + .0390625
	else
	r = r - .21875
	end if

	! initialize k

	30 k = i

	! select a central element of x and save it in t

	ij = i + r*(j-i)
	it = ind(ij)
	t = x(it)

	! if the first element of the array is greater than t,
	!   interchange it with t

	indx = ind(i)
	if (x(indx) > t) then
	ind(ij) = indx
	ind(i) = it
	it = indx
	t = x(it)
	end if

	! initialize l

	l = j

	! if the last element of the array is less than t,
	!   interchange it with t

	indx = ind(j)
	if (x(indx) >= t) go to 50
	ind(ij) = indx
	ind(j) = it
	it = indx
	t = x(it)

	! if the first element of the array is greater than t,
	!   interchange it with t

	indx = ind(i)
	if (x(indx) <= t) go to 50
	ind(ij) = indx
	ind(i) = it
	it = indx
	t = x(it)
	go to 50

	! interchange elements k and l

	40 itt = ind(l)
	ind(l) = ind(k)
	ind(k) = itt

	! find an element in the upper part of the array which is
	!   not larger than t

	50 l = l - 1
	indx = ind(l)
	if (x(indx) > t) go to 50

	! find an element in the lower part of the array whcih is not smaller than t

	60 k = k + 1
	indx = ind(k)
	if (x(indx) < t) go to 60

	! if k <= l, interchange elements k and l

	if (k <= l) go to 40

	! save the upper and lower subscripts of the portion of the
	!   array yet to be sorted

	if (l-i > j-k) then
	il(m) = i
	iu(m) = l
	i = k
	m = m + 1
	go to 80
	end if

	il(m) = k
	iu(m) = j
	j = l
	m = m + 1
	go to 80

	! begin again on another unsorted portion of the array

	70 m = m - 1
	if (m == 0) return
	i = il(m)
	j = iu(m)

	80 if (j-i >= 11) go to 30
	if (i == 1) go to 20
	i = i - 1

	! sort elements i+1,...,j.  note that 1 <= i < j and j-i < 11.

	90 i = i + 1
	if (i == j) go to 70
	indx = ind(i+1)
	t = x(indx)
	it = indx
	indx = ind(i)
	if (x(indx) <= t) go to 90
	k = i

	100 ind(k+1) = ind(k)
	k = k - 1
	indx = ind(k)
	if (t < x(indx)) go to 100

	ind(k+1) = it
	go to 90
	end subroutine qsortd

! ---------------------------------------------------------
! Given x(:) and y(:) which tabulate a function and the derivative at the boundary points
! this function returns the second derivative of the spline at each point
! ---------------------------------------------------------
		subroutine splin1(x,y,yp1,ypn,y2)
		real(kind=8), INTENT(IN) :: x(:), y(:), yp1, ypn
		real(kind=8), INTENT(INOUT) :: y2(size(x))
		integer :: n, i, k
		real(kind=8) :: p, qn, sig, un, u(size(x))

			n = size(x)

			if (yp1 > .99e30) then
				y2(1) = 0.e0
				u(1) = 0.e0
			else
				y2(1) = -0.5e0
				u(1) = (3.e0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
			endif

			do i = 2, n-1
				sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
				p = sig * y2(i-1)+2.e0
				y2(i) = (sig-1.e0)/p
				u(i) = (6.e0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/&
					(x(i+1)-x(i-1))-sig*u(i-1))/p
			enddo
			if (ypn > .99e30) then
				qn = 0.e0
				un = 0.e0
			else
				qn = 0.5e0
				un = (3.e0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
			endif

			y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.e0)

			do k = n-1, 1, -1
				y2(k) = y2(k)*y2(k+1)+u(k)
			enddo

		end subroutine splin1

! ---------------------------------------------------------
! Given xa(:) and ya(:) which tabulate a function, returns the interpolation using
! splines of vector x(:) in y(:)
! ---------------------------------------------------------
		subroutine spline(xa,ya,x,y)
		real(kind=8), INTENT(INOUT) :: y(:)
		real(kind=8), INTENT(IN) :: xa(:), ya(:), x(:)
		real(kind=8) :: y2a(size(xa))
		integer :: n_x, n, i, k, khi, klo
		real(kind=8) :: a, b, h, extrap

			n = size(xa)
			n_x = size(x)
			call splin1(xa,ya,huge(1.d0),huge(1.d0),y2a)

			do i = 1, n_x

! Downward extrapolation
				if (x(i) < xa(1)) then
!					y(i) = ya(1)
					y(i) = ya(1) + (ya(1)-ya(2))/(xa(1)-xa(2)) * (xa(1) - x(i))
				else

! Upward extrapolation
				if (x(i) > xa(n)) then
!					y(i) = ya(n)
					y(i) = ya(n) + (ya(n)-ya(n-1)) / (xa(n)-xa(n-1)) * (x(i) - xa(n))
				else
! In range
						klo = 1
						khi = n
1						if(khi-klo > 1) then
							k = (khi+klo)/2
							if (xa(k) > x(i)) then
								khi = k
							else
								klo = k
							endif
							go to 1
						endif

						h = xa(khi)-xa(klo)

						if (h == 0.d0) then
							print *, 'bad xa input in spline'
							stop
						endif
						a = (xa(khi)-x(i))/h
						b = (x(i)-xa(klo))/h

						y(i) = a*ya(klo)+b*ya(khi)+((a**3.e0-a)*y2a(klo)+(b**3.e0-b)*y2a(khi))*(h**2.e0)/6.e0
					endif
				endif
			enddo

		end subroutine spline
		
	subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
	integer n
	double precision a(n,n), c(n,n)
	double precision L(n,n), U(n,n), b(n), d(n), x(n)
	double precision coeff
	integer i, j, k

	! step 0: initialization for matrices L and U and b
	! Fortran 90/95 aloows such operations on matrices
	L=0.0
	U=0.0
	b=0.0

	! step 1: forward elimination
	do k=1, n-1
	do i=k+1,n
		coeff=a(i,k)/a(k,k)
		L(i,k) = coeff
		do j=k+1,n
			a(i,j) = a(i,j)-coeff*a(k,j)
		end do
	end do
	end do

	! Step 2: prepare L and U matrices 
	! L matrix is a matrix of the elimination coefficient
	! + the diagonal elements are 1.0
	do i=1,n
	L(i,i) = 1.0
	end do
	! U matrix is the upper triangular part of A
	do j=1,n
	do i=1,j
		U(i,j) = a(i,j)
	end do
	end do

	! Step 3: compute columns of the inverse matrix C
	do k=1,n
	b(k)=1.0
	d(1) = b(1)
	! Step 3a: Solve Ld=b using the forward substitution
	do i=2,n
		d(i)=b(i)
		do j=1,i-1
		d(i) = d(i) - L(i,j)*d(j)
		end do
	end do
	! Step 3b: Solve Ux=d using the back substitution
	x(n)=d(n)/U(n,n)
	do i = n-1,1,-1
		x(i) = d(i)
		do j=n,i+1,-1
		x(i)=x(i)-U(i,j)*x(j)
		end do
		x(i) = x(i)/u(i,i)
	end do
	! Step 3c: fill the solutions x(n) into column k of C
	do i=1,n
		c(i,k) = x(i)
	end do
	b(k)=0.0
	end do
	end subroutine inverse


end module maths