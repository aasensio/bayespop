module chainAnalysis_m
use params
use maths, only : int_tabulated, qsortd, spline
! use maths_mod
implicit none

contains

!------------------------------------------------------------
! Read the chains
!------------------------------------------------------------
	subroutine read_chains_multinest(fich_chain, chain_analysis)
	type(markov_t) :: chain_analysis
	character(len=100) :: fich_chain
	integer :: nparams, chain_length, i, j
	logical :: exists

		i = 0

		inquire(file=trim(adjustl(fich_chain))//'post_equal_weights.dat', exist=exists)
		if (.not.exists) then
			write(*,*) 'File not found : ', trim(adjustl(fich_chain))//'post_equal_weights.dat'
			stop
		endif

		open(unit=12,file=trim(adjustl(fich_chain))//'post_equal_weights.dat',action='read',status='old')

		do while (.true.)
   		read (12, *, end=999)
   		i = i + 1
		enddo
999   continue
		close(12)

		chain_analysis%nlength = i
		chain_analysis%nparam = nest_nPar
		print *, 'Length of posterior samples : ', chain_analysis%nlength
		print *, 'Number of parameters : ', chain_analysis%nparam

		allocate(chain_analysis%chain(chain_analysis%nparam+1,chain_analysis%nlength))
		allocate(chain_analysis%best_parameters(chain_analysis%nparam))
		allocate(chain_analysis%most_probable(chain_analysis%nparam))
		allocate(chain_analysis%param_names(chain_analysis%nparam))

		chain_analysis%param_names = (/'B   ', 'thB ', 'chiB' /)
		
		open(unit=12,file=trim(adjustl(fich_chain))//'post_equal_weights.dat',action='read',&
			status='old')
		do i = 1, chain_analysis%nlength
			read(12,*) (chain_analysis%chain(j,i),j=1,chain_analysis%nparam+1)
		enddo
		close(12)

! Read evidence
		inquire(file=trim(adjustl(fich_chain))//'stats.dat', exist=exists)
		if (.not.exists) then
			write(*,*) 'File not found : ', trim(adjustl(fich_chain))//'stats.dat'
			stop
		endif

		open(unit=12,file=trim(adjustl(fich_chain))//'stats.dat',action='read',status='old')
		read(12,FMT='(16X,E20.12)') chain_analysis%evidence
		close(12)

! Read the MAP parameters
		open(unit=12,file=trim(adjustl(fich_chain))//'.map.pars',action='read',status='old')
		read(12,*) chain_analysis%best_parameters
		close(12)


	end subroutine read_chains_multinest

!------------------------------------------------------------
! Get the optimal binning for building a histogram of a variable
! Ref: Freedman, D, & Diaconis, P. (1981). "On the histogram as a density
!      estimator: L2 theory". Zeitschrift für Wahrscheinlichkeitstheorie und
!      verwandte Gebiete 57 (4): 453–476
!  BIN=2*IQR(var)/n^(1/3)   with n the number of elements of the array 'var'
!------------------------------------------------------------
	function optbin(var)
	real(kind=8) :: optbin, var(:)
	integer :: n, quart1, quart3
	integer, allocatable :: ind(:)

		n = size(var)
		allocate(ind(n))

		call qsortd(var,ind,n)

		quart1 = n*0.25
		quart3 = n*0.75

		optbin = 2.0*(var(ind(quart3)) - var(ind(quart1))) / (n*1.0)**(1.0/3.0)

	end function optbin

!------------------------------------------------------------
! 1D histograms
!------------------------------------------------------------
	subroutine oned_histogram(chain, x, yGauss, yStep)
	real(kind=8) :: chain(:)
	real(kind=8) :: step, xmin, xmax
	real(kind=8), pointer :: x(:), yGauss(:), yStep(:)
	real(kind=8), allocatable :: xx(:)
	integer :: i, n, nn, ncut
	real(kind=8) :: wei, norm, error_norm, sig

! Bins
		step = optbin(chain)

		xmin = minval(chain)
		xmax = maxval(chain)
		n = (xmax-xmin) / step
		step = (xmax-xmin) / float(n)

! Variables
		allocate(x(n))
		allocate(yGauss(n))
		allocate(yStep(n))

		do i = 1, n
			x(i) = xmin + step * (i-1.0)
		enddo

! Variables for the normalization
		if (n <= 50) then
			nn = 50
			allocate(xx(nn))
			do i = 1, nn
				xx(i) = xmin + (xmax-xmin)/nn * i
			enddo
		else
			allocate(xx(n))
			xx = x
		endif

! Doing the histogram
		sig = 1.2d0
		do i = 1, n
			ncut = count(chain >= x(i)-step/2.d0 .and. chain < x(i)+step/2.0)

! Smoothing kernel
			wei = sum(exp(-0.50 * (chain-x(i))**2 / (sig*step)**2.0))
			norm = int_tabulated(xx, exp(-0.50 * (xx-x(i))**2 / (sig*step)**2.0))
			yGauss(i) = wei / norm
			yStep(i) = ncut
		enddo

		deallocate(xx)

		yGauss = yGauss / maxval(yGauss)
		yStep = yStep / maxval(yStep)

	end subroutine oned_histogram

!------------------------------------------------------------
! Confidence intervals
! INPUT
!   xdata : x position of the histogram
!   ydata : y position of the histogram
! OUTPUT
!   est : estimated value of the parameter
! errup : value of x at the upper confidence interval
! errdown : value of x at the lower confidence interval
! conf : desired confidence level
!------------------------------------------------------------
	subroutine conf_limits_ccf(xdata, ydata, conf, est, errup, errdown)
	real(kind=8) :: xdata(:), ydata(:)
	real(kind=8) :: conf, est, errup, errdown
	real(kind=8) :: xpdf_low, xpdf_up, xpdf_mid, xmin, xmax, lower, upper, norm
	real(kind=8), allocatable :: xx(:), yy(:), x(:), y(:), xF(:), F(:)
	integer, allocatable :: ind(:)
	integer :: n, npto, loc(1), i

		n = size(xdata)

		xpdf_low = 0.50 - conf/200.
		xpdf_mid = 0.5
		xpdf_up  = 0.5 + conf/200.

		if (xpdf_low < 0.0 .or. xpdf_up > 1) then
			print *,'wrong value for CONF'
			stop
		endif

		allocate(xx(n))
		allocate(yy(n))
		allocate(ind(n))

		xx = xdata
		yy = ydata / maxval(ydata)

! Sorting
		call qsortd(xx,ind,n)
		xx = xx(ind)
		yy = yy(ind)
		deallocate(ind)

! Interpolation by splines
		npto = 300
		xmin = xx(1)
		xmax = xx(n)

		if (npto > n) then
			allocate(x(npto))
			allocate(y(npto))

			do i = 1, npto
				x(i) = xmin + (xmax-xmin)/npto * i
			enddo
			call spline(xx,yy,x,y)
			n = npto
		else
			allocate(x(n))
			allocate(y(n))
			x = xx
			y = yy
		endif

! Peak normalization
		y = y / maxval(y)

! Computing the cumulative distribution function
		norm = int_tabulated(x,y)

		allocate(xF(0:n-1))
		allocate(F(0:n-1))
		xF = x(1:n)

		do i = 5, n
			F(i-1) = int_tabulated(x(1:i-1),y(1:i-1)) / norm
		enddo

		loc = minloc(abs(F - xpdf_low))
		lower = xF(loc(1))

		loc = minloc(abs(F - xpdf_mid))
		est = xF(loc(1))

		loc = minloc(abs(F - xpdf_up))
		upper = xF(loc(1))

		errup = upper - est
		errdown = est - lower

		deallocate(x,y)
		deallocate(xF,F)
		deallocate(xx,yy)

	end subroutine conf_limits_ccf

!------------------------------------------------------------
! Analyze the chains
!------------------------------------------------------------
	subroutine analyze_chains_multinest(chain_analysis)
	type(markov_t) :: chain_analysis
! 	type(prior_type) :: prior
	logical :: suggest
	real(kind=8), pointer :: x(:), yGauss(:), yStep(:)
	integer :: i, j, n
	character(len=120) :: str_parameter
	real(kind=8) :: est, errup1s, errup2s, errdown1s, errdown2s, chi2_mean, chi2_max, pars(3), pctg
	integer, pointer :: indx(:)
	real(kind=8), allocatable :: chain_function(:)

! 		chain_analysis%filename = nest_root
! 		
! 		call read_chains_multinest(chain_analysis%filename, chain_analysis)
! 
! 		open(unit=12,file=trim(adjustl(chain_analysis%filename))//'.hist1D',&
! 			action='write',status='replace')
! 		write(12,*) chain_analysis%nparam
! 
! 		open(unit=13,file=trim(adjustl(chain_analysis%filename))//'.confidence',&
! 			action='write',status='replace')
! 
! ! Do histograms for all inverted variables
! 		do i = 1, chain_analysis%nparam
! 
! 			if (prior%typ(i) /= 'D') then
! 				call oned_histogram(chain_analysis%chain(i,:), x, yGauss, yStep)
! 
! ! 1sigma confidence intervals
! 				call conf_limits_ccf(x, yStep, one_sigma, est, errup1s, errdown1s)
! 				call conf_limits_ccf(x, yStep, two_sigma, est, errup2s, errdown2s)
! 			else
! 				allocate(x(10))
! 				x = prior%mu(i)
! 				allocate(yGauss(10))
! 				yGauss = 1.0
! 				allocate(yStep(10))
! 				yStep = 1.0
! 				est = prior%mu(i)
! 				errup1s = 0.d0
! 				errup2s = 0.d0
! 				errdown1s = 0.d0
! 				errdown2s = 0.d0
! 			endif
! 
! 			chain_analysis%most_probable(i) = est
! 			write(13,*) est, chain_analysis%best_parameters(i), errdown1s, errup1s, errdown2s, errup2s
! 			str_parameter = chain_analysis%param_names(i)
! 			write(*,FMT='(A,A,F9.4,A7,F9.4)') trim(adjustl(str_parameter)), ' : E(x)=', est, ' - MAP=', &
! 				chain_analysis%best_parameters(i)
! 			write(*,FMT='(A,F9.4,A,F9.4)') '1-sigma : -', errdown1s, '  +', errup1s
! 			write(*,FMT='(A,F9.4,A,F9.4)') '2-sigma : -', errdown2s, '  +', errup2s
! 			n = size(x)
! 			write(12,*) i, n
! 			do j = 1, n
! 				write(12,*) x(j), yGauss(j), yStep(j)
! 			enddo
! 			deallocate(x)
! 			deallocate(yGauss)
! 			deallocate(yStep)
! 
! 		enddo
! 
! ! Calculate average ln L for calculating Kullback-Leibler distance
! 		chain_analysis%avg_lnL = sum(chain_analysis%chain(chain_analysis%nparam+1,:)) / chain_analysis%nlength
! 
! 		write(13,*) 'Evidence, <ln L> and Kullback-Leibler divergence'
! 		write(13,*) chain_analysis%evidence, chain_analysis%avg_lnL, &
! 			-chain_analysis%evidence + chain_analysis%avg_lnL
! 
! 		write(*,*) 'Evidence, <ln L> and Kullback-Leibler divergence'
! 		write(*,*) chain_analysis%evidence, chain_analysis%avg_lnL, &
! 			-chain_analysis%evidence + chain_analysis%avg_lnL
! 
! 		close(12)
! 		close(13)
! 
! 		open(unit=12,file=trim(adjustl(chain_analysis%filename))//'.stokesSamples',&
! 			action='write',status='replace',form='unformatted')
! 
! ! The two last are the MAP profiles and the estimated one
! 		write(12) chain_analysis%nlength
! 
! ! Save file with evaluated profiles at the posterior samples
! 		do i = 1, chain_analysis%nlength
! 			pars = chain_analysis%chain(1:3,i)
! 
! 			call lininterpolDatabase(pars, database%model)
! 			
! 			write(12) database%model
! 		enddo
! 
! 		close(12)

	end subroutine analyze_chains_multinest

end module chainAnalysis_m