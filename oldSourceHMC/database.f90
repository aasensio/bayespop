module database_m
use maths, only : hunt
use params, only : database_t, database
use chk_g95_m
implicit none

contains

!-------------------------------------------------------------
! Read the database and precalculate some quantities
!-------------------------------------------------------------
	subroutine readDatabase
	integer :: i, nr, c1, j, k, c2, nmax
	character(len=20) :: access, form

! Test if the database is the original one or the convolved one
		if (index(database%filename, 'convolved') /= 0) then
			database%nlambda = 41
		else
			database%nlambda = 62
		endif
		
		database%ndim = 3
		allocate(database%indi(database%ndim+1))
		allocate(database%ee(database%ndim,2**database%ndim))
		allocate(database%wrk(2*database%nlambda,2**database%ndim))

! Do some precomputations for accelerating the linear interpolation routines on the database
		do i = 1, database%ndim
			database%indi(i) = database%ndim - i + 1
		enddo

		do i = 1, database%ndim
			nr = 2**(i-1)
			c1 = 1
			do j = 1, nr
				c2 = 1
				do k = 1, 2**database%ndim/nr
					database%ee(database%ndim-database%indi(i)+1,c1) = (c2-1) / 2**(database%ndim-i)
					c1 = c1 + 1
					c2 = c2 + 1
				enddo
			enddo
		enddo

! Make array of parameters
		nmax = 73
		allocate(database%params(database%ndim,nmax))
		allocate(database%B(51))
		allocate(database%thB(37))
		allocate(database%chiB(73))

		database%params = 1.d10

! B
		do i = 1, 51
			database%params(1,i) = (i-1.d0) * 5
			database%B(i) = (i-1.d0) * 5
		enddo

! thB
		do i = 1, 37
			database%params(2,i) = (i-1.d0) * 5
			database%thB(i) = (i-1.d0) * 5
		enddo

! chiB
		do i = 1, 73
			database%params(3,i) = (i-1.d0) * 5
			database%chiB(i) = (i-1.d0) * 5
		enddo

! The database
		allocate(database%lambda(database%nlambda))
		allocate(database%stokesI(database%nlambda))
		allocate(database%stokesQU(2*database%nlambda,73,37,51))
		allocate(database%model(2*database%nlambda))

! Read the database
		write(*,*) 'Reading database ', database%filename

! 
		if(is_g95()) then
	      access = 'stream'
	      form = 'unformatted'
	   else                    ! Use values appropriate for lf95 express 7.10.02 or
                           ! ifort Package ID: W_FC_C_9.0.029
	      access = 'sequential'
	      form = 'binary'
	   endif
   
		open(unit=12,file=database%filename,action='read',status='old',form=form,access=access)
		read(12) database%lambda
		read(12) database%stokesI
		read(12) database%stokesQU
		close(12)
		write(*,*) 'Done.'

! Compute Q/I and U/I
		do i = 1, 73
			do j = 1, 37
				do k = 1, 51
					database%stokesQU(1:database%nlambda,i,j,k) = database%stokesQU(1:database%nlambda,i,j,k) / database%stokesI
					database%stokesQU(database%nlambda+1:2*database%nlambda,i,j,k) = &
						database%stokesQU(database%nlambda+1:2*database%nlambda,i,j,k) / database%stokesI
				enddo
			enddo
		enddo

	end subroutine readDatabase
	
!-------------------------------------------------------------------
! Linear interpolation in database
! Note that the parameters in the code are in a different order to
! those on the database
! DB: chiB, thB, B
! code: B, thB, chiB
! and this modification has to be done before calling this routine
!-------------------------------------------------------------------
	subroutine lininterpolDatabase(pars, stokesOut)
	real(kind=8) :: pars(:), stokesOut(:), delta
	integer :: i, j, near(3), indices(3), ind

! Find the indices of the hypercube around the desired value
		call hunt(database%B, 51, pars(1), near(1))
		call hunt(database%thB, 37, pars(2), near(2))
		call hunt(database%chiB, 73, pars(3), near(3))

! Extract the values of the function that will be used
		do i = 1, 2**database%ndim
			indices = near + database%ee(:,i)
			database%wrk(:,i) = database%stokesQU(:,indices(3),indices(2),indices(1))
		enddo

! Do the actual linear interpolation
		do i = 1, database%ndim
			ind = database%indi(i)

			delta = -(pars(ind) - database%params(ind,near(ind))) / &
				(database%params(ind,near(ind)) - database%params(ind,near(ind)+1))

			do j = 1, 2**(database%ndim-i)

				database%wrk(:,j) = (database%wrk(:,2*j) - database%wrk(:,2*j-1)) * delta + database%wrk(:,2*j-1)

			enddo

		enddo

! Return Q/I and U/I
		stokesOut = database%wrk(:,1)

	end subroutine lininterpolDatabase

end module database_m