module mpi_mod
use params
implicit none

contains
 	
!------------------------------------------------------------
! Send a new computation to a slave
!------------------------------------------------------------
	subroutine sendNewTask(slave, indComputation, indCase, indGalaxy)
	integer :: slave, indComputation, indCase, indGalaxy
	integer :: ierr, i
	include 'mpif.h'
	
! Send this to continue computing
		call MPI_Send(0, 1, MPI_INTEGER, slave, 13, MPI_COMM_WORLD, ierr)
		
		call MPI_Send(indComputation, 1, MPI_INTEGER, slave, 13, MPI_COMM_WORLD, ierr)
		call MPI_Send(indCase, 1, MPI_INTEGER, slave, 13, MPI_COMM_WORLD, ierr)
		call MPI_Send(indGalaxy, 1, MPI_INTEGER, slave, 13, MPI_COMM_WORLD, ierr)
		
! Send the priors
		do i = 1, 6
			call MPI_Send(priors(i)%typ, 1, MPI_INTEGER, slave, 13, MPI_COMM_WORLD, ierr)
			call MPI_Send(priors(i)%lower, 1, MPI_DOUBLE_PRECISION, slave, 13, MPI_COMM_WORLD, ierr)
			call MPI_Send(priors(i)%upper, 1, MPI_DOUBLE_PRECISION, slave, 13, MPI_COMM_WORLD, ierr)
			call MPI_Send(priors(i)%mu, 1, MPI_DOUBLE_PRECISION, slave, 13, MPI_COMM_WORLD, ierr)
			call MPI_Send(priors(i)%sigma, 1, MPI_DOUBLE_PRECISION, slave, 13, MPI_COMM_WORLD, ierr)
		enddo
		
		call MPI_Send(library%nPixAdapted, 1, MPI_INTEGER, slave, 13, MPI_COMM_WORLD, ierr)
		call MPI_Send(library%nSpec, 1, MPI_INTEGER, slave, 13, MPI_COMM_WORLD, ierr)
		call MPI_Send(library%velScale, 1, MPI_DOUBLE_PRECISION, slave, 13, MPI_COMM_WORLD, ierr)
		
		call MPI_Send(library%specAdapted, library%nPixAdapted * library%nSpec, MPI_DOUBLE_PRECISION, slave, 13, MPI_COMM_WORLD, ierr)
		call MPI_Send(library%age, library%nSpec, MPI_DOUBLE_PRECISION, slave, 13, MPI_COMM_WORLD, ierr)
		call MPI_Send(library%metallicity, library%nSpec, MPI_DOUBLE_PRECISION, slave, 13, MPI_COMM_WORLD, ierr)
		call MPI_Send(library%imfSlope, library%nSpec, MPI_DOUBLE_PRECISION, slave, 13, MPI_COMM_WORLD, ierr)
		call MPI_Send(library%mgfe, library%nSpec, MPI_DOUBLE_PRECISION, slave, 13, MPI_COMM_WORLD, ierr)
		
		call MPI_Send(galaxy%snr, 1, MPI_DOUBLE_PRECISION, slave, 13, MPI_COMM_WORLD, ierr)
		
		call MPI_Send(galaxy%spec, galaxy%nPix, MPI_DOUBLE_PRECISION, slave, 13, MPI_COMM_WORLD, ierr)
		call MPI_Send(galaxy%maskIndex, galaxy%nPix, MPI_INTEGER, slave, 13, MPI_COMM_WORLD, ierr)
		call MPI_Send(galaxy%fixVLOS, 1, MPI_INTEGER, slave, 13, MPI_COMM_WORLD, ierr)
		
 	end subroutine sendNewTask
 	
!------------------------------------------------------------
! Receive a new computation from the master
!------------------------------------------------------------
	subroutine receiveNewTask(killFlag, indComputation, indCase, indGalaxy)
	integer :: killFlag, indComputation, indCase, indGalaxy
	include 'mpif.h'
	integer :: ierr, i, status(MPI_STATUS_SIZE)	
		
		call MPI_Recv(killFlag, 1, MPI_INTEGER, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)		
		
		if (killFlag == 0) then
		
			call MPI_Recv(indComputation, 1, MPI_INTEGER, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
			call MPI_Recv(indCase, 1, MPI_INTEGER, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
			call MPI_Recv(indGalaxy, 1, MPI_INTEGER, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
		
! Send the priors
			do i = 1, 6
				call MPI_Recv(priors(i)%typ, 1, MPI_INTEGER, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
				call MPI_Recv(priors(i)%lower, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
				call MPI_Recv(priors(i)%upper, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
				call MPI_Recv(priors(i)%mu, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
				call MPI_Recv(priors(i)%sigma, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
			enddo
			
			call MPI_Recv(library%nPixAdapted, 1, MPI_INTEGER, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
			call MPI_Recv(library%nSpec, 1, MPI_INTEGER, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
			call MPI_Recv(library%velScale, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
			
			galaxy%nPix = library%nPixAdapted
			
			if (associated(library%specAdapted)) deallocate(library%specAdapted)
			if (associated(library%age)) deallocate(library%age)
			if (associated(library%metallicity)) deallocate(library%metallicity)
			if (associated(library%imfSlope)) deallocate(library%imfSlope)
			if (associated(library%mgfe)) deallocate(library%mgfe)
			
			allocate(library%specAdapted(library%nPixAdapted,library%nSpec))
			allocate(library%age(library%nSpec))
			allocate(library%metallicity(library%nSpec))
			allocate(library%imfSlope(library%nSpec))
			allocate(library%mgfe(library%nSpec))
			
			call MPI_Recv(library%specAdapted, library%nPixAdapted * library%nSpec , MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
			call MPI_Recv(library%age, library%nSpec , MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
			call MPI_Recv(library%metallicity, library%nSpec , MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
			call MPI_Recv(library%imfSlope, library%nSpec , MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
			call MPI_Recv(library%mgfe, library%nSpec , MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
			
			call MPI_Recv(galaxy%snr, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
			
			if (associated(galaxy%spec)) deallocate(galaxy%spec)
			if (associated(galaxy%maskIndex)) deallocate(galaxy%maskIndex)
			if (associated(galaxy%synth)) deallocate(galaxy%synth)
			if (associated(galaxy%synth2)) deallocate(galaxy%synth2)
			if (associated(galaxy%convolvedSpectrum)) deallocate(galaxy%convolvedSpectrum)
			
			allocate(galaxy%spec(galaxy%nPix))
			allocate(galaxy%maskIndex(galaxy%nPix))
			allocate(galaxy%synth(galaxy%nPix))
			allocate(galaxy%synth2(galaxy%nPix))
			allocate(galaxy%convolvedSpectrum(galaxy%nPix))
			
			call MPI_Recv(galaxy%spec, galaxy%nPix , MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
			call MPI_Recv(galaxy%maskIndex, galaxy%nPix , MPI_INTEGER, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
			call MPI_Recv(galaxy%fixVLOS, 1 , MPI_INTEGER, MPI_ANY_SOURCE, 13, MPI_COMM_WORLD, status, ierr)
			
! Fix all velocities to a common value
			if (galaxy%fixVLOS == 0) then
				sdim = 3*library%nSpec
			else
				sdim = library%nSpec+2
			endif
			
			nest_nPar = sdim

! Allocate memory for gradient of spectrum used in the calculation of the derivative of the log-posterior
			if (associated(galaxy%synthGrad)) deallocate(galaxy%synthGrad)
			allocate(galaxy%synthGrad(galaxy%nPix,sdim))

			if (allocated(map_pars)) deallocate(map_pars)
			allocate(map_pars(sdim))		

! For periodic boundary conditions
			if (allocated(pWrap)) deallocate(pWrap)
			allocate(pWrap(sdim))
			pWrap = 0

			if (associated(galaxy%trialWeight)) deallocate(galaxy%trialWeight)
			if (associated(galaxy%trialVelocity)) deallocate(galaxy%trialVelocity)
			if (associated(galaxy%trialDispersion)) deallocate(galaxy%trialDispersion)
			
			allocate(galaxy%trialWeight(library%nSpec))
			allocate(galaxy%trialVelocity(library%nSpec))
			allocate(galaxy%trialDispersion(library%nSpec))
					
			if (associated(galaxy%trial)) deallocate(galaxy%trial)
			allocate(galaxy%trial(sdim))

			maxslhood = -1.d100		
		endif
			
 	end subroutine receiveNewTask

!------------------------------------------------------------
! Receive result from slave
!------------------------------------------------------------
	subroutine receiveNewResult(slave, computationOK, indCase, indGalaxy)
	integer :: slave, indCase, indGalaxy
	logical :: computationOK
	include 'mpif.h'
	integer :: ierr, i, status(MPI_STATUS_SIZE)
	character(len=4) :: strGalaxy
	integer :: nSamples, nPix, nPars
	real(kind=8), allocatable :: samples(:,:)
	
		call MPI_Recv(slave, 1, MPI_INTEGER, MPI_ANY_SOURCE, 14, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(computationOK, 1, MPI_LOGICAL, MPI_ANY_SOURCE, 15, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(indCase, 1, MPI_INTEGER, MPI_ANY_SOURCE, 16, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(indGalaxy, 1, MPI_INTEGER, MPI_ANY_SOURCE, 17, MPI_COMM_WORLD, status, ierr)
		
		call MPI_Recv(nSamples, 1, MPI_INTEGER, MPI_ANY_SOURCE, 18, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(nPars, 1, MPI_INTEGER, MPI_ANY_SOURCE, 19, MPI_COMM_WORLD, status, ierr)
						
		if (allocated(samples)) deallocate(samples)
		allocate(samples(nPars,nSamples))
				
		call MPI_Recv(samples, nPars*nSamples, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 20, MPI_COMM_WORLD, status, ierr)				
		
 		write(strGalaxy, FMT='(I4)') indGalaxy
 		open(unit=12,file='output/'//trim(adjustl(fileCases(indCase)))//"_"//trim(adjustl(strGalaxy))//".samples",form='unformatted',action='write',access='stream')
 		write(12) nSamples, nPars
 		do i = 1, nSamples			
 			write(12) samples(:,i)
 		enddo
 		close(12)
		 				
		call MPI_Recv(nPix, 1, MPI_INTEGER, MPI_ANY_SOURCE, 21, MPI_COMM_WORLD, status, ierr)		
				
		if (allocated(samples)) deallocate(samples)
		allocate(samples(nPix,nSamples))
				
		call MPI_Recv(samples, nPix*nSamples, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 22, MPI_COMM_WORLD, status, ierr)
		open(unit=12,file='output/'//trim(adjustl(fileCases(indCase)))//"_"//trim(adjustl(strGalaxy))//".synth",form='unformatted',action='write',access='stream')
		write(12) nSamples, nPix
		do i = 1, nSamples			
			write(12) samples(:,i)
		enddo
		close(12)
		
! Best fit
		call MPI_Recv(library%nSpec, 1, MPI_INTEGER, MPI_ANY_SOURCE, 23, MPI_COMM_WORLD, status, ierr)
		
		if (associated(galaxy%trialWeight)) deallocate(galaxy%trialWeight)
		if (associated(galaxy%trialVelocity)) deallocate(galaxy%trialVelocity)
		if (associated(galaxy%trialDispersion)) deallocate(galaxy%trialDispersion)
		if (associated(galaxy%synth)) deallocate(galaxy%synth)
		
		allocate(galaxy%trialWeight(library%nSpec))
		allocate(galaxy%trialVelocity(library%nSpec))
		allocate(galaxy%trialDispersion(library%nSpec))
		allocate(galaxy%synth(nPix))
		
		call MPI_Recv(galaxy%trialWeight, library%nSpec, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 24, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(galaxy%trialVelocity, library%nSpec, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 25, MPI_COMM_WORLD, status, ierr)
		call MPI_Recv(galaxy%trialDispersion, library%nSpec, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 26, MPI_COMM_WORLD, status, ierr)
				
		call MPI_Recv(galaxy%synth, nPix, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 27, MPI_COMM_WORLD, status, ierr)
				
		open(unit=12,file='output/'//trim(adjustl(fileCases(indCase)))//"_"//trim(adjustl(strGalaxy))//".MAP",form='unformatted',action='write',access='stream')
		write(12) library%nSpec, nPix		
		write(12) galaxy%trialWeight
		write(12) galaxy%trialVelocity
		write(12) galaxy%trialDispersion
		write(12) galaxy%synth
		close(12)
		
		call MPI_Recv(nItersLM, 1, MPI_INTEGER, MPI_ANY_SOURCE, 28, MPI_COMM_WORLD, status, ierr)
		if (allocated(samples)) deallocate(samples)
		allocate(samples(2,nItersLM))
		
		call MPI_Recv(samples, 2*nItersLM, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, 29, MPI_COMM_WORLD, status, ierr)
		open(unit=12,file='output/'//trim(adjustl(fileCases(indCase)))//"_"//trim(adjustl(strGalaxy))//".LM",form='unformatted',action='write',access='stream')
		write(12) nItersLM
		write(12) samples
		close(12)
						
	end subroutine receiveNewResult
	
!------------------------------------------------------------
! Send result to master
!------------------------------------------------------------
	subroutine sendNewResult(slave, computationOK, indCase, indGalaxy)
	integer :: slave, indCase, indGalaxy
	logical :: computationOK
	include 'mpif.h'
	integer :: ierr, i, status(MPI_STATUS_SIZE), j
	integer :: nSamples, nPix
	real(kind=8), allocatable :: samples(:,:), temp(:)
	
		call MPI_Send(slave, 1, MPI_INTEGER, 0, 14, MPI_COMM_WORLD, ierr)
		call MPI_Send(computationOK, 1, MPI_LOGICAL, 0, 15, MPI_COMM_WORLD, ierr)
		call MPI_Send(indCase, 1, MPI_INTEGER, 0, 16, MPI_COMM_WORLD, ierr)
		call MPI_Send(indGalaxy, 1, MPI_INTEGER, 0, 17, MPI_COMM_WORLD, ierr)
		
! Send the samples for the parameters
		flPfx = 'temp/test'//trim(adjustl(myrankStr))
		open(unit=20,file= (trim(flPfx)//".samples"),form='unformatted',action='read',access='stream')
		read(20) nSamples, library%nSpec
		
		if (allocated(samples)) deallocate(samples)
		if (allocated(temp)) deallocate(temp)
		allocate(samples(3*library%nSpec,nSamples))
		allocate(temp(3*library%nSpec))
		
		do i = 1, nSamples
			read(20) temp
			samples(:,i) = temp
		enddo
		close(20, status="delete")
				
		call MPI_Send(nSamples, 1, MPI_INTEGER, 0, 18, MPI_COMM_WORLD, ierr)
		call MPI_Send(library%nSpec, 1, MPI_INTEGER, 0, 19, MPI_COMM_WORLD, ierr)
		call MPI_Send(samples, 3*library%nSpec*nSamples, MPI_DOUBLE_PRECISION, 0, 20, MPI_COMM_WORLD, ierr)
						
		nPix = galaxy%nPix
		
! Send the samples for the spectra
		open(unit=20,file= (trim(flPfx)//".spectra"),form='unformatted',action='read',access='stream')
		read(20) nPix
		
		if (allocated(samples)) deallocate(samples)
		if (allocated(temp)) deallocate(temp)
		allocate(samples(nPix,nSamples))
		allocate(temp(nPix))
		do i = 1, nSamples
			read(20) temp
			samples(:,i) = temp
		enddo
		close(20, status="delete")
				
		call MPI_Send(nPix, 1, MPI_INTEGER, 0, 21, MPI_COMM_WORLD, ierr)
		call MPI_Send(samples, galaxy%nPix*nSamples, MPI_DOUBLE_PRECISION, 0, 22, MPI_COMM_WORLD, ierr)
		
! Best fit
		open(unit=15,file='temp/bestFitPars'//trim(adjustl(myrankStr))//'.dat',action='read',status='old')
  		do i = 1, library%nSpec
			read(15,*) galaxy%trialWeight(i), galaxy%trialVelocity(i), galaxy%trialDispersion(i)
		enddo
		close(15, status="delete")
		
		call MPI_Send(library%nSpec, 1, MPI_INTEGER, 0, 23, MPI_COMM_WORLD, ierr)		
		call MPI_Send(galaxy%trialWeight, library%nSpec, MPI_DOUBLE_PRECISION, 0, 24, MPI_COMM_WORLD, ierr)
		call MPI_Send(galaxy%trialVelocity, library%nSpec, MPI_DOUBLE_PRECISION, 0, 25, MPI_COMM_WORLD, ierr)
		call MPI_Send(galaxy%trialDispersion, library%nSpec, MPI_DOUBLE_PRECISION, 0, 26, MPI_COMM_WORLD, ierr)
		
		open(unit=15,file='temp/bestFitSpec'//trim(adjustl(myrankStr))//'.dat',action='read',status='old')
		read(15,*) nPix
		do i = 1, nPix
			read(15,*) galaxy%synth(i)
		enddo
		close(15, status="delete")
		
		call MPI_Send(galaxy%synth, galaxy%nPix, MPI_DOUBLE_PRECISION, 0, 27, MPI_COMM_WORLD, ierr)
		
		if (allocated(samples)) deallocate(samples)
		allocate(samples(2,nItersLM))
		call MPI_Send(nItersLM, 1, MPI_INTEGER, 0, 28, MPI_COMM_WORLD, ierr)
		
		open(unit=15,file='temp/LM_'//trim(adjustl(myrankStr))//'.iterations',action='read',status='old')
		do i = 1, 7
			read(15,*)
		enddo
		do i = 1, nItersLM
			read(15,*) (samples(j,i),j=1,2)			
		enddo
		call MPI_Send(samples, 2*nItersLM, MPI_DOUBLE_PRECISION, 0, 29, MPI_COMM_WORLD, ierr)
		close(15, status="delete")
				
	end subroutine sendNewResult

!------------------------------------------------------------
! Send the kill signal to all slaves
!------------------------------------------------------------
	subroutine killSlave(slave)
	integer :: slave
	integer :: ierr, i
	include 'mpif.h'

! Send a message with the killing flag activated
		call MPI_Send(1, 1, MPI_INTEGER, slave, 13, MPI_COMM_WORLD, ierr)

 	end subroutine killSlave
 	
end module mpi_mod