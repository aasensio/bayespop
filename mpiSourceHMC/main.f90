program main

use init, only : setProblem
use params, only : galaxy, sdim, nCases, nGalaxies, whichCase, whichGalaxy, fileCases, priors, myrankStr, myrank
use sampling_m, only : hmcSample
use maximumlikelihood_m, only : getMaximumLikelihood
use maths, only : init_random_seed
use mpi_mod, only : sendNewTask, sendNewResult, receiveNewTask, receiveNewResult, killSlave

implicit none

	integer :: mpi_status, nprocs, ierr, i, killFlag, slave, packagesize, j, loop
	integer, allocatable :: slaveActive(:)
	logical :: computationOK
	character(len=8) :: date
	character(len=10) :: time
	character(len=5) :: zone
	integer :: values(8), indCase, indGalaxy
	
	include 'mpif.h'

	call MPI_INIT(mpi_status)
	call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, mpi_status)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, mpi_status)

	if (myrank == 0) then		
		write(*,FMT='(A,I3)') 'Number of nodes: ', nprocs
		write(*,FMT='(A,I3,A,I3,A)') 'Node ', myrank, '/', nprocs, ' - MASTER'
	else	
		write(*,FMT='(A,I3,A,I3,A)') 'Node ', myrank, '/', nprocs, ' - SLAVE'
		write(myrankStr,FMT='(I4)') myrank		
	endif
	
	allocate(slaveActive(nprocs-1))
	slaveActive = 1
	
	call init_random_seed
	
! Read the configuration file
	if (myrank == 0) then
		open(unit=12,file='preproc_data/fileList.dat',action='read',status='old')
		read(12,*) nCases
		allocate(nGalaxies(nCases))
		allocate(fileCases(nCases))
		do i = 1, nCases
			read(12,*) nGalaxies(i)
			read(12,FMT='(A)') fileCases(i)
			fileCases(i) = fileCases(i)//'.input'
		enddo
		close(12)
		allocate(whichCase(sum(nGalaxies)))
		allocate(whichGalaxy(sum(nGalaxies)))
		loop = 1
		do i = 1, nCases
			do j = 1, nGalaxies(i)
				whichCase(loop) = i
				whichGalaxy(loop) = j
				loop = loop + 1
			enddo
		enddo
	endif
	
! Broadcast the number of cases and the names of the files
	call MPI_Bcast(nCases,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	if (myrank /= 0) then
		allocate(nGalaxies(nCases))
		allocate(fileCases(nCases))
	endif
	call MPI_Bcast(nGalaxies,nCases,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)	
	call MPI_Bcast(fileCases,nCases*len(fileCases(1)),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
	
		
	call MPI_Barrier(MPI_COMM_WORLD, ierr)
	
	
	allocate(priors(6))
	
! If running in parallel mode and at least one slave is available
	if (nprocs >= 2) then
	
! Initial distribution of works
		if (myrank == 0) then

			do loop = 1, min(sum(nGalaxies), nprocs-1)
				slave = loop				
				call setProblem(fileCases(whichCase(loop)), whichGalaxy(loop))
				call sendNewTask(slave, loop, whichCase(loop), whichGalaxy(loop))
				call date_and_time(date, time, zone, values)
				write(*,FMT='(A10,A,I5,A,I5,A,I5)') time, ' MASTER -> SLAVE ', slave, ' --- Galaxy ', loop,' /', sum(nGalaxies)
			enddo
			
		endif
		
		killFlag = 0

! Work with the remaining chunks
		do while(killFlag == 0)

!-------------------
! MASTER
!-------------------
			if (myrank == 0) then

! Receive new result and save it				
				call receiveNewResult(slave, computationOK, indCase, indGalaxy)
								
! New computation
				if (loop <= sum(nGalaxies)) then					
					call setProblem(fileCases(whichCase(loop)), whichGalaxy(loop))
					call sendNewTask(slave, loop, whichCase(loop), whichGalaxy(loop))					
					call date_and_time(date, time, zone, values)
					write(*,FMT='(A10,A,I5,A,I5,A,I5)') time, ' MASTER -> SLAVE ', slave, ' --- Galaxy ', loop,' /', sum(nGalaxies)
					loop = loop + 1
				else					
					call killSlave(slave)
					slaveActive(slave) = 0
					if (sum(slaveActive) == 0) then
						killFlag = 1
					endif
					call date_and_time(date, time, zone, values)
					write(*,FMT='(A10,A,I5,A)') time, ' MASTER -> SLAVE ', slave, ' --- Killed'
				endif
			endif

!-------------------
! SLAVE
!-------------------
			if (myrank /= 0) then

! Receive new chunk
				call receiveNewTask(killFlag, loop, indCase, indGalaxy)

				if (killFlag == 0) then

! Do the computation					
					call getMaximumLikelihood
 					call hmcSample
															
					call sendNewResult(myrank, .TRUE., indCase, indGalaxy)
					
					call date_and_time(date, time, zone, values)
					write(*,FMT='(A10,A,I5,A,I5,A,I5)') time, ' SLAVE ', myrank,' -> MASTER --- Finished ', loop,' /', sum(nGalaxies)
				endif
				
			endif
			
		enddo

	
	endif

! 	call setProblem
	

!	do galaxy%whichComputing = 1, galaxy%nSpec
! 		galaxy%whichComputing = 1
	

! Find the maximum a-posteriori values
!   		call getMaximumLikelihood

! Hamiltonian Metropolis
! 		call hmcSample
		
!	enddo

! Finalize MPI
	call MPI_FINALIZE(mpi_status)
   
end
