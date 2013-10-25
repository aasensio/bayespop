program main

use init, only : setProblem
use params, only : galaxy, sdim
use sampling_m, only : hmcSample
use maximumlikelihood_m, only : getMaximumLikelihood
use maths, only : init_random_seed

implicit none

	integer :: i
	
	call init_random_seed

	call setProblem
	

!	do galaxy%whichComputing = 1, galaxy%nSpec
		galaxy%whichComputing = 1
	

! Find the maximum a-posteriori values
 		call getMaximumLikelihood

! Hamiltonian Metropolis
		call hmcSample
		
!	enddo
   
end
