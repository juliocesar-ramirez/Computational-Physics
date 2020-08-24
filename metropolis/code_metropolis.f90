PROGRAM code_metropolis
  USE utilities
  IMPLICIT NONE

  INTEGER :: L,iteration
  INTEGER,ALLOCATABLE, DIMENSION(:,:) :: spin
  REAL :: energy, magnatize
  WRITE(*,*) "Enter the number of lattice points in one dimension"
  READ (*,*) L

  WRITE(*,*) "Enter the number of iteration"
  READ (*,*) iteration

  ALLOCATE(spin(L,L))
  
  CALL initial_values(L,spin)
  

  CALL i_magnetization_and_energy(spin,L,energy,magnatize)
  CALL mt_metropolis(iteration,spin,L,energy,magnatize)
  
  
END PROGRAM code_metropolis
