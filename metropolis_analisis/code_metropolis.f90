PROGRAM code_metropolis
  USE utilities
  IMPLICIT NONE

  
  INTEGER :: L,iteration
  INTEGER,ALLOCATABLE, DIMENSION(:,:) :: spin
  REAL :: Esum, M

  
  WRITE(*,*) "Enter the number of lattice points in one dimension"
  READ (*,*) L

  WRITE(*,*) "Enter the number of iteration"
  READ (*,*) iteration

  ALLOCATE(spin(L,L))
  
  CALL initial_values(L,spin)
  
  Esum=0
  M=0
  CALL i_magnetization_and_energy(spin,L,Esum,M)
  CALL mt_metropolis(iteration,spin,L,Esum,M)
  
  
END PROGRAM code_metropolis
