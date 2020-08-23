PROGRAM main
  USE utilities_random
  IMPLICIT NONE
  REAL :: r,sum=0,sumsq=0,avg,sigma
  INTEGER  :: i, n 



  !  CALL RANDOM_SEED
  WRITE (*,*) "La cantidad de numeros random a crear"
  READ (*,*) n

  CALL random1(n,sum,sumsq)
  avg = average(sum,n)
  sigma=standard_derivation(sum,sumsq,n)
  WRITE(*,*) n
  WRITE(*,*) avg, sigma
END PROGRAM main
  
  
