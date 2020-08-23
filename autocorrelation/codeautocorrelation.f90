
PROGRAM main
  USE utilities_autocorrelation
  IMPLICIT NONE
  REAL :: r,sum=0,sumsq=0,avg,sigma,ck
  INTEGER  :: i, n 
  REAL ,ALLOCATABLE, DIMENSION(:) :: q




  !  CALL RANDOM_SEED
  WRITE (*,*) "La cantidad de numeros random a crear"
  READ (*,*) n

  ALLOCATE(q(n))

  CALL RANDOM_NUMBER(q)
  !    CALL random1(q)

  avg = average(q,n)
  sigma=standard_derivation(q,n)
  WRITE(*,*) q
  WRITE(*,*) avg, sigma

  CALL autocorrelation(q,n,avg,sigma,ck)
  WRITE(*,*) ck

  END PROGRAM main


