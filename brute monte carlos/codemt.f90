
PROGRAM main
  USE utilities_mt
  IMPLICIT NONE
  REAL :: avg,sigma,ck,a1,b1,integral
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
   WRITE(*,*)"ck: ", ck

   a1=0
   b1=3
   CALL montecarlos(q,n,a1,b1,integral,fx)
   WRITE(*,*) "Monte CArlos", integral

  END PROGRAM main


