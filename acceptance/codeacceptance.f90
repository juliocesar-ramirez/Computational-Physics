
PROGRAM main
  USE utilities_acceptance
  IMPLICIT NONE
  REAL :: avg,sigma,ck,a1,b1,integral,integral1
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
   !   WRITE(*,*) q
   WRITE(*,*) avg, sigma

   CALL autocorrelation(q,n,avg,sigma,ck)
   WRITE(*,*)"ck: ", ck

   a1=0
   b1=3
   CALL montecarlos(q,n,a1,b1,integral,fx)
   WRITE(*,*) "Monte CArlos", integral

   CALL acceptance_rejection_method(q,n,b1,integral1)
   WRITE(*,*) "acceptancex",integral1
 END PROGRAM main


