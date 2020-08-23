MODULE utilities_random
  IMPLICIT NONE
CONTAINS
  SUBROUTINE random1(n,sum,sumsq)
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(OUT) :: sum, sumsq
    
    INTEGER :: i
    REAL :: ram
    DO i=1,n
       CALL RANDOM_NUMBER(ram)
       sum=sum+ram
       sumsq=sumsq+ram*ram
    END DO
  END SUBROUTINE random1
    
  REAL FUNCTION average(sum,n)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: sum
    !    REAL, INTENT(OUT) :: average
    
    average=sum/REAL(n)

  END FUNCTION average

  REAL FUNCTION standard_derivation(sum,sumsq,n)
    IMPLICIT NONE
    REAL , INTENT(IN) :: sum, sumsq
    INTEGER, INTENT(IN) :: n
    !  REAL, INTENT(OUT) :: standard_derivation
    standard_derivation=sqrt((1/REAL(n))*sumsq - (1/REAL(n)*sum)**2)
    
  END FUNCTION STANDARD_DERIVATION
END MODULE utilities_random
