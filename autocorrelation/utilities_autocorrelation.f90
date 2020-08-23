MODULE utilities_autocorrelation
  IMPLICIT NONE
CONTAINS
  
  REAL FUNCTION average(q,n)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    REAL,DIMENSION(n), INTENT(IN) :: q

    INTEGER :: i
    REAL :: sum=0

    DO i=i, n
       sum=sum+q(i)
    END DO
    average=sum/REAL(n)

  END FUNCTION average
  !------------------------------------------------------------------------
  REAL FUNCTION standard_derivation(q,n)
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: n
    REAL, DIMENSION(n), INTENT(IN) :: q
    
    REAL  :: sum, sumsq
    INTEGER :: i
    DO i=0,n
       sum=sum+q(i)
       sumsq=sumsq+q(i)**2
    END DO
    standard_derivation=sqrt((1/REAL(n))*sumsq - (1/REAL(n)*sum)**2)

  END FUNCTION STANDARD_DERIVATION
  !------------------------------------------------------------------------
  SUBROUTINE autocorrelation(q, n,avg,sigma,ck)

    INTEGER, INTENT(IN) :: n
    REAL, DIMENSION(n), INTENT(IN) :: q
    REAL :: ck
    REAL :: avg,sigma
    INTEGER :: i,k,j
    REAL :: sum1
    
    OPEN(UNIT=8, FILE="CK",ACTION="WRITE",STATUS="OLD")
    DO j=1, 30
       k=j-1
       sum1=0
       DO i=1,n-k
          sum1=sum1 + (q(i)*q(i+k))
       END DO
       sum1=sum1/(real(n-k))
       ck=(sum1-(avg*avg))/(sigma*sigma)
       WRITE(8,*) k,ck
 END DO
END SUBROUTINE autocorrelation



END MODULE utilities_autocorrelation
