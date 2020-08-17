PROGRAM composite_trapezoidal
  IMPLICIT NONE
  REAL :: a=0,b=1,x=0,h=0
  REAL :: height,sum=0,fun,integral
  INTEGER :: i,n

  
  WRITE(*,*) "Númefro de divisiones"
  READ(*,*) n
  CALL imprimir

  h = height(n,a,b)
  DO i=1,n-1 ,1
     x=a+h*i
     sum=sum+fun(x)
  END DO
  integral=h/2 *(2*(sum)+fun(a)+fun(b))

  WRITE(*,*) integral
  
  
END PROGRAM composite_trapezoidal
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
SUBROUTINE imprimir
  WRITE(*,100)
  WRITE(*,100)
  WRITE(*,110) "COMPOSITE TRAPEZOIDAL RULE"
  WRITE(*,100)

100 FORMAT ("======================================================================")  
110 FORMAT (A40)  

END SUBROUTINE imprimir
!---------------------------------------------------------------------------------------------
REAL FUNCTION fun(x)
  IMPLICIT NONE
  REAL, INTENT(IN):: x
  fun=4/(1+x**2)
END FUNCTION fun
!---------------------------------------------------------------------------------------------
REAL FUNCTION height(n,a,b)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL, INTENT(IN) :: a,b
  height= (b-a)/REAL(n)
END FUNCTION height
