PROGRAM rk4
  USE utilities_rungekutta
  IMPLICIT NONE
  REAL :: dx=0.02,x=0,x0=0,v=0,v0=0.1
  INTEGER :: iteraciones,i

  iteraciones=INT(20/dx)
  OPEN(UNIT=1,FILE="rk4",STATUS="OLD")
  DO i=1, iteraciones

     CALL runge_kutta(x,x0,v,v0,dx)
     WRITE(1,*) FLOAT(i)/100, x
     
     x0=x
     v0=v
  END DO
    
END PROGRAM rk4
