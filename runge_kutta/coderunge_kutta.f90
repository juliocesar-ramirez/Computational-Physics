PROGRAM rk4
  USE utilities_rungekutta
  IMPLICIT NONE
  REAL :: dx=0.01,y=0.0,y0=0.0
  INTEGER :: iteraciones,i

  iteraciones=INT(1.55/dx)
  OPEN(UNIT=1,FILE="rk4",STATUS="OLD")
  DO i=1, iteraciones

     CALL runge_kutta(y,y0,dx)
     WRITE(1,*) FLOAT(i)/100, y
     y0=y
  END DO
    
END PROGRAM rk4
