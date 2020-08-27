PROGRAM rk4
  USE utilities_rungekutta
  IMPLICIT NONE
  REAL(KIND=8) :: dx=0.001D0,x=0.0D0,x0=0.0D0,v=0.0D0,v0=1.996D0
  INTEGER :: iteraciones,i

  iteraciones=INT(60/dx)
  OPEN(UNIT=1,FILE="x",STATUS="OLD")
  OPEN(UNIT=2,FILE="v",STATUS="OLD")
  DO i=1, iteraciones

     CALL runge_kutta(x,x0,v,v0,dx)
     WRITE(1,*) FLOAT(i)/100, x
     WRITE(2,*) FLOAT(i)/100, v
     
     x0=x
     v0=v
  END DO
    
END PROGRAM rk4
