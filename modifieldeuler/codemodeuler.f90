PROGRAM euler                                                                                                                                                                                                      
  USE mdeuler
  IMPLICIT NONE
  REAL :: dx=0.01,y=0.0,y0=0.0
  INTEGER :: iteraciones,i
  iteraciones=INT(1.55/dx)

  OPEN(UNIT=1,FILE="euler",STATUS="OLD")
  DO i=1, iteraciones
     CALL modified(y,y0,dx/2)

     CALL modified(y,y0,dx)
     y0=y
     WRITE(1,*) FLOAT(i)/100, y
  END DO

END PROGRAM euler
