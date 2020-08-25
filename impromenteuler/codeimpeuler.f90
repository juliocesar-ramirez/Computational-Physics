PROGRAM euler
  USE mdeuler
  IMPLICIT NONE
  REAL :: dx=0.01,y=0.0,y0=0.0,f0,fe,fxyav
  INTEGER :: iteraciones,i
  iteraciones=INT(1.55/dx)

  OPEN(UNIT=1,FILE="euler",STATUS="OLD")
  DO i=1, iteraciones
     f0=fxy(y)
     CALL modified(y,y0,dx)
     fe=fxy(y)
     fxyav=(f0+fe)/2
     
     CALL improved(y,y0,dx,fxyav)
     y0=y
     WRITE(1,*) FLOAT(i)/100, y
  END DO
    
END PROGRAM euler
