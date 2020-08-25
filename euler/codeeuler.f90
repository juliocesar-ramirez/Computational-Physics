PROGRAM euler
  IMPLICIT NONE
  REAL :: dx=0.01,fxy,y=0.0
  INTEGER :: iteraciones,i
  iteraciones=INT(1.55/dx)

  OPEN(UNIT=1,FILE="euler",STATUS="OLD")
  DO i=1, iteraciones


     fxy=1.0+y**2
     y=y+dx*fxy

     WRITE(1,*) REAL(i)/100,y
  END DO
     
END PROGRAM euler
