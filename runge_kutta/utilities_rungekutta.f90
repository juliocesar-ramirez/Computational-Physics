MODULE utilities_rungekutta
  IMPLICIT NONE
CONTAINS
  SUBROUTINE euler(y,y0,dx,fy)
    IMPLICIT NONE
    REAL, INTENT(IN) :: dx
    REAL, INTENT(INOUT) :: y,y0
    REAL, INTENT(OUT) :: fy
    fy=fxy(y)
    y=y0+dx*fy
  END SUBROUTINE euler
  !--------------------------------------------------------------------------------------------
  SUBROUTINE euler1(y,y0,dx,fy)
    IMPLICIT NONE
    REAL, INTENT(IN) :: dx
    REAL, INTENT(INOUT) :: y,y0
    REAL, INTENT(IN) :: fy
   y=y0+dx*fy
 END SUBROUTINE euler1
 !--------------------------------------------------------------------------------------------
 REAL FUNCTION fxy(y)
    IMPLICIT NONE
    REAL, INTENT(IN) :: y
    fxy=1+y**2
  END FUNCTION fxy
  !--------------------------------------------------------------------------------------------
  SUBROUTINE runge_kutta(y,y0,dx)
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: y,y0
    REAL, INTENT(IN) :: dx
    REAL :: f0,f1,f2,f3,fav

    CALL euler(y,y0,dx/2,f0)
    CALL euler(y,y0,dx/2,f1)
    CALL euler(y,y0,dx,f2)
    f3=fxy(y)
    fav=(f0+2*f1+2*f2+f3)
    CALL euler1(y,y0,dx/6,fav)
  END SUBROUTINE runge_kutta
END MODULE utilities_rungekutta
