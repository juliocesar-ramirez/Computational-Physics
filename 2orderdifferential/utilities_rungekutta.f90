MODULE utilities_rungekutta
  IMPLICIT NONE
CONTAINS
  SUBROUTINE eulera(x,x0,v0,dt,fy,fxy)
    IMPLICIT NONE
    REAL, INTENT(IN) :: dt,v0
    REAL, INTENT(INOUT) :: x,x0
    REAL, INTENT(OUT) :: fy
    REAL, EXTERNAL :: fxy
    
    fy=fxy(v0)
    x=x0+dt*fy
  END SUBROUTINE eulera
  !--------------------------------------------------------------------------------------------
  SUBROUTINE eulerb(v,x0,v0,dt,fy,fxy)
    IMPLICIT NONE
    REAL, INTENT(IN) :: dt,v0
    REAL, INTENT(INOUT) :: v,x0
    REAL, INTENT(OUT) :: fy
    REAL, EXTERNAL :: fxy

    fy=fxy(x0)
    v=v0+dt*fy
  END SUBROUTINE eulerb

  !--------------------------------------------------------------------------------------------
  SUBROUTINE euler1(x,x0,dt,fy)
    IMPLICIT NONE
    REAL, INTENT(IN) :: dt
    REAL, INTENT(INOUT) :: x,x0
    REAL, INTENT(IN) :: fy
    x=x0+dt*fy
  END SUBROUTINE euler1
  !--------------------------------------------------------------------------------------------
  REAL FUNCTION fvt(x)
    IMPLICIT NONE
    REAL, INTENT(IN) :: x
    REAL, PARAMETER :: k=1,m=1
    fvt=-k/m*x
  END FUNCTION fvt
  !--------------------------------------------------------------------------------------------   
  REAL FUNCTION fxt(v)
    IMPLICIT NONE
    REAL, INTENT(IN) :: v
    fxt=v
  END FUNCTION fxt
  !--------------------------------------------------------------------------------------------
  SUBROUTINE runge_kutta(x,x0,v,v0,dt)
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: x,x0
    REAL, INTENT(INOUT) :: v,v0
    REAL, INTENT(IN) :: dt
    REAL :: f0,f1,f2,f3,fav
    REAL :: fv0,fv1,fv2,fv3,fvav
    CALL eulera(x,x0,v0,dt/2,f0,fxt)
    CALL eulerb(v,x0,v0,dt/2,fv0,fvt)

    CALL eulera(x,x0,v0,dt/2,f1,fxt)
    CALL eulerb(v,x0,v0,dt/2,fv1,fvt)


    CALL eulera(x,x0,v0,dt/2,f2,fxt)
    CALL eulerb(v,x0,v0,dt/2,fv2,fvt)


    f3=fxt(v)
    fv3=fvt(x)

    fav=(f0+2*f1+2*f2+f3)
    fvav=(fv0+2*fv1+2*fv2+fv3)

    CALL euler1(x,x0,dt/6,fav)
    CALL euler1(v,v0,dt/6,fvav)
  END SUBROUTINE runge_kutta
END MODULE utilities_rungekutta
