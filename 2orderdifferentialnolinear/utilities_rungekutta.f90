MODULE utilities_rungekutta
  IMPLICIT NONE
CONTAINS
  SUBROUTINE eulera(x,x0,v0,dt,fy,fxy)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: dt,v0
    REAL(KIND=8), INTENT(INOUT) :: x,x0
    REAL(KIND=8), INTENT(OUT) :: fy
    REAL(KIND=8), EXTERNAL :: fxy
    
    fy=fxy(v0)
    x=x0+dt*fy
  END SUBROUTINE eulera
  !--------------------------------------------------------------------------------------------
  SUBROUTINE eulerb(v,x0,v0,dt,fy,fxy)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: dt,v0
    REAL(KIND=8), INTENT(INOUT) :: v,x0
    REAL(KIND=8), INTENT(OUT) :: fy
    REAL(KIND=8), EXTERNAL :: fxy

    fy=fxy(x0)
    v=v0+dt*fy
  END SUBROUTINE eulerb

  !--------------------------------------------------------------------------------------------
  SUBROUTINE euler1(x,x0,dt,fy)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), INTENT(INOUT) :: x,x0
    REAL(KIND=8), INTENT(IN) :: fy
    x=x0+dt*fy
  END SUBROUTINE euler1
  !--------------------------------------------------------------------------------------------
  REAL(KIND=8) FUNCTION fvt(x)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: x
    fvt=-SIN(x)
  END FUNCTION fvt
  !--------------------------------------------------------------------------------------------   
  REAL(KIND=8) FUNCTION fxt(v)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: v
    fxt=v
  END FUNCTION fxt
  !--------------------------------------------------------------------------------------------
  SUBROUTINE runge_kutta(x,x0,v,v0,dt)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(INOUT) :: x,x0
    REAL(KIND=8), INTENT(INOUT) :: v,v0
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8) :: f0,f1,f2,f3,fav
    REAL(KIND=8) :: fv0,fv1,fv2,fv3,fvav
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
