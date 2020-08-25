MODULE mdeuler
CONTAINS
  SUBROUTINE modified(y,y0,dx)
    REAL, INTENT(IN) :: dx
    REAL, INTENT(INOUT) :: y,y0
    y=y0+dx*fxy(y)
  END SUBROUTINE modified
  !--------------------------------------------------------------------------------------------                                                                                                                    
  REAL FUNCTION fxy(y)
    REAL, INTENT(IN) :: y
    fxy=1+y**2
  END FUNCTION fxy
END MODULE mdeuler


