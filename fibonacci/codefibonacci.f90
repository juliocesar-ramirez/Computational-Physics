PROGRAM main
  IMPLICIT NONE
  REAL, EXTERNAL :: sequence
  INTEGER,PARAMETER :: nterminos =20 
  INTEGER :: count
  REAL, DIMENSION(nterminos),SAVE :: terminos
  terminos(1)=1
  terminos(2)=1

  DO count=1,nterminos,1
     CALL avanzar(sequence,count,terminos,terminos(count),terminos(count+1))
END DO

CALL imprimir
WRITE(*,*) terminos
END PROGRAM main
!------------------------------------------------------------

REAL FUNCTION sequence (a,b)
  IMPLICIT NONE
  REAL, INTENT (IN) :: a,b
  sequence=a+b
  RETURN
END FUNCTION sequence
!------------------------------------------------------------
REAL FUNCTION ratio (a,b)
  IMPLICIT NONE
  REAL, INTENT(IN) :: a,b

  ratio=a/b
END FUNCTION ratio
!  ----------------------------------------------------------
SUBROUTINE imprimir
  WRITE(*,100)
  WRITE(*,100)
  WRITE(*,110) "FIBONNACI"
  WRITE(*,100)

100 FORMAT ("======================================================================")  
110 FORMAT (A40)  

END SUBROUTINE imprimir
!------------------------------------------------------------
SUBROUTINE avanzar(sequence,count,terminos,a,b)
  REAL, EXTERNAL :: sequence
  INTEGER, INTENT(IN) :: count
  REAL, INTENT(INOUT), DIMENSION(nterminos) :: terminos

  terminos(count+2)=sequence(a,b)
END SUBROUTINE avanzar
!------------------------------------------------------------
