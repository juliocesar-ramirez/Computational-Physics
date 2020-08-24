MODULE utilities
  IMPLICIT NONE
CONTAINS
  !---------------------------------------------------------------------------------------------
  SUBROUTINE initial_values(L,spin)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: L
    INTEGER, INTENT(OUT), DIMENSION(L,L) :: spin
    INTEGER :: i,j
    REAL :: r
    
    DO i=1,L
       DO j=1,L
          CALL RANDOM_NUMBER(r)
          IF(r<0.5) THEN
             spin(i,j)=-1
          ELSE
             spin(i,j)=1
          END IF
       END DO
    END DO
  END SUBROUTINE initial_values
  !---------------------------------------------------------------------------------------------
  SUBROUTINE periodic_boundary_condition(L,i,j,a,b,c,d)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,j,L
    INTEGER, INTENT(OUT) :: a,b,c,d

    IF(i==1) b=L
    IF(i==L) a=1
    IF(j==1) d=L
    IF(j==L) c=1
  END SUBROUTINE PERIODIC_BOUNDARY_CONDITION
  !---------------------------------------------------------------------------------------------
  REAL FUNCTION E1(spin,L,i,j,a,b,c,d)
    INTEGER, INTENT(IN) :: L
    INTEGER, INTENT(IN), DIMENSION(L,L) :: spin
    INTEGER, INTENT(IN) :: i,j,a,b,c,d
    REAL :: J_ising=1.0
    
    E1=J_ising*(FLOAT(spin(i,j)) * ( FLOAT(spin(a,j))+FLOAT(spin(b,j))+FLOAT(spin(i,c))+FLOAT(spin(i,d)) ) )
  END FUNCTION E1
  !---------------------------------------------------------------------------------------------
  SUBROUTINE i_magnetization_and_energy(spin,L,energy,magnatize)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: L
    INTEGER, INTENT(IN), DIMENSION(L,L) :: spin
    REAL, INTENT(OUT) :: energy, magnatize
    
    INTEGER :: a,b,c,d
    INTEGER :: i,j
    REAL :: Esum=0,M=0,N,mag
    DO i=1,L
       DO j=1,L
          
          CALL neighbour(i,j,a,b,c,d)
          
          CALL periodic_boundary_condition(L,i,j,a,b,c,d)
          Esum=Esum+E1(spin,L,i,j,a,b,c,d)
          M=M+FLOAT(spin(i,j ))
       END DO
    END DO
    N=FLOAT(L)**2
    mag=M/N
    Esum=Esum*0.5
    energy=Esum
    magnatize=M
    WRITE(*,*) "Initial energy E, E per spin", Esum,Esum/N
    WRITE(*,*) "Inital magnetization M, M per spin", M, mag 
    
  END SUBROUTINE i_magnetization_and_energy
  !---------------------------------------------------------------------------------------------
  SUBROUTINE neighbour(i,j,a,b,c,d)
    INTEGER, INTENT(IN) :: i,j
    INTEGER, INTENT(OUT) :: a,b,c,d
    a=i+1
    b=i-1
    c=j+1
    d=j-1      

  END SUBROUTINE neighbour
  !---------------------------------------------------------------------------------------------
  SUBROUTINE mt_metropolis(iteration,spin,L,energy, magnatize)
    INTEGER, INTENT(IN) :: iteration
    INTEGER, INTENT(IN) :: L
    REAL, INTENT(IN) :: energy, magnatize
    INTEGER, INTENT(INOUT), DIMENSION(L,L) :: spin
    !    REAL, EXTERNAL :: E
    
    REAL :: r,Ei,Ef,dE,E,M,u,T=2,h,N
    INTEGER :: time, mm , nn, i,j,a,b,c,d

    N= FLOAT(L)**2
    E=energy
    M=magnatize
    DO time=1, iteration

       DO mm=1, L
          DO nn=1,L
             CALL RANDOM_NUMBER(r)
             i=INT(r*FLOAT(L))+1
             CALL RANDOM_NUMBER(r)
             j=INT(r*FLOAT(L))+1
             !             WRITE(*,*) i,j
             CALL neighbour(i,j,a,b,c,d)
             CALL periodic_boundary_condition(L,i,j,a,b,c,d)
             Ei=E1(spin,L,i,j,a,b,c,d)
             spin(i,j)=-spin(j,i)
             Ef=E1(spin,L,i,j,a,b,c,d)
             dE=Ef-Ei
             IF (dE<=0.0) THEN
                E=E+dE
                M=M+(2.0*FLOAT(spin(i,j)))
             ELSE
                u=EXP(-dE/(T))
                CALL RANDOM_NUMBER(h)
                IF(h<u) THEN
                   E=E+dE
                   M=M+(2.0*FLOAT(spin(i,j)))
                ELSE
                   spin(i,j)=-spin(i,j)
                END IF
             END IF

          END DO
       END DO

       WRITE (*,*) time ,M/N, E/N
    END DO

    


  END SUBROUTINE mt_metropolis
END MODULE utilities

