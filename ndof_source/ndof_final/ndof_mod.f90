MODULE ndof_mod




  CONTAINS 

!=====================================================================!

   SUBROUTINE Hermite_Coeff(degree,A,B)  
      IMPLICIT NONE
      INTEGER i,j,degree
      REAL*8  A(0:10), B(0:10,0:10)

      B(0,0)=1.D0 ; B(1,0)=0.D0 ; B(1,1)=2.D0

      IF (degree == 0) THEN
        A(0) = 1
        RETURN 

      ELSE IF (degree == 1) THEN
        A(0) = 0
        A(1) = 2
        RETURN 

      ELSE IF (degree>1) THEN 
        DO i = 2, degree
          B(i,0)=-2.D0*(i-1)*B(i-2,0)
          DO j = 1, i
           B(i,j)=2.D0*B(i-1,j-1)-2.D0*(i-1)*B(i-2,j) !Recursion Relation
          END dO 
        END DO 
        DO i = 0, degree
          A(i)=B(degree,i)
        END DO
      END IF

      RETURN 

    END SUBROUTINE Hermite_Coeff 


!=====================================================================!


END MODULE ndof_mod
