! Please see the test directory for documentation. 
! Code is written for 1 dimension only, would need to modify for more. 
! Generate norm sobol seq and evaluate each hermite polynomial
! Then calculate the matrix elements for the HO wavefunction

PROGRAM main
  USE sobol
  IMPLICIT NONE

  INTEGER :: d, Nsobol
  INTEGER :: i, j, k, m, n, o, p      
  INTEGER, PARAMETER :: deg = 10       !polynomial we want to calculate to (3-10)
  INTEGER*8 :: skip
  DOUBLE PRECISION, ALLOCATABLE:: norm(:,:) !matrix of normal sobol points
  DOUBLE PRECISION, DIMENSION(1:10) :: herm, coef 
  DOUBLE PRECISION, DIMENSION(1:10, 1:10) :: A  !matrix elements <H_k|e^-x^2/2|H_m>

  d = 1                           
  Nsobol = 10000000
  skip = 1000

  ALLOCATE (norm(d,Nsobol))

  OPEN(UNIT=10, FILE='data.dat')
  A=0d0                             

  DO n = 1, Nsobol                
    CALL sobol_stdnormal(d,skip,norm(:,n))
  END DO
  norm=norm/sqrt(2.)

  coef(1) = 1
  coef(2) = 1.0 / (SQRT(2.0))
  DO p = 3,deg
     coef(p) = coef(p-1)*(1 /SQRT(2.0*REAL(p-1)))
  END DO

  DO i = 1, Nsobol              
        herm(1) = 1.0             
        herm(2) = 2.0*norm(1,i)       
        DO j = 3,deg      
             herm(j) = (2.0*norm(d,i)*herm(j-1)) - (2.0*(j-2)*herm(j-2))
        END DO
      
        DO k = 1, deg
             DO m = 1, deg
                 A(k,m) = A(k,m) + coef(k)*herm(k)*coef(m)*herm(m)
             END DO
        END DO
  END DO
  A = A / Nsobol            
  
  DO o=1,deg
     WRITE(10,*) A(1:deg,o)
     END DO
  
  CLOSE(10)

END PROGRAM main
