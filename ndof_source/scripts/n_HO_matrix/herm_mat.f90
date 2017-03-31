!Code to calculate HO matrix elements using normal sobol sequence with convergence analysis (for a single dimension only). 
 
!User should set what degree (deg) polynomial they would like to calculate up to
!Then set the Nsobol number of points to use and the skip
!Code generates a sequence of normally distributed sobol points 'norm' valid for 1 dimension (using sobol.f90 and sobol_stdnormal.f90 functions+subroutines). 
!Then calculates using recursion the HO wavefunction coeficients
!Then the hermite polynomials are calculated recursively
!Finally the matrix elements are calculated by the product of the hermite polynomials and coeficients. 

!Everything is working as of 3-20-17   -Shane 
PROGRAM mat_eval
  
  USE sobol
  IMPLICIT NONE
  INTEGER :: d, Nsobol                  ! dimension and number of points
  INTEGER :: n, i, j, k, m, o, p         !iterative variables
  INTEGER, PARAMETER :: deg = 10       !polynomial we want to calculate up to
  INTEGER*8 :: skip                    !set similar to Nsobol
  DOUBLE PRECISION, ALLOCATABLE:: norm(:,:)   !vector of normal sobol points
  DOUBLE PRECISION, DIMENSION(1:10) :: herm, coef    !hermite polynomial and coef vectors
  DOUBLE PRECISION, DIMENSION(1:10, 1:10) :: A     ! matrix elements

  d = 1                           
  Nsobol = 100
  skip = 100
  ALLOCATE (norm(d,Nsobol))
  A=0d0
!=========================Get each sobol pointpoint=========================!
  DO n = 1, Nsobol                
    CALL sobol_stdnormal(d,skip,norm(:,n))
  END DO
  norm=norm/SQRT(2.)      
!=========================Get each sobol pointpoint=========================!

!=========================Get each wavefn coef.=========================!
  coef(1) = 1
  coef(2) = 1.0 / (SQRT(2.0))
  DO p = 3,deg
    coef(p) = coef(p-1)*(1 /SQRT(2.0*REAL(p-1)))
  END DO
  !=========================Get each wavefn coef.=========================!
  OPEN(UNIT=9, FILE='converge.dat')
!=========================evaluate each sobol point=========================!
DO i = 1, Nsobol              
        herm(1) = 1.0             
        herm(2) = 2.0*norm(1,i)       
!=============evaluate each herm for a single sobol point===================!
        DO j = 3,deg      
             herm(j) =(2.0*norm(d,i)*herm(j-1)) - (2.0*(j-2)*herm(j-2))
        END DO
!=========evaluate each herm for a single sobol point======================!

!==============evaluate each matrix element for a single point===================!
        DO k = 1, deg
             DO m = 1, deg
                 A(k,m) = A(k,m) + coef(k)*herm(k)*coef(m)*herm(m)
             END DO
        END DO
 !=============evaluate each matrix element for a single point======================!

! add in function to print results as a function of N
  IF (mod(i,10)==0) THEN
      WRITE(9,*) (A) / REAL(i)
  END IF

  END DO

  CLOSE(9)
!=========================evaluate each sobol point=========================!
  A = A / Nsobol
!=========================write out matrix elements=========================!
  OPEN(UNIT=10, FILE='final_matrix.dat')
  do o=1,deg
     Write(10,*) A(1:deg,o)
  enddo
  CLOSE(10)
!=========================write out matrix elements=========================!

END PROGRAM mat_eval
