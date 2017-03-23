! Please see the test directory for documentation. 
! Code is written for 1 dimension only, would need to modify for more. 
! Generate norm sobol seq and evaluate each hermite polynomial
! Then calculate the matrix elements for the HO wavefunction
! Optional Analysis calculates convergence as a function of Nsobol

PROGRAM main
  USE sobol
  IMPLICIT NONE
  INTEGER :: d, Nsobol
  INTEGER :: n, i, j, k, m, o, p
  INTEGER, PARAMETER :: deg = 10       !polynomial we want to calculate up to
  INTEGER*8 :: skip
  DOUBLE PRECISION, ALLOCATABLE:: norm(:,:)
  DOUBLE PRECISION, DIMENSION(1:10) :: herm, coef
  DOUBLE PRECISION, DIMENSION(1:10, 1:10) :: A

  d = 1                           
  Nsobol = 100
  skip = 100
  ALLOCATE (norm(d,Nsobol))
  A=0d0

!=============Get each sobol point (normal distribution)====================!
  DO n = 1, Nsobol                
    CALL sobol_stdnormal(d,skip,norm(:,n))
  END DO
  norm=norm/SQRT(2.)      

!===========Recursively Calculate HO Coefficients===========================!
  coef(1) = 1
  coef(2) = 1.0 / (SQRT(2.0))
  DO p = 3,deg
    coef(p) = coef(p-1)*(1 /SQRT(2.0*REAL(p-1)))
  END DO

!=========================evaluate each sobol point=========================!
  OPEN(UNIT=9, FILE='converge.dat')
  DO i = 1, Nsobol              
        herm(1) = 1.0             
        herm(2) = 2.0*norm(1,i)       
        DO j = 3,deg      
             herm(j) =(2.0*norm(d,i)*herm(j-1)) - (2.0*(j-2)*herm(j-2))
        END DO

!=========evaluate each matrix element using the HO wavefunction==============!
        DO k = 1, deg
             DO m = 1, deg
                 A(k,m) = A(k,m) + coef(k)*herm(k)*coef(m)*herm(m)
             END DO
        END DO

!=========Convergence Analysis at n points, remove if not needed==========!
  IF (mod(i,10)==0) THEN
      WRITE(9,*) A / REAL(i)
  END IF

  END DO
  CLOSE(9)
  
  A = A / Nsobol
  OPEN(UNIT=10, FILE='mat_elem.dat')
  do o=1,deg
     Write(10,*) A(1:deg,o)
  enddo
  CLOSE(10)

END PROGRAM main
