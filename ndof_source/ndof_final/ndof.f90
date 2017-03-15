! Please see the test directory for documentation. 
! Generate normal sobol sequence and evaluate hermite polynomials at each point
! Code is written for 1 dimension only, would need to modify for more. 

PROGRAM main
  USE sobol
  IMPLICIT NONE

  INTEGER :: d, Nsobol, i, j, k
  INTEGER*8 :: skip
  DOUBLE PRECISION, ALLOCATABLE:: norm(:,:)
  DOUBLE PRECISION, DIMENSION(1:10) :: herm

  d = 1                           
  Nsobol = 100                   
  skip = 10                        

  OPEN(UNIT=10, FILE='data.dat')

  ALLOCATE (norm(d,Nsobol))

  DO i = 1, Nsobol                
    CALL sobol_stdnormal(d,skip,norm(:,i))
    WRITE(10,*) norm(:,i)
  END DO

  DO j = 1, Nsobol              
    herm(1) = 1.0             !initialize first 2 polynomials for recursion
    herm(2) = 2.0*norm(d,j)       

    WRITE(10,*), 'next point', norm(d,j)
    WRITE(10,*), herm(1)
    WRITE(10,*), herm(2)

    DO k = 3,10       
      herm(k) = (2.0*norm(d,j)*herm(k-1)) - (2.0*(k-2)*herm(k-2))
      WRITE(10,*), herm(k)
    END DO
  END DO
  CLOSE(10)

END PROGRAM main
