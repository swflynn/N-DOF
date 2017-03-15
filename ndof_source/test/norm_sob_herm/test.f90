! Program uses sobol.f90 module, and the sobol_stdnormal.f90 file
! 'norm_a' is one array containing all of the sobol points (so we can shuffle later on)
! 'herm' takes each entry in 'norm_a' and calculates: H_0 - H_9

! the code can generate any dimension of sobol points 'norm', but 'norm_a' and 'herm' only works for 1 currently
! As is the code is working 3/14/17  -Shane

PROGRAM test
  USE sobol
  IMPLICIT NONE

  INTEGER :: d, Nsobol, i, j, k, l
  INTEGER*8 :: skip
  DOUBLE PRECISION, ALLOCATABLE:: norm(:,:) 
  DOUBLE PRECISION, DIMENSION(1:10) :: herm

  d = 1                           
  Nsobol = 10                   
  skip = 10                        ! skip is not defined at 0

  OPEN(UNIT=10, FILE='data.dat')

  ALLOCATE (norm(d, Nsobol))
!================generate norm_a===============================!
  DO i = 1, Nsobol                

    CALL sobol_stdnormal(d,skip,norm(:,i))
    WRITE(10,*) norm(:,i) 

  END DO
!  PRINT *, norm
!================generate norm_a===============================!

!================generate herm===============================!
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
!================generate herm===============================!

  CLOSE(10)

END PROGRAM test
