! Code to calculate an integral for 2 hermite polnomials and a gaussian
! easiest possible calculation not dimension, no loop over points 
!Generate sobol points (1 dimension) loop over points, evaluate <0|g|0>
!Now calculate for a row of hermite polynomials

PROGRAM simp
  USE sobol
  IMPLICIT NONE

  INTEGER :: d, Nsobol, i, j, k
  INTEGER*8 :: skip
  DOUBLE PRECISION, ALLOCATABLE:: norm(:,:) 
  DOUBLE PRECISION, DIMENSION:: herm(:) 
  DOUBLE PRECISION, DIMENSION(1:10) :: mat


!============ variables for integral ==========================!
  REAL :: gauss, herm1, integ 
!============ variables for integral ==========================!

  d = 1                           
  Nsobol = 10
  skip = 10
  integ = 0.0


  ALLOCATE(norm(d, Nsobol))
!========= generate nsobol array=========================!

  DO i = 1, Nsobol                
    CALL sobol_stdnormal(d,skip,norm(:,i))
!    WRITE(*,*) norm(:,i) 
  END DO
!========= generate nsobol array=========================!

  DO j = 1, Nsobol !Loop each sobol point
  DO k = 1, 10      !Loop across row need to define herm(k) forall k

          herm1 = 1.0 
          gauss = val_op(norm(d,j))
          integ = herm1*gauss*herm(k)
          mat = mat+integ

          WRITE(*,*) 'this is the sobol point', norm(d,j)
          WRITE(*,*) 'this is the hermite', herm
          WRITE(*,*) 'this is the second hermite', herm
          WRITE(*,*) 'this is our gauss', gauss
          WRITE(*,*) 'this is our matrix element', integ

       END DO
  
  END DO
  
  
  WRITE(*,*) 'This is the integral', mat/SIZE(norm(d,:))




!STOP





CONTAINS 

FUNCTION  val_op(x) RESULT (value_fn)
    IMPLICIT NONE

    DOUBLE PRECISION :: x
    DOUBLE PRECISION :: value_fn

    value_fn = EXP(-(X**2))

END FUNCTION val_op


END PROGRAM simp 
