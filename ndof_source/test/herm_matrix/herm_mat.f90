PROGRAM mat_eval
  USE sobol
  IMPLICIT NONE

  INTEGER :: d, Nsobol, n, i, j
  INTEGER :: k, m, o
  INTEGER, PARAMETER :: deg = 3       ! set which polynomial we want to calculate to (3-10)
  INTEGER*8 :: skip
  DOUBLE PRECISION, ALLOCATABLE:: norm(:,:) !matrix of normal sobol points
  DOUBLE PRECISION, DIMENSION(1:10) :: herm !matrix evalue harmite polynomials
  DOUBLE PRECISION, DIMENSION(1:10, 1:10) :: A  !matrix elements <H_k|H_m>

  d = 1                           
  Nsobol = 100000
  skip = 10                        

  ALLOCATE (norm(d,Nsobol))

OPEN(UNIT=10, FILE='data.dat')
A=0d0
!=========================Get each sobol pointpoint=========================!
  DO n = 1, Nsobol                
    CALL sobol_stdnormal(d,skip,norm(:,n))
  END DO
norm=norm/sqrt(2.)
!=========================Get each sobol pointpoint=========================!

!=========================evaluate each sobol point=========================!
  DO i = 1, Nsobol              
        herm(1) = 1.0             
        herm(2) = 2.0*norm(1,i)       

        
      !==================evaluate each polynomial for a single sobol point=======================!
        DO j = 3,deg      
             herm(j) = (2.0*norm(d,i)*herm(j-1)) - (2.0*(j-2)*herm(j-2))
        END DO
      !=====================evaluate each polynomial for a single sobol point=========================!
      
      !=====================evaluate each matrix element for a single point=========================!
        DO k = 1, deg
             DO m = 1, deg
                 A(k,m) = A(k,m) + herm(k)*herm(m)
             END DO
        END DO
      !=====================evaluate each matrix element for a single point=========================!

  END DO
!=========================evaluate each sobol point=========================!
  A = A / Nsobol            !normalize integrals

  
!=========================write out matrix elements in nice format=========================!
  do o=1,deg
     Write(10,*) A(1:deg,o)
  enddo
!=========================write out matrix elements in nice format=========================!
  
  CLOSE(10)

END PROGRAM mat_eval
