PROGRAM int_test
  USE sobol
  IMPLICIT NONE

  INTEGER :: d, Nsobol, n, i, j
  INTEGER :: k, m
  INTEGER*8 :: skip
  DOUBLE PRECISION, ALLOCATABLE:: norm(:,:)
  DOUBLE PRECISION, DIMENSION(1:10) :: herm
  DOUBLE PRECISION, DIMENSION(1:10, 1:10) :: A_km

  d = 1                           
  Nsobol = 100000
  skip = 10                        

  ALLOCATE (norm(d,Nsobol))

OPEN(UNIT=10, FILE='data.dat')
!=========================Get each sobol pointpoint=========================!
  DO n = 1, Nsobol                
    CALL sobol_stdnormal(d,skip,norm(:,n))
  END DO
!=========================Get each sobol pointpoint=========================!

!=========================evaluate each sobol point=========================!
  DO i = 1, Nsobol              
        herm(1) = 1.0             
        herm(2) = 2.0*norm(d,i)       

        
      !======================evaluate each herm for a single sobol point=========================!
        DO j = 3,10       
             herm(j) = (2.0*norm(d,i)*herm(j-1)) - (2.0*(j-2)*herm(j-2))
        END DO
      !=====================evaluate each herm for a single sobol point=========================!
      
      !=====================evaluate each matrix element for a single point=========================!
        DO k = 1, 10
             DO m = 1, 10
                 A_km(k,m) = A_km(k,m) + herm(k)*herm(m)
             END DO
        END DO
      !=====================evaluate each matrix element for a single point=========================!

  END DO
!=========================evaluate each sobol point=========================!
  A_km = A_km / Nsobol

  Write(10,*) A_km
  
  CLOSE(10)

END PROGRAM int_test
