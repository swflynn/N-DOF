! Code to calculate an integral for 2 hermite polnomials and a gaussian

PROGRAM int_test

  USE sobol
  IMPLICIT NONE

  INTEGER :: d, Nsobol, i, j, k
  INTEGER*8 :: skip
  DOUBLE PRECISION, ALLOCATABLE:: norm(:,:)
  DOUBLE PRECISION, DIMENSION(1:10) :: herm

! need to calculate matrix 
  !  DOUBLE PRECISION, DIMENSION(1:10, 1:10) :: mat

  d = 1                           
  Nsobol = 1                   
  skip = 10                        

  ALLOCATE (norm(d,Nsobol))





  DO i = 1, Nsobol                
    CALL sobol_stdnormal(d,skip,norm(:,i))
    WRITE(*,*) norm(:,i)
  END DO






  DO j = 1, Nsobol              


    herm(1) = 1.0             
    herm(2) = 2.0*norm(d,j)       

    WRITE(*,*), 'next point', norm(d,j)
    WRITE(*,*), herm(1)
    WRITE(*,*), herm(2)






    DO k = 3,10       
      herm(k) = (2.0*norm(d,j)*herm(k-1)) - (2.0*(k-2)*herm(k-2))
      WRITE(*,*), herm(k)








    END DO





  END DO







!CLOSE(10)

END PROGRAM int_test
