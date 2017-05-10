PROGRAM s_mat_eval

 IMPLICIT NONE
 INTEGER :: Nsobol
 INTEGER, PARAMETER :: deg = 4                  !highest degree polynomial
 INTEGER, PARAMETER :: d = 3                    !spatial dimension
 INTEGER, PARAMETER :: vmax = 9                 !Highest Excitation to consider
 INTEGER, PARAMETER :: jmax = 220               ! Permutations from subroutine
 INTEGER :: i, j, k, m
 REAL :: initial_time, final_time
 DOUBLE PRECISION, ALLOCATABLE :: scrambled_u(:,:), scrambled_z(:), herm(:,:), coef(:,:), A(:,:,:)
! Integer, Allocatable :: v(:,:)


  CALL CPU_TIME(initial_time)
 
 Nsobol=1
 !Nsobol=100000
 ALLOCATE(scrambled_u(d, Nsobol), scrambled_z(d), herm(deg, d), coef(deg, d), A(deg,deg, d))
!ALLOCATE(v(d,jmax))


 A = 0d0
 
 !=========================Read in Scrambled Sequence=========================!
 OPEN(UNIT=70, FILE='s_sobol_unif.dat', STATUS='OLD', ACTION='READ')
 READ(70,*) scrambled_u
 CLOSE(UNIT=70)
 !=========================Read in Scrambled Sequence=========================!

 open(unit=90, file='testcoef.dat')
!============ Generates a coef for each spatial dim (deg,s_dim)=============!
  coef(1,:) = 1.0
  coef(2,:) = 1.0 / SQRT(2.)
  DO i = 3,deg
    coef(i,:) = coef(i-1,:) * (1 / SQRT(2.*(i-1)))
  END DO
!============ Generates a coef for each spatial dim (deg,s_dim)=============!


!=========================evaluate each sobol point each dimension=========================!
  DO i = 1, Nsobol              
    CALL scrambled_sobol_stdnormal(d, scrambled_u(:,i), scrambled_z(:))
         scrambled_z = scrambled_z/SQRT(2.)

        herm(1,:) = 1.0             
        herm(2,:) = 2.0*scrambled_z(:)       
        DO j = 3,deg      
             herm(j,:) =(2.*scrambled_z(:)*herm(j-1,:)) - (2.*(j-2)*herm(j-2,:))
        END DO
        herm(:,:)=herm(:,:)*coef(:,:)
!=========================evaluate each sobol point each dimension=========================!


!=============evaluate product herm(deg)*herm(deg) for all deg, spatial===========!
!for each dimension take herm(i,:)*herm(j,:) and loop i and j from 1 to deg
        Do k=1,deg
            DO m=1,deg
            A(k,m,:) = herm(k,:)*herm(m,:)
            END DO
        END DO
!=============evaluate product herm(deg)*herm(deg) for all deg, spatial===========!
! I am here 4/18/17
  
  END DO

!  write(90,*) A
  close(90) 


!make sure new subroutine works
! CALL vofj(d,vmax, v)

 
  DEALLOCATE(scrambled_u, scrambled_z, herm, coef, A)
!  DEALLOCATE(v)
   CALL CPU_TIME(final_time)
   WRITE(*,*) 'TOTAL TIME: ', final_time - initial_time

END PROGRAM s_mat_eval
