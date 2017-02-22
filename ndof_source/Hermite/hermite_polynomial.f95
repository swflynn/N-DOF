! Need to generate the hermite polynomials for the harmonic oscillator wavefunction
! Have found code to generate their coeficient, we then know they are jsut powers of x based on how many integers there are (prints 0 if no power)

PROGRAM Hermite

REAL*8 :: A(0:10)
REAL*8 :: B(0:10,0:10)
INTEGER :: degree,k

   PRINT*,' '
   degree = 6
     WRITE(*,25) degree

     CALL Hermite_Coeff(degree,A,B)

     DO k = 0, degree
       WRITE(*,50) k, A(k)  
     END DO 
     PRINT *,' '
   print *,' '

   stop
25 format(' Hermite polynomial coefficients for Hermite Polynomial ',i2/)
50 format('   A(',i2,') = ',f10.0)

END PROGRAM HERMITE

SUBROUTINE Hermite_Coeff(degree,A,B)  
  INTEGER i,j,degree
  REAL*8  A(0:10), B(0:10,0:10)
  !Establish l0 and l1 coefficients
  B(0,0)=1.d0 ; B(1,0)=0.d0 ; B(1,1)=2.d0
  !Return if order is less than two
  IF (degree>1) THEN 
    DO i = 2, degree
      B(i,0)=-2.d0*(i-1)*B(i-2,0)
      do j = 1, i
        !Basic recursion relation
        B(i,j)=2.d0*B(i-1,j-1)-2.d0*(i-1)*B(i-2,j)
      end do
    end do
    do i = 0, degree
      A(i)=B(degree,i)
    end do
  end if
  return
END SUBROUTINE Hermite_Coeff
