! This is my first attempt at computing something (2/20/17)
! This program intends to calculate a simplified form of the <0|G(x)|0>
! Harmonic oscialator basis functions. I have already included the
! hermite polynomial and simplified the total expression for this
! evaluation. 
! I am evaluating integral with trapezoid rule:
! See http://f90in15minutes.wikidot.com/numerical-integration


PROGRAM basis00

  IMPLICIT NONE
  REAL, PARAMETER :: force_const = 1.0, reduc_mass = 1.0, hbar = 1.0
  REAL, PARAMETER :: upper_limit=3.0
  REAL :: a, coef, u
  REAL, PARAMETER :: pi = 4*atan(1.0)
  INTEGER, PARAMETER :: step_size=10
  INTEGER :: i

  a = ((force_const*reduc_mass)**0.5)/hbar
  coef= (1/pi)* (a/2)**0.5

  !  WRITE (*,*) 'The alpha paramater is ', a 
!    WRITE (*,*) 'The constant of integration is', coef

! evaluate our function at #step_size points, to make sure things are working
! I may need to insert a cutoff for f(x) = 0 at an upper limit of 10 get
!f(x) = 10^-60

  DO i = 0,step_size
     u = upper_limit*i/step_size
     WRITE (*,*) u, integrand(u)
  END DO

  STOP

  CONTAINS

! function to calculation x, f(x) values 
    FUNCTION integrand(x) RESULT (values)
      IMPLICIT NONE
      REAL :: x, values

      values = EXP(-(x**2)*(a+0.5)) 
    END FUNCTION integrand

END PROGRAM basis00
