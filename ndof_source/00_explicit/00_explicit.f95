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

  DO i = 0,step_size
     u = upper_limit*i/step_size
     WRITE (*,*) u, integrand(u)
  END DO

  STOP

  CONTAINS

    FUNCTION integrand(x) RESULT (values)
      IMPLICIT NONE
      REAL :: x, values

      values = EXP(-(x**2)*(a+0.5)) 
    END FUNCTION integrand

END PROGRAM basis00
