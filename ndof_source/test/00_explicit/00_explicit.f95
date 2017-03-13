PROGRAM basis00

  IMPLICIT NONE
  REAL, PARAMETER :: force_const = 1.0, reduc_mass = 1.0, hbar = 1.0
  REAL, PARAMETER :: upper_limit=10.0
  REAL :: a, coef, increment
  REAL, PARAMETER :: pi = 4*atan(1.0)
  INTEGER, PARAMETER :: step_size=1000
  INTEGER :: i

  a = ((force_const*reduc_mass)**0.5)/hbar
  coef= (1/pi)* (a/2)**0.5

  CALL trapezoid_integration(step_size,upper_limit)

  CONTAINS

    SUBROUTINE trapezoid_integration(step_size,upper_limit)
      IMPLICIT NONE
      INTEGER :: step_size, i
      REAL :: upper_limit, integral, increment, h

      integral = 0.0

      do i=0,step_size
        increment = (upper_limit*i) / step_size

        IF ((i.EQ.0) .OR. (i.EQ.step_size)) THEN
          integral = integral+(2.0*integrand(increment))
        END IF
      ENd DO

      h = upper_limit/step_size
      integral = (h/2.0)*integral*coef

      WRITE (*,*) 'Trapezoid Integration give : ', integral
    END SUBROUTINE trapezoid_integration

    FUNCTION integrand(x) RESULT (values)
      IMPLICIT NONE
      REAL :: x, values

      values = EXP(-(x**2)*(a+0.5)) 
    END FUNCTION integrand

END PROGRAM basis00
