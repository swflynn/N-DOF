! develop code to evaluate any hermite polynomial for a value of x. 

PROGRAM herm_fn_8
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 10
	REAL :: u
  INTEGER :: i

!Check a few values to make sure our function is evaluating correctly
! Does not work numerical instability

  DO i = 0, n
    u = 2.0*i/n
    WRITE (*,*) u, herm_fn_8_eval(u)
  END DO

STOP

CONTAINS


   FUNCTION  herm_fn_8_eval(x) RESULT (value_herm_fn_8)
    IMPLICIT NONE

    REAL :: x
    REAL :: value_herm_fn_8

    value_herm_fn_8 = 256.0*(x**8) - 3584.0*(x**6) + 13440.0*(x**4) - 13440.0*(x**2) &
      & + 1680.0

   END FUNCTION herm_fn_8_eval

END PROGRAM herm_fn_8
