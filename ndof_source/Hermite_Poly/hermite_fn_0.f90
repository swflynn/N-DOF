! develop code to evaluate any hermite polynomial for a value of x. 
! The first function is simply 1

PROGRAM herm_fn_0
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 10
	REAL :: u
  INTEGER :: i

!Check a few values to make sure our function is evaluating correctly
! This is working fine 2/28/17 all answers are 1

  DO i = 0, n
    u = 3.0*i/n
    WRITE (*,*) u, herm_fn_eval(u)
  END DO

STOP

CONTAINS


   FUNCTION  herm_fn_eval(x) RESULT (value_herm_fn)
    IMPLICIT NONE

    REAL :: x
    REAL :: value_herm_fn

    value_herm_fn = 1.0*(x**0)

   END FUNCTION herm_fn_eval

END PROGRAM herm_fn_0 
