! Evaluate H_9 for a given x

PROGRAM herm_fn_9
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 10
	REAL :: u
  INTEGER :: i

!Check a few values to make sure our function is evaluating correctly
! Does not work numerical instability 3/2/17

  DO i = 0, n
    u = 2.0*i/n
    WRITE (*,*) u, herm_fn_9_eval(u)
  END DO

STOP

CONTAINS


   FUNCTION  herm_fn_9_eval(x) RESULT (value_herm_fn_9)
    IMPLICIT NONE

    REAL :: x
    REAL :: value_herm_fn_9

    value_herm_fn_9 = 512.0*(x**9) - 9216.0*(x**7) + 48384.0*(x**5) - 80640.0*(x**3) &
      & + 30240.0*(x)

   END FUNCTION herm_fn_9_eval

END PROGRAM herm_fn_9
