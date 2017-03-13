! code to write function for normalized gaussian function to be used in
! <v|G(x)|v>

PROGRAM expec_fn
  IMPLICIT NONE

  INTEGER, PARAMETER :: n = 10
	REAL :: u
  INTEGER :: i

!Check a few values to make sure our function is evaluating correctly
! This is working fine 2/28/17

  DO i = 0, n
    u = 3.0*i/n
    WRITE (*,*) u, e_fn_eval(u)
  END DO

STOP

CONTAINS


   FUNCTION  e_fn_eval(x) RESULT (value_e_fn)
    IMPLICIT NONE

    REAL :: x
    REAL :: value_e_fn

    value_e_fn = EXP(-X**2)

   END FUNCTION e_fn_eval

END PROGRAM expec_fn
