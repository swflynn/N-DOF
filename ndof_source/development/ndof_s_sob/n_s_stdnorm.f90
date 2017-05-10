!> Computes the inverse cumulative density function (CDF), i.e., the quantile,
! of the standard normal distribution given u uniform on the unit hypercube.
FUNCTION beasley_springer_moro(u) result(x)
    IMPLICIT NONE

    INTEGER :: i, j
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: u

    DOUBLE PRECISION :: r
    DOUBLE PRECISION, DIMENSION(SIZE(u)) :: x, y


    DOUBLE PRECISION, PARAMETER, DIMENSION(0:3) :: a = (/ &
            2.50662823884, &
            -18.61500062529, &
            41.39119773534, &
            -25.44106049637 /)

    DOUBLE PRECISION, PARAMETER, DIMENSION(0:3) :: b = (/ &
            -8.47351093090, &
            23.08336743743, &
            -21.06224101826, &
            3.13082909833 /)

    DOUBLE PRECISION, PARAMETER, DIMENSION(0:8) :: c = (/ &
            0.3374754822726147, &
            0.9761690190917186, &
            0.1607979714918209, &
            0.0276438810333863, &
            0.0038405729373609, &
            0.0003951896511919, &
            0.0000321767881768, &
            0.0000002888167364, &
            0.0000003960315187 /)


    y = u - 0.5D0

    DO j = 1, SIZE(u)
        IF (ABS(y(j)) < 0.42) THEN
            r = y(j)*y(j)
            x(j) = y(j)*(((a(3)*r + a(2))*r + a(1))*r + a(0))/((((b(3)*r + b(2))*r + b(1))*r + b(0))*r + 1)
        ELSE
            IF (y(j) > 0) THEN
                r = LOG(-LOG(1-u(j)))
            ELSE IF (y(j) < 0) THEN
                r = LOG(-LOG(u(j)))
            END IF
            x(j) = c(0) + r*(c(1) + r*(c(2) + r*(c(3) + r*(c(4) + r*(c(5) + r*(c(6) + r*(c(7) + r*c(8))))))))
            IF (y(j) < 0) THEN
                x(j) = -x(j)
            END IF
        END IF
    END DO

END FUNCTION beasley_springer_moro

!==============================================================================!

! commented out because this requires the sobol module, only needed for non-scrambled
!> Returns a d-dimensional Sobol sequence of p points following a standard
!  normal distribution
!subroutine sobol_stdnormal(d, skip, x_stdnormal)
!    use sobol
!    implicit none
    !> dimension
!    INTEGER(kind = 4), INTENT(IN) :: d
 
    !> number of initial points to be skipped
!    INTEGER(kind = 8), INTENT(IN) :: skip   

    !> return an array of doubles, standard normal
!    DOUBLE PRECISION, DIMENSION(d), INTENT(OUT) :: x_stdnormal     

!    interface
!        FUNCTION beasley_springer_moro(u)
!            double precision :: u(:)
!            double precision :: beasley_springer_moro(size(u))
!        end function beasley_springer_moro
!    end interface

!    x_stdnormal = beasley_springer_moro(i8_sobol(int(d, 8), skip))

!END subroutine sobol_stdnormal

!==============================================================================!


!> Returns a d-dimensional Sobol sequence of p points following a standard
!  normal distribution
subroutine scrambled_sobol_stdnormal(d, scrambled_u, x_stdnormal)

    implicit none

    !> dimenstion
    INTEGER(kind = 4), INTENT(IN) :: d

    !> return an array of doubles, standard normal
    DOUBLE PRECISION, DIMENSION(d), INTENT(INOUT) :: x_stdnormal     

    DOUBLE PRECISION, DIMENSION(d), INTENT(IN) :: scrambled_u

    interface
        FUNCTION beasley_springer_moro(u)
            double precision :: u(:)
            double precision :: beasley_springer_moro(size(u))
        end function beasley_springer_moro
    end interface


    x_stdnormal = beasley_springer_moro(scrambled_u)


    END subroutine scrambled_sobol_stdnormal


!==============================================================================!
! Subroutine to calculate permutations for a 9 spatial dimension, 
! variable maximum excitation. 
subroutine vofj(spatial_dim, vmax)

IMPLICIT NONE
Integer, intent(in) :: spatial_dim, vmax
integer :: j,Vm(Spatial_dim),vv(Spatial_dim),k,v1,v2,v3,v4,v5,v6,v7,v8,v9
integer, parameter :: Jmax=220
!INTEGER, v(spatial_dim, Jmax)

j=0
 Vm(1)=Vmax
 if(Spatial_dim>9) stop 'Spatial_dim>9'
  do v1=0,Vm(1)
     vv(1)=v1
     if(Spatial_dim >= 2) then
        Vm(2)=Vm(1)-vv(1)
        do v2=0,Vm(2)
           vv(2)=v2
           if(Spatial_dim >= 3) then
              Vm(3)=Vm(2)-vv(2)
              do v3=0,Vm(3)
                 vv(3)=v3
                 if(Spatial_dim >= 4) then
                    Vm(4)=Vm(3)-vv(3)
                    do v4=0,Vm(4)
                       vv(4)=v4
                       if(Spatial_dim >= 5) then
                          Vm(5)=Vm(4)-vv(4)
                          do v5=0,Vm(5)
                             vv(5)=v5
                             if(Spatial_dim >= 6) then
                                Vm(6)=Vm(5)-vv(5)
                                do v6=0,Vm(6)
                                   vv(6)=v6
                                   if(Spatial_dim >= 7) then
                                      Vm(7)=Vm(6)-vv(6)
                                      do v7=0,Vm(7)
                                         vv(7)=v7
                                         if(Spatial_dim >= 8) then
                                            Vm(8)=Vm(7)-vv(7)
                                            do v8=0,Vm(8)
                                               vv(8)=v8
                                               if(Spatial_dim == 9) then
                                                  Vm(9)=Vm(8)-vv(8)
                                                  do v9=0,Vm(9)
                                                     vv(9)=v9
                                                     j=j+1
!                                                  write(*,*) 'j=',j,vv !here
                                                     do k=1,Spatial_dim
!                                                        v(k,j)=vv(k) !
                                                     enddo
                                                  enddo
                                               else
                                                  j=j+1
                                                  do k=1,Spatial_dim
!                                                    v(k,j)=vv(k) !
                                                  enddo
                                               endif
                                            enddo
                                         else
                                            j=j+1
                                            do k=1,Spatial_dim
!                                               v(k,j)=vv(k) !
                                            enddo
                                         endif
                                      enddo
                                   else
                                      j=j+1
                                      do k=1,Spatial_dim
!                                         v(k,j)=vv(k) !
                                      enddo
                                   endif
                                enddo
                             else
                                j=j+1
                                do k=1,Spatial_dim
!                                   v(k,j)=vv(k)  !
                                enddo
                             endif
                          enddo
                       else 
                          j=j+1 
                          do k=1,Spatial_dim
!                             v(k,j)=vv(k) ! 
                          enddo
                       endif
                    enddo
                 else
                    j=j+1
                    do k=1,Spatial_dim
!                       v(k,j)=vv(k)  ! 
                    enddo
                 endif
              enddo
           else
              j=j+1
              do k=1,Spatial_dim
!                 v(k,j)=vv(k)   !
              enddo
           endif
        enddo
     else
        j=j+1
        do k=1,Spatial_dim
!           v(k,j)=vv(k)    !
        enddo
     endif
  enddo 

 write(*,*) 'Jmax=', j
END subroutine vofj
