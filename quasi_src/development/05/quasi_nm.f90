! look into lapack for matrix diagonalization
Module quasi_nm
  implicit none
!============================================================================================!
!======================================Global paramaters=====================================!
!============================================================================================!
  double Precision, Parameter :: deg=180/dacos(-1d0)
  double Precision, Parameter :: pi=dacos(-1d0)
  double precision, parameter :: bohr = 0.52917721092
  double precision, parameter :: autocm = 2.194746313D5
  double precision, parameter :: autokcalmol = 627.5096 
  double precision, parameter :: melectron=1822.88839
  double precision, parameter :: Hmass = 1.00782503223*melectron
  double precision, parameter :: Omass = 15.99491461957*melectron
  double precision, allocatable :: sqrt_mass(:),mass(:)
  character(len=2), allocatable :: atom_type(:)
  character(len=5) potential
  INTEGER, PARAMETER :: Vmax(9) = [0,0,0,0,0,0,3,3,3]
!============================================================================================!
!=======================================Global Variables=====================================!
!============================================================================================!
  INTEGER:: Dim, Natoms, Jmax, Vtot
  INTEGER, ALLOCATABLE :: v(:,:)
  double precision, allocatable :: omegak(:)
  double precision :: enum
  double precision, allocatable :: kd_mat(:,:)
!============================================================================================!
CONTAINS
  
  subroutine CM(Natoms,q)
    implicit none
    integer :: Natoms,i
    double precision :: q(3,Natoms),vec(3)
    vec=0
    do i=1,Natoms
       vec(:)=vec(:)+mass(i)*q(:,i)
    enddo
    vec=vec/sum(mass(1:Natoms))
    do i=1,Natoms
       q(:,i)=q(:,i)-vec(:)
    enddo
  end subroutine CM
  
  function Atom_mass(atom)
    implicit none
    double precision :: Atom_mass
    character(len=2), intent(in) :: atom
    if (atom=='O') then
       Atom_mass=Omass
    else if (atom=='H') then
       Atom_mass=Hmass
    else
       write(*,*) 'atom ', atom, ' is not recognized'
       stop 
    endif
  end function Atom_mass

  subroutine frequencies_from_Hess(Dim,Hess,omega,U)
    implicit none
    integer :: i,IERR,Dim
    double precision :: Hess(Dim,Dim),A(Dim,Dim),omega(Dim),FV1(Dim),FV2(Dim),U(Dim,Dim)
    A=Hess
    call RS(Dim,Dim,A,omega,1,U,FV1,FV2,IERR) 
          write(*,*) 'Frequencies from the Hessian:'
          write(*,*) 
    do i=Dim,1,-1
       omega(i)=sign(dsqrt(dabs(omega(i))),omega(i))
       write(*,*) omega(i)*autocm, 'norm= 1?',sum(U(:,i)**2)
    enddo
    
  end subroutine frequencies_from_Hess
  
  subroutine Get_Hessian(q,H)
    implicit none
    integer :: i,j
    double precision ::  H(Dim,Dim),q(Dim),r(Dim),E,force(Dim),force0(Dim)
    double precision, parameter :: s=1d-7
    
    r=q
    call water_potential(Natoms/3, r, E, force0)
    do  i=1,Dim
       r(i)=q(i)+s
       call water_potential(Natoms/3, r, E, force)
       r(i)=q(i)
       do j=1,Dim
          H(i,j)=(force0(j)-force(j))/s
       enddo
    enddo
    ! symmetrize
    do i=1,Dim
       do j=1,i
          if(i.ne.j) H(i,j)=(H(i,j)+H(j,i))/2
          H(i,j)=H(i,j)/(sqrt_mass(i)*sqrt_mass(j)) ! mass-scaled Hessian    \tilde{K} in atomic units
          if(i.ne.j)  H(j,i)=H(i,j)
       enddo
    enddo
   ! write(*,*) 'Here is the mass scaled Hessian', H
    
  end subroutine Get_Hessian
  
  subroutine water_potential(NO,q,energy,force)
    USE iso_c_binding
    USE TIP4P_module
    IMPLICIT NONE
    
    INTEGER, INTENT(IN) :: NO                              ! number of water molecules
    DOUBLE PRECISION, DIMENSION(9*NO), INTENT(IN) :: q   ! coordinates
    DOUBLE PRECISION, DIMENSION(9*NO), INTENT(INOUT) :: force
    DOUBLE PRECISION, INTENT(INOUT) :: energy
    
    if(potential=='tip4p') then
       call TIP4P(NO, q, energy, force)
    else if(potential=='mbpol') then
       call calcpotg(NO, energy, q*bohr, force)
       force = -force*bohr/autokcalmol
       energy = energy/autokcalmol
    else
       stop 'cannot identify the potential'
    endif
    
  end subroutine water_potential

SUBROUTINE permutation(d,Vtot,Vmax)

    Implicit none
    INTEGER :: d, Vtot,Vmax(d)
    INTEGER :: j,vv(d),k,v1,v2,v3,v4,v5,v6,v7,v8,v9

    j=0
    if(d>9) stop 'Spatial_dim>9'
    do v1=0,Vmax(1)
       vv(1)=v1
       if(d>= 2) then
          do v2=0,min(Vtot-vv(1),Vmax(2))
             vv(2)=v2
             if(d>= 3) then
                do v3=0,min(Vtot-sum(vv(1:2)),Vmax(3))
                   vv(3)=v3
                   if(d>= 4) then
                      do v4=0,min(Vtot-sum(vv(1:3)),Vmax(4))
                         vv(4)=v4
                         if(d>= 5) then
                            do v5=0,min(Vtot-sum(vv(1:4)),Vmax(5))
                               vv(5)=v5
                               if(d>= 6) then
                                  do v6=0,min(Vtot-sum(vv(1:5)),Vmax(6))
                                     vv(6)=v6
                                     if(d>= 7) then
                                        do v7=0,min(Vtot-sum(vv(1:6)),Vmax(7))
                                           vv(7)=v7
                                           if(d>= 8) then
                                              do v8=0,min(Vtot-sum(vv(1:7)),Vmax(8))
                                                 vv(8)=v8
                                                 if(d== 9) then
                                                    do v9=0,min(Vtot-sum(vv(1:8)),Vmax(9))
                                                       vv(9)=v9
                                                       j=j+1
                                                    enddo
                                                 else
                                                    j=j+1
                                                 endif
                                              enddo
                                           else
                                              j=j+1
                                           endif
                                        enddo
                                     else
                                        j=j+1
                                     endif
                                  enddo
                               else
                                  j=j+1
                               endif
                            enddo
                         else 
                            j=j+1 
                         endif
                      enddo
                   else
                      j=j+1
                   endif
                enddo
             else
                j=j+1
             endif
          enddo
       else
          j=j+1
       endif
    enddo
    Jmax = j
!===================================With Jmax, Run again=================================!
    !WRITE(*,*) 'Jmax = ', Jmax
    ALLOCATE(v(d,Jmax))
    
    j=0
    if(d>9) stop 'Spatial_dim>9'
    do v1=0,Vmax(1)
       vv(1)=v1
       if(d>= 2) then
          do v2=0,min(Vtot-vv(1),Vmax(2))
             vv(2)=v2
             if(d>= 3) then
                do v3=0,min(Vtot-sum(vv(1:2)),Vmax(3))
                   vv(3)=v3
                   if(d>= 4) then
                      do v4=0,min(Vtot-sum(vv(1:3)),Vmax(4))
                         vv(4)=v4
                         if(d>= 5) then
                            do v5=0,min(Vtot-sum(vv(1:4)),Vmax(5))
                               vv(5)=v5
                               if(d>= 6) then
                                  do v6=0,min(Vtot-sum(vv(1:5)),Vmax(6))
                                     vv(6)=v6
                                     if(d>= 7) then
                                        do v7=0,min(Vtot-sum(vv(1:6)),Vmax(7))
                                           vv(7)=v7
                                           if(d>= 8) then
                                              do v8=0,min(Vtot-sum(vv(1:7)),Vmax(8))
                                                 vv(8)=v8
                                                 if(d== 9) then
                                                    do v9=0,min(Vtot-sum(vv(1:8)),Vmax(9))
                                                       vv(9)=v9
                                                       j=j+1
                                                       v(:,j)=vv 
                                                    enddo
                                                 else
                                                    j=j+1
                                                    v(:,j)=vv           
                                                 endif
                                              enddo
                                           else
                                              j=j+1
                                              v(:,j)=vv           
                                           endif
                                        enddo
                                     else
                                        j=j+1
                                        v(:,j)=vv           
                                     endif
                                  enddo
                               else
                                  j=j+1
                                  v(:,j)=vv                            
                               endif
                            enddo
                         else 
                            j=j+1 
                            v(:,j)=vv
                         endif
                      enddo
                   else
                      j=j+1
                      v(:,j)=vv           
                   endif
                enddo
             else
                j=j+1
                v(:,j)=vv           
             endif
          enddo
       else
          j=j+1
          v(:,j)=vv           
       endif
    enddo
!    write(*,*) v
  END SUBROUTINE permutation

subroutine kd_matrix(kd_mat)
  implicit none
  integer :: i, i1, qnum
  double precision :: kd_mat(Jmax,Jmax)
!  double precision :: omegak(dim)
!  integer :: Jmax

  qnum=0
  enum=0d0
    DO i=1,Jmax
      DO i1=1,Jmax
        IF(i.NE.i1) THEN
          kd_mat(i,i1) = 0.
        ELSE
            write(*,*) 'here we go getting kd_mat'
            write(*,*) 'we pass in v, omegak'
            write(*,*) v(:,i)
            write(*,*) omegak
            call kd_ene(v(:,i), omegak)
            write(*,*) 'here is enum'
            write(*,*) enum
            kd_mat(i,i1) = enum
            write(*,*) 'are these the same?'
            write(*,*) kd_mat(i,i1)
        END IF
      END DO
    END DO
  write(*,*) 'kd_mat'
  write(*,*) kd_mat
end subroutine kd_matrix

subroutine kd_ene(v, omegak)
  implicit none
  double precision :: omegak(dim)
  integer :: i, k, i1, v(dim), qnum

  qnum = 0
  enum = 0d0

  do i = 1,dim
    qnum = qnum + v(i)
  end do
  write(*,*) 'Here is qnum'
  write(*,*) qnum
!  enum=0d0
  do k=1,qnum
    enum = enum + omegak(k)*(k+.5)
    write(*,*) 'here is enum'
    write(*,*) enum
  end do

end subroutine kd_ene

END MODULE quasi_nm
!===========================================================================================!
!===========================================================================================!
!===========================================================================================!
!===========================================================================================!
PROGRAM main
  USE quasi_nm
  
  IMPLICIT NONE
  real :: initial_time, final_time  
  integer(kind=8) :: skip
  integer :: Nsobol, i,j
!  double precision, allocatable :: omega(:), omegak(:), Hess(:,:), U(:,:), H1(:,:), q(:),q0(:)
  double precision, allocatable :: omega(:), Hess(:,:), U(:,:), H1(:,:), q(:),q0(:)
!  double precision, allocatable :: Hess(:,:), U(:,:), H1(:,:), q(:),q0(:)
  double precision, allocatable :: q1(:), force(:)
  double precision :: freq_cutoff, E, Ener, V0, potdif, harpot
  character(len=50) coord_in,frequencies_out
!===========================================================================================!
!====================================Variables My Code======================================! 
! There is a U defined in the Get_Hessian function U(dim,dim)
! I will define my code U => Umat, contains the PE matrix elements Umat(Jmax,Jmax)
! Replace scrambled_z with y for normal mode points. 
!===========================================================================================!
  DOUBLE PRECISION, ALLOCATABLE :: y(:), herm(:,:), A(:,:,:), Umat(:,:), U1mat(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: C(:,:), FV1(:), FV2(:), eigenvalues(:)
  INTEGER:: m, n, IERR, j1, k
  DOUBLE PRECISION :: B
!ittm  external calc_pot_link2f90
!ittm  external calc_pot_link2f90_g
!==========================================================================================!
!===============================Variables to run ssobol.f==================================!
!==========================================================================================!
INTEGER :: TAUS, IFLAG, max_seq, SAM
LOGICAL :: FLAG(2)
!==========================================================================================!
!===============================Variables to run ssobol.f==================================!
!==========================================================================================!
!  allocatable :: kd_mat(:,:)
!  double precision :: enum
  integer :: qnum, i1
!===========================================================================================!
!====================================Read Input File========================================! 
!===========================================================================================! 
  CALL CPU_TIME(initial_time)
  OPEN(60,FILE='input.dat')
  READ(60,*) Vtot
  READ(60,*) potential                
  READ(60,*) coord_in
  READ(60,*) freq_cutoff
  READ(60,*) frequencies_out
  READ(60,*) Nsobol
  READ(60,*) skip
  CLOSE(60)
  freq_cutoff = freq_cutoff /autocm
!==========================================================================================!
!===============================Variables to run ssobol.f==================================!
  SAM = 1
  max_seq = 30
  IFLAG = 1
!===========================================================================================!
!======================================Read xyz File========================================! 
!===========================================================================================!
  OPEN(61,File=coord_in)
  READ(61,*) Natoms
  READ(61,*)
  Dim= 3*Natoms ! cartesian coordinates
!===========================================================================================!
!======================================(vlad allocations)===================================!
!===========================================================================================! 
  ALLOCATE(omega(Dim), omegak(Dim), atom_type(Natoms), sqrt_mass(Dim), mass(Natoms), &
       q(Dim), q0(Dim), force(Dim), Hess(Dim,Dim),  H1(Dim,Dim), &
       U(Dim,Dim))
!===========================================================================================!
!============================Read Atom Type, get mass=======================================! 
!======================q0 contains x,y,z coordinates of equ config==========================!
!===========================================================================================!
  DO i=1,Natoms
     READ(61,*) atom_type(i), q0(3*i-2:3*i)   ! coordinates in Angstroms
     mass(i)=Atom_mass(atom_type(i))
     sqrt_mass(3*i-2:3*i)=SQRT(mass(i))
  END DO
 CLOSE(61) 
 q0=q0/bohr            ! convert coordinates to atomic units
!===========================================================================================!
!=========================Let user see what they defined====================================! 
!=========================For convenience, can be removed===================================!
  WRITE(*,*) 'Natoms=              ',Natoms
  WRITE(*,*) 'Dim=                 ',Dim
  WRITE(*,*) 'potential=           ',potential
  WRITE(*,*) 'Maximum Excitation=  ',Vtot
  WRITE(*,*) 'frequency cutoff=    ',freq_cutoff*autocm
  WRITE(*,*) 'Nsobol=              ',Nsobol
  WRITE(*,*) 'skip=                ',skip
  WRITE(*,*) 'molecule:            ', atom_type
!===========================================================================================!
!====================================Equ Config=============================================! 
!===========================================================================================!
  CALL water_potential(Natoms/3, q0, E, force)
  WRITE(*,*) 'E0 =',E*autokcalmol, 'kcal/mol'
  WRITE(*,*) 'force ='
  DO i=1,Natoms
     WRITE(*,*) force(3*i-2:3*i)
  END DO
!===========================================================================================!
!====================================Hessian Frequencies====================================! 
!===========================================================================================!
  CALL Get_Hessian(q0,Hess)          
  OPEN(62,file=frequencies_out)
  CALL frequencies_from_Hess(Dim,Hess,omega,U)
! call  quasi_MC(q,omega,U,Nsobol,skip,freq_cutoff)
!===========================================================================================!
!==================================Normal Modes U(dim,dim)==================================! 
!==============================sqrt(.5) comes from std normal dist==========================!
!===========================================================================================!
!  write(*,*) 'This is omega', omega*autocm    ! testing remove after
  DO k=1,Dim
     omegak(k) = ABS(omega(k))
     IF(omegak(k)<freq_cutoff) THEN
        omegak(k) = freq_cutoff
     !   write(*,*) omegak*autocm   ! testing 
     !   write(*,*) 'Use cutoff'    
     ELSE
        omegak(k) = omegak(k) 
     !   WRITE(*,*) omegak*autocm
     !   write(*,*) 'Use the original frequency'
     END IF
     DO i = 1, Dim
        U(i,k) = SQRT(0.5/omegak(k)) / sqrt_mass(i)*U(i,k)
!        write(*,*) 'Here are the normal modes', U(i,k)
     END DO 
   END DO ! loop over eigenvalues
!        write(*,*) 'Here are the normal modes'
!        write(*,*)  U
!===========================================================================================!
!============================Determine Jmax and Permutations================================!
!===========================================================================================!
  CALL permutation(dim,Vtot,Vmax)
  ALLOCATE(Umat(Jmax,Jmax), U1mat(Jmax,Jmax), C(Jmax,Jmax), FV1(Jmax), FV2(Jmax))
  ALLOCATE(eigenvalues(Jmax)) 
  ALLOCATE(kd_mat(Jmax,Jmax))
  ALLOCATE(y(dim), herm(0:Vtot,dim), A(0:Vtot,0:Vtot,dim))
!===========================================================================================!
!=================================Kronicker-Delta Matrix====================================!
!===========================================================================================!
!    kd_mat = 0d0 ????
!  call kd_matrix(kd_mat)  ?????
  write(*,*) 'did i make the kd matrix?'
!===========================================================================================!
  Umat=0d0  ! initialize PE matrix
  U1mat=0d0 
  OPEN(UNIT=80, FILE='matrix.dat')
  OPEN(UNIT=81, FILE='eigenvalues.dat')
!===========================================================================================!
!================================Loop over Sobol Points=====================================!
!===============================Normalize with Gaussian=====================================!
!=============Calculate (our normalization) Hermite's, up to Vtot, ini poly=0 ==============!
!===========================================================================================!
  
  CALL INSSOBL(FLAG,dim,Nsobol,TAUS,y,max_seq,IFLAG)
  DO i = 1, Nsobol
     CALL GOSSOBL(y)
     CALL sobol_stdnormal(dim, y)
     y = y/SQRT(2.)    ! factor from definition of normal distribution
     herm(0,:) = 1.0                       ! Re-normalized hermite polynomial now
     herm(1,:) = SQRT(2.)*y(:)       
     DO j = 2,Vtot      
        herm(j,:) = (SQRT(2./j)*y(:)*herm(j-1,:)) - (SQRT(((j-1d0)/j))*herm(j-2,:))
     END DO
!===================================Evaluate Herm * Herm ===================================!
!==========================Matrix A: herm(deg)*herm(deg), deg(0,Vtot)=======================!
!===========================================================================================!
     DO m=0,Vtot
        DO n=0,Vtot
           A(m,n,:) = herm(m,:)*herm(n,:)
        END DO
     END DO
!===========================================================================================!
!==========================Difference Potential=============================================!
!===========================================================================================!
     q1=q
     DO j = 1, Dim
        q1(:) = q1(:) + y(j)*U(:,j)
     END DO
!     write(*,*) 'Here is q1'
!     write(*,*) q1
     CALL water_potential(Natoms/3, q1, potdif, force)
     do j = 1,dim
        potdif = potdif - 0.5*omegak(j)*omegak(j)*y(j)*y(j) 
     END DO
!     write(*,*) 'Here is the  Actual Potential: ', V0
!     write(*,*) 'Here is the  Harmonic Potential: ', harpot
!     write(*,*) 'Here is the potential difference: ', potdif
!================================Evaluate Matrix Elements ==================================!
!===============Matrix Umat: PE matrix elements. U1mat for partial average==================!
!======================V0 has the difference potential======================================!
     DO j=1,Jmax         ! Loop over jmax
        DO j1=j,Jmax 
           B=potdif
           DO k=1,dim
              B=B*A(v(k,j),v(k,j1),k)
           END DO
           U1mat(j1,j) = U1mat(j1,j) + B
        END DO
     END DO
! Write out partial average and flush
! Add pe diff matrix to kdelta matrix 
     IF(MOD(i,1)==0) THEN
        Umat = Umat + U1mat
        U1mat = 0d0
        C=Umat/i
        do j=1,Jmax
!           C(j,j)=C(j,j)+ kd_ener
        enddo
        WRITE(80,*) i, C !Writes entire matrix out
        CALL RS(Jmax,Jmax,C,eigenvalues,1,B,FV1,FV2,IERR) 
        WRITE(81,*) i, ABS(1-eigenvalues(:))
        FLUSH(80)
        FLUSH(81)
      END IF

  END DO ! end loop over sobol points
  CLOSE(UNIT=80)

  CLOSE(UNIT=81)

 CALL CPU_TIME(final_time)
 WRITE(*,*) 'TOTAL TIME: ', final_time - initial_time
 CLOSE(62)
!===========================================================================================!
!========================================output.dat=========================================! 
!===========================================================================================!
 OPEN(90,FILE='output.dat')
 WRITE(90,*) 'Here is the output file for your calculation'
 WRITE(90,*) 'Natoms=                     ', Natoms
 WRITE(90,*) 'Dim=                        ', Dim
 WRITE(90,*) 'Vtot =                      ', Vtot
 WRITE(90,*) 'Jmax =                      ', Jmax
 WRITE(90,*) 'potential=                  ', potential
 WRITE(90,*) 'Nsobol=                     ', Nsobol
 WRITE(90,*) 'skip=                       ', skip
 WRITE(90,*) 'coord_in=                   ', coord_in
 WRITE(90,*) 'freq_out=                   ', frequencies_out
 WRITE(90,*) 'Calculation total time (s)  ', final_time - initial_time
 CLOSE(90)

END PROGRAM main
