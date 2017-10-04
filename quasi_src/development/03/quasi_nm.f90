! This is a stripped down version of quasi_nm.f90
! use this for testing, I have removed all the subroutines I did not use to get PE matrix
! Take this code and copy it over to scratch when you need it

Module quasi_nm
  implicit none
  double Precision, Parameter :: deg=180/dacos(-1d0)
  double Precision, Parameter :: pi=dacos(-1d0)
  double precision, parameter :: hbar = 1d0
  double precision, parameter :: bohr = 0.52917721092
  double precision, parameter :: autoeV = 27.211385
  double precision, parameter :: eVtoau = 3.674932379D-2
  double precision, parameter :: autocm = 2.194746313D5
  double precision, parameter :: autokcalmol = 627.5096 
  double precision, parameter :: Ktoau = 3.1668114D-6
  double precision, parameter :: cmtoau = 4.5563352527D-6
  double precision, parameter :: melectron=1822.88839
  double precision, parameter :: Hmass = 1.00782503223*melectron
  double precision, parameter :: Omass = 15.99491461957*melectron
  double precision, allocatable :: sqrt_mass(:),mass(:)
  character(len=2), allocatable :: atom_type(:)
  character(len=5) potential
  integer :: Dim, Natoms, Lmon  
!==============================Global paramaters from old code===============================!
  INTEGER, ALLOCATABLE :: v(:,:)
  INTEGER :: Jmax
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
       write(*,*) omega(i)*autocm
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
  
  SUBROUTINE permutation(d,Vmax)
    Implicit none
    INTEGER :: d, Vmax
    INTEGER :: j,Vm(d),vv(d),k,v1,v2,v3,v4,v5,v6,v7,v8,v9

    j=0
       Vm(1)=Vmax
       if(d>9) stop 'Spatial_dim>9'
        do v1=0,Vm(1)
           vv(1)=v1
           if(d>= 2) then
              Vm(2)=Vm(1)-vv(1)
              do v2=0,Vm(2)
                 vv(2)=v2
                 if(d>= 3) then
                    Vm(3)=Vm(2)-vv(2)
                    do v3=0,Vm(3)
                       vv(3)=v3
                       if(d>= 4) then
                          Vm(4)=Vm(3)-vv(3)
                          do v4=0,Vm(4)
                             vv(4)=v4
                             if(d>= 5) then
                                Vm(5)=Vm(4)-vv(4)
                                do v5=0,Vm(5)
                                   vv(5)=v5
                                   if(d>= 6) then
                                      Vm(6)=Vm(5)-vv(5)
                                      do v6=0,Vm(6)
                                         vv(6)=v6
                                         if(d>= 7) then
                                            Vm(7)=Vm(6)-vv(6)
                                            do v7=0,Vm(7)
                                               vv(7)=v7
                                               if(d>= 8) then
                                                  Vm(8)=Vm(7)-vv(7)
                                                  do v8=0,Vm(8)
                                                     vv(8)=v8
                                                     if(d== 9) then
                                                        Vm(9)=Vm(8)-vv(8)
                                                        do v9=0,Vm(9)
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
!================With Jmax, Run again to determine v(d,Jmax)===========================!
!WRITE(*,*) 'Jmax = ', Jmax
ALLOCATE(v(d,Jmax))

j=0
 Vm(1)=Vmax
 if(d>9) stop 'Spatial_dim>9'
  do v1=0,Vm(1)
     vv(1)=v1
     if(d>= 2) then
        Vm(2)=Vm(1)-vv(1)
        do v2=0,Vm(2)
           vv(2)=v2
           if(d>= 3) then
              Vm(3)=Vm(2)-vv(2)
              do v3=0,Vm(3)
                 vv(3)=v3
                 if(d>= 4) then
                    Vm(4)=Vm(3)-vv(3)
                    do v4=0,Vm(4)
                       vv(4)=v4
                       if(d>= 5) then
                          Vm(5)=Vm(4)-vv(4)
                          do v5=0,Vm(5)
                             vv(5)=v5
                             if(d>= 6) then
                                Vm(6)=Vm(5)-vv(5)
                                do v6=0,Vm(6)
                                   vv(6)=v6
                                   if(d>= 7) then
                                      Vm(7)=Vm(6)-vv(6)
                                      do v7=0,Vm(7)
                                         vv(7)=v7
                                         if(d>= 8) then
                                            Vm(8)=Vm(7)-vv(7)
                                            do v8=0,Vm(8)
                                               vv(8)=v8
                                               if(d== 9) then
                                                  Vm(9)=Vm(8)-vv(8)
                                                  do v9=0,Vm(9)
                                                     vv(9)=v9
                                                     j=j+1
                                                     do k=1,d
                                                        v(k,j)=vv(k)       
                                                     enddo
                                                  enddo
                                               else
                                                  j=j+1
                                                  do k=1,d
                                                    v(k,j)=vv(k)           
                                                  enddo
                                               endif
                                            enddo
                                         else
                                            j=j+1
                                            do k=1,d
                                               v(k,j)=vv(k)                
                                            enddo
                                         endif
                                      enddo
                                   else
                                      j=j+1
                                      do k=1,d
                                         v(k,j)=vv(k)                      
                                      enddo
                                   endif
                                enddo
                             else
                                j=j+1
                                do k=1,d
                                   v(k,j)=vv(k)                            
                                enddo
                             endif
                          enddo
                       else 
                          j=j+1 
                          do k=1,d
                             v(k,j)=vv(k)                                  
                          enddo
                       endif
                    enddo
                 else
                    j=j+1
                    do k=1,d
                       v(k,j)=vv(k)                                       
                    enddo
                 endif
              enddo
           else
              j=j+1
              do k=1,d
                 v(k,j)=vv(k)                                              
              enddo
           endif
        enddo
     else
        j=j+1
        do k=1,d
           v(k,j)=vv(k)                                                    
        enddo
     endif
  enddo 
!WRITE(*,*) v! (This contains all of the permutations given Vmax, dim). 
END SUBROUTINE permutation

END MODULE quasi_nm
!===========================================================================================!
!===========================================================================================!
!===========================================================================================!
PROGRAM main
  USE quasi_nm
  
  IMPLICIT NONE
  real :: initial_time, final_time  
  integer(kind=8) :: skip
  integer :: Nsobol, i,j, Vmax
  double precision, allocatable :: omega(:), omegak(:), Hess(:,:), U(:,:), H1(:,:), q(:),q0(:)
  double precision, allocatable :: q1(:), force(:)
  double precision :: freq_cutoff, E, Ener, V0, potdif, harpot
  character(len=50) coord_in,frequencies_out
!===========================================================================================!
  DOUBLE PRECISION, ALLOCATABLE :: y(:), herm(:,:), A(:,:,:), Umat(:,:), U1mat(:,:), kron(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: C(:,:), FV1(:), FV2(:), eigenvalues(:)
  INTEGER:: m, n, IERR, j1, k
  DOUBLE PRECISION :: B
!ittm  external calc_pot_link2f90
!ittm  external calc_pot_link2f90_g
!===============================Variables to run ssobol.f==================================!
INTEGER :: TAUS, IFLAG, max_seq, SAM
LOGICAL :: FLAG(2)
!====================================Read Input File========================================! 
  CALL CPU_TIME(initial_time)
  OPEN(60,FILE='input.dat')
  READ(60,*) Vmax
  READ(60,*) potential                
  READ(60,*) coord_in
  READ(60,*) freq_cutoff
  READ(60,*) frequencies_out
  READ(60,*) Nsobol
  READ(60,*) skip
  CLOSE(60)
  freq_cutoff = freq_cutoff /autocm
!===============================Variables to run ssobol.f==================================!
  SAM = 1
  max_seq = 30
  IFLAG = 1
!======================================Read xyz File========================================! 
  OPEN(61,File=coord_in)
  READ(61,*) Natoms
  READ(61,*)
  Dim= 3*Natoms ! cartesian coordinates
!======================================(vlad allocations)===================================!
  ALLOCATE(omega(Dim), omegak(Dim), atom_type(Natoms), sqrt_mass(Dim), mass(Natoms), &
       q(Dim), q0(Dim), force(Dim), Hess(Dim,Dim),  H1(Dim,Dim), &
       U(Dim,Dim))
!============================Read Atom Type, get mass=======================================! 
!======================q0 contains x,y,z coordinates of equ config==========================!
  DO i=1,Natoms
     READ(61,*) atom_type(i), q0(3*i-2:3*i)   ! coordinates in Angstroms
     mass(i)=Atom_mass(atom_type(i))
     sqrt_mass(3*i-2:3*i)=SQRT(mass(i))
  END DO
 CLOSE(61) 
 q0=q0/bohr            ! convert coordinates to atomic units
!=========================Let user see what they defined====================================! 
!=========================For convenience, can be removed===================================!
  WRITE(*,*) 'Natoms=              ',Natoms
  WRITE(*,*) 'Dim=                 ',Dim
  WRITE(*,*) 'potential=           ',potential
  WRITE(*,*) 'Maximum Excitation=  ',Vmax
  WRITE(*,*) 'frequency cutoff=    ',freq_cutoff*autocm
  WRITE(*,*) 'Nsobol=              ',Nsobol
  WRITE(*,*) 'skip=                ',skip
  WRITE(*,*) 'molecule:            ', atom_type
!====================================Equ Config=============================================! 
  CALL water_potential(Natoms/3, q0, E, force)
  WRITE(*,*) 'E0 =',E*autokcalmol, 'kcal/mol'
  WRITE(*,*) 'force ='
  DO i=1,Natoms
     WRITE(*,*) force(3*i-2:3*i)
  END DO
!====================================Hessian Frequencies====================================! 
  CALL Get_Hessian(q0,Hess)          
  OPEN(62,file=frequencies_out)
  CALL frequencies_from_Hess(Dim,Hess,omega,U)
! call  quasi_MC(q,omega,U,Nsobol,skip,freq_cutoff)
!==================================Normal Modes U(dim,dim)==================================! 
!==============================sqrt(.5) comes from std normal dist==========================!
!  write(*,*) 'This is omega', omega*autocm    ! testing remove after
!  write(*,*)
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
        write(*,*) 'Here are the normal modes'
        write(*,*)  U
!==========================Start code for PE matrix elements================================!
  CALL permutation(dim,Vmax)
  ALLOCATE(Umat(Jmax,Jmax), U1mat(Jmax,Jmax), C(Jmax,Jmax), FV1(Jmax), FV2(Jmax))
  ALLOCATE(kron(Jmax,Jmax))
  ALLOCATE(eigenvalues(Jmax)) 
  ALLOCATE(y(dim), herm(0:Vmax,dim), A(0:Vmax,0:Vmax,dim))
!===========================================================================================!
  Umat=0d0  ! initialize PE matrix
  U1mat=0d0 

  OPEN(UNIT=80, FILE='matrix.dat')
  OPEN(UNIT=81, FILE='eigenvalues.dat')
!==============================Compute Kronicker delta matrix===============================!
kron = 0d0
write(*,*) 'here is kron', kron
  do j = 1, Jmax
   do j1 = 1,Jmax
   if(j.ne.j1) Then
     kron(j,j1) = 0d0
   else
     kron(j,j1) = 5d0
   end if
   end do
  end do
write(*,*) 'Here it is again'
write(*,*) kron
!================================Loop over Sobol Points=====================================!
!===============================Normalize with Gaussian=====================================!
!=============Calculate (our normalization) Hermite's, up to Vmax, ini poly=0 ==============!
  V0 = 0d0
  potdif = 0d0
  harpot = 0d0
  CALL INSSOBL(FLAG,dim,Nsobol,TAUS,y,max_seq,IFLAG)
  DO i = 1, Nsobol
     potdif = 0d0   ! potential difference
     V0 = 0d0        ! actual potential
     harpot = 0d0   ! harmonic potential
     CALL GOSSOBL(y)
     CALL sobol_stdnormal(dim, y)
     y = y/SQRT(2.)    ! factor from definition of normal distribution
     herm(0,:) = 1.0                       ! Re-normalized hermite polynomial now
     herm(1,:) = SQRT(2.)*y(:)       
     DO j = 2,Vmax      
        herm(j,:) = (SQRT(2./j)*y(:)*herm(j-1,:)) - (SQRT(((j-1d0)/j))*herm(j-2,:))
     END DO
!===================================Evaluate Herm * Herm ===================================!
!==========================Matrix A: herm(deg)*herm(deg), deg(0,Vmax)=======================!
     DO m=0,Vmax
        DO n=0,Vmax
           A(m,n,:) = herm(m,:)*herm(n,:)
        END DO
     END DO
!==========================Difference Potential=============================================!
     q1=q
     DO j = 1, Dim
        q1(:) = q1(:) + y(j)*U(:,j)
     END DO
!     write(*,*) 'Here is q1'
!     write(*,*) q1
     CALL water_potential(Natoms/3, q1, Ener, force)
     V0 = V0 + Ener
     do j = 1,dim
        harpot = harpot + 0.5*omegak(j)*omegak(j)*y(j)*y(j) 
     END DO
     potdif = V0 - harpot
!     write(*,*) 'Here is the  Actual Potential: ', V0
!     write(*,*) 'Here is the  Harmonic Potential: ', harpot
!     write(*,*) 'Here is the potential difference: ', potdif
! testing out harmonic terms add in this code later
!================================Evaluate Matrix Elements ==================================!
!===============Matrix Umat: PE matrix elements. U1mat for partial average==================!
!======================V0 has the difference potential======================================!
     DO j=1,Jmax         ! Loop over jmax
        DO j1=j,Jmax 
           B=A(v(1,j),v(1,j1),1)
           DO k=2,dim
              B=B*A(v(k,j),v(k,j1),k)
           END DO
           B = B* potdif
           Umat(j1,j) = Umat(j1,j) + B
        END DO
     END DO
! Write out partial average and flush
     IF(MOD(i,1)==0) THEN
        Umat = Umat + U1mat
        U1mat = 0d0
        C=Umat/i
        WRITE(80,*) i, C !Writes entire matrix out
        CALL RS(Jmax,Jmax,C,eigenvalues,0,B,FV1,FV2,IERR) 
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
!========================================output.dat=========================================! 
 OPEN(90,FILE='output.dat')
 WRITE(90,*) 'Here is the output file for your calculation'
 WRITE(90,*) 'Natoms=                     ', Natoms
 WRITE(90,*) 'Dim=                        ', Dim
 WRITE(90,*) 'Vmax =                      ', Vmax
 WRITE(90,*) 'Jmax =                      ', Jmax
 WRITE(90,*) 'potential=                  ', potential
 WRITE(90,*) 'Nsobol=                     ', Nsobol
 WRITE(90,*) 'skip=                       ', skip
 WRITE(90,*) 'coord_in=                   ', coord_in
 WRITE(90,*) 'freq_out=                   ', frequencies_out
 WRITE(90,*) 'Calculation total time (s)  ', final_time - initial_time
 CLOSE(90)

END PROGRAM main
