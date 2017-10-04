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
  double precision, parameter :: Dmass = 2.01410177812*melectron
  double precision, parameter :: Cmass = 12.0000000*melectron
  double precision, parameter :: Omass = 15.99491461957*melectron
  double precision, parameter :: Fmass = 18.99840316273*melectron
  double precision, parameter :: Clmass = 34.968852682*melectron
  double precision, parameter :: Brmass = 78.9183376*melectron
  double precision, parameter :: Imass = 126.9044719*melectron
  
  double precision, allocatable :: sqrt_mass(:),mass(:)
  character(len=2), allocatable :: atom_type(:)
  character(len=5) potential
  integer :: Dim, Natoms, Lmon, Vtot, Vmax(9)
!============================================================================================!
!==============================Global paramaters from old code===============================!
!============================================================================================!
  INTEGER, ALLOCATABLE :: v(:,:)
  INTEGER :: Jmax
!============================================================================================!
!============================================================================================!
!============================================================================================!
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
    else if (atom=='C') then
       Atom_mass=Cmass
    else if (atom=='H') then
       Atom_mass=Hmass
    else if (atom=='D') then
       Atom_mass=Dmass
    else  if (atom=='Cl') then
       Atom_mass=Clmass
    else  if (atom=='Br') then
       Atom_mass=Brmass
    else  if (atom=='I') then
       Atom_mass=Imass
    else  if (atom=='F') then
       Atom_mass=Fmass
    else
       write(*,*) 'atom ', atom, ' is not recognized'
       stop 
    endif
  end function Atom_mass
!========================================================================================!
! I will want ot put my code here, for now I am going to write it in the actual program
!========================================================================================!
!  subroutine quasi_MC(q,omega,U,Nsobol,skip,freq_cutoff)

!    implicit none
!    integer :: i,k,j,Nsobol
!    integer(kind=8) :: skip
!    double precision ::  q(Dim),q1(Dim),G(Dim,Dim),G1(Dim,Dim),U(Dim,Dim),y(Dim)
!   double precision ::  r(Dim),omegak,d(Dim),omega(Dim),force(Dim),Ener,freq_cutoff,V0  ! added V0 
! removed index in omega, since we do not repeat calculation like in SCP
       
!ittm    external calc_pot_link2f90
!ittm    external calc_pot_link2f90_g

! issue with enddo here 

!    do k=1,Dim      ! just run from 1 to Dim
!       omegak=dabs(omega(k))
!       if(omegak<freq_cutoff) omegak=freq_cutoff
!       do i=1,Dim
!          U(i,k)=sqrt(0.5/omegak)/sqrt_mass(i)*U(i,k)   ! the sqrt(0.5) factor comes from using the standard normal distr exp(-y^2/2)  
          ! to sample points with the distribution exp(-y^2)
!       enddo
!     enddo ! add en enddo that was missing here

!    V0=0d0
    !     compute the averages
!    do k=1,Nsobol
!       call sobol_stdnormal(Dim,skip,y)
!       q1=q
!       do i=1,Dim  ! sum over Dim largest eigenvalues, change index from k to i inside this loop
!          q1(:)=q1(:)+y(i)*U(:,i)
!       enddo
!       call water_potential(Natoms/3, q1, Ener, force)
!       V0=V0+Ener
!    enddo
!    V0=V0/Nsobol
!WRITE(*,*) 'V0 from quasi MC', V0   ! Testing
       
!  end subroutine quasi_MC

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
!       if(dabs(omega(i))>1d-8) write(*,*) omega(i)*autocm
       write(*,*) omega(i)*autocm
    enddo
    
  end subroutine frequencies_from_Hess
  
  subroutine Get_Hessian(q,H)

    implicit none
    integer :: i,j
    double precision ::  H(Dim,Dim),q(Dim),r(Dim),E,force(Dim),force0(Dim)
    double precision, parameter :: s=1d-7
    
!ittm    external calc_pot_link2f90
!ittm    external calc_pot_link2f90_g
    
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
  
  subroutine project(L,Hess,q,Dproj)
    implicit none
    integer :: i,j,k,Dproj,L
    double precision ::  Hess(L,L),vec(L,Dproj),&
         Projector(L,L),q(L),B(L,L)
    ! 3 translational vectors
    vec=0d0
    do k=1,3
       do i=1,L/3
          vec(k+(i-1)*3,k)=sqrt_mass(k+(i-1)*3)
       enddo
    enddo
    ! 3 rotational vectors
    if(Dproj==6) then
       do i=1,L/3
          vec(2+(i-1)*3,4)=q(3+(i-1)*3)*sqrt_mass(2+(i-1)*3)
          vec(3+(i-1)*3,4)=-q(2+(i-1)*3)*sqrt_mass(3+(i-1)*3)
          vec(3+(i-1)*3,5)=q(1+(i-1)*3)*sqrt_mass(3+(i-1)*3)
          vec(1+(i-1)*3,5)=-q(3+(i-1)*3)*sqrt_mass(1+(i-1)*3)
          vec(1+(i-1)*3,6)=q(2+(i-1)*3)*sqrt_mass(1+(i-1)*3)
          vec(2+(i-1)*3,6)=-q(1+(i-1)*3)*sqrt_mass(2+(i-1)*3)
       enddo
    endif
    !     Gram-Schmidt 
    call gram_schmidt(L,Dproj,vec)
    do i=1,L
       Projector(i,i)=1
       do j=1,L
          if(i.ne.j) Projector(i,j)=0d0
       enddo
    enddo
    
    do k=1,Dproj
       do i=1,L
          do j=1,L
             Projector(i,j)=Projector(i,j)-vec(i,k)*vec(j,k)
          enddo
       enddo
    enddo
    call rmatmat(L,Hess,Projector,B)
    call rmatmat(L,Projector,B,Hess)
    
  end subroutine project
   
  subroutine rmatvec(N,A,b,c,flag)
    implicit none
    logical flag
    integer :: N,i,j
    double precision ::  A(N,N),b(N),c(N)
    c=0d0
    if(flag) then   ! transpose A
       do j=1,N
          do i=1,N
             c(i)=c(i)+A(i,j)*b(j)
          enddo
       enddo
    else
       do i=1,N
          c(i)=dot_product(A(1:N,i),b)
       enddo
    endif
    
  end subroutine rmatvec

  subroutine rmatmat(N,A,B,C)
    !
    !     careful: C(j,i)=sum_k A(k,j)*B(k,i)
    !
    implicit none
    integer :: N,i,j
    double precision ::  A(N,N),B(N,N),C(N,N),AT(N)
    do i=1,N
       do j=1,N
          AT(j)=A(i,j)
       enddo
       do j=1,N
          C(i,j)=dot_product(AT,B(:,j))
          !            C(j,i)=rdot(N,A(1,j),B(1,i))   ! careful: B is transposed
       enddo
    enddo
    
  end subroutine rmatmat

  subroutine rmatmat1(N,A,B,C)
    !
    !     careful: C(j,i)=sum_k A(k,j)*B(k,i)
    !
    implicit none
    integer :: N,i,j
    double precision ::  A(N,N),B(N,N),C(N,N)!,BT(N)
    do i=1,N
       do j=1,N
          !            BT(j)=B(i,j)
          !         enddo
          !            C(j,i)=rdot(N,A(1,j),BT)
          C(j,i)=dot_product(A(:,j),B(:,i))   ! careful: B is transposed
       enddo
    enddo
    
  end subroutine rmatmat1
  
  subroutine gram_schmidt(N,k,vec)
    implicit none
    integer :: k,N,l,l1
    double precision ::  vec(N,k)
    do l=1,k
       do l1=1,l-1
          vec(:,l)=vec(:,l)-dot_product(vec(:,l),vec(:,l1))*vec(:,l1)
       enddo
       vec(:,l)=vec(:,l)/dsqrt(dot_product(vec(:,l),vec(:,l)))
    enddo    
  end subroutine gram_schmidt
  
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
    !================With Jmax, Run again to determine v(d,Jmax)===========================!
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
    !WRITE(*,*) v! (This contains all of the permutations given Vtot, dim). 
  END SUBROUTINE permutation

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
  integer :: Nsobol, i,j, Vtot
  double precision, allocatable :: omega(:), Hess(:,:), U(:,:), H1(:,:), q(:),q0(:),force(:)
  double precision :: freq_cutoff, E
  character(len=50) coord_in,frequencies_out
!===========================================================================================!
!====================================Variables My Code======================================! 
! There is a U defined in the Get_Hessian function U(dim,dim)
! I will define my code U => Umat, contains the PE matrix elements Umat(Jmax,Jmax)
!===========================================================================================!
  DOUBLE PRECISION, ALLOCATABLE :: scrambled_z(:), herm(:,:), A(:,:,:), Umat(:,:), U1(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: C(:,:), FV1(:), FV2(:), eigenvalues(:)
  INTEGER:: m, n, IERR, j1, k
  DOUBLE PRECISION :: B, omegak
!ittm  external calc_pot_link2f90
!ittm  external calc_pot_link2f90_g
!===========================================================================================!
!====================================Read Input File========================================! 
!===========================================================================================! 
  CALL CPU_TIME(initial_time)
  OPEN(60,FILE='input.dat')
  READ(60,*) potential
  READ(60,*) Vtot
  READ(60,*) Vmax(1:9)

  READ(60,*) coord_in
  READ(60,*) freq_cutoff
  READ(60,*) frequencies_out
  READ(60,*) Nsobol
  READ(60,*) skip
  CLOSE(60)
  freq_cutoff = freq_cutoff /autocm
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
  ALLOCATE(omega(Dim), atom_type(Natoms), sqrt_mass(Dim), mass(Natoms), &
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
!===========================================================================================!
  write(*,*) 'This is omega', omega*autocm
  omegak = 0.
  DO k=1,Dim
     omegak = ABS(omega(k))
     IF(omegak<freq_cutoff) THEN
        omegak = freq_cutoff
!        write(*,*) omegak*autocm   ! testing 
!        write(*,*) 'Use cutoff'    
     ELSE
        omegak = omegak 
        WRITE(*,*) omegak*autocm
        write(*,*) 'Use the original frequency'
     END IF
   END DO
!===========================================================================================!  
!==========================Start old code for PE matrix elements============================!
!======================make quasi_mc fn once it works, for now do it below==================!
!===========================================================================================!
  CALL permutation(dim,Vtot,Vmax)
  ALLOCATE(Umat(Jmax,Jmax), U1(Jmax,Jmax), C(Jmax,Jmax), FV1(Jmax), FV2(Jmax))
  ALLOCATE(eigenvalues(Jmax)) 
  ALLOCATE(scrambled_z(dim), herm(0:Vtot,dim), A(0:Vtot,0:Vtot,dim))
!===========================================================================================!
  Umat=0d0  ! initialize PE matrix
  U1=0d0 
  OPEN(UNIT=80, FILE='matrix.dat')
  OPEN(UNIT=81, FILE='eigenvalues.dat')
!===========================================================================================!
!================================Loop over Sobol Points=====================================!
!=============================Normalize with the Gaussian===================================!
!================Calculate (our normalized) Hermite's, up to Vmax, ini poly=0 ==============!
!===========================================================================================!
  DO i = 1, Nsobol
     CALL sobol_stdnormal(dim,skip,scrambled_z)
     scrambled_z = scrambled_z/SQRT(2.)    ! factor from definition of normal distribution
     herm(0,:) = 1.0                       ! Re-normalized hermite polynomial now
     herm(1,:) = SQRT(2.)*scrambled_z(:)       
     DO j = 2,Vtot      
        herm(j,:) = (SQRT(2./j)*scrambled_z(:)*herm(j-1,:)) - (SQRT(((j-1d0)/j))*herm(j-2,:))
     END DO
!===================================Evaluate Herm * Herm ===================================!
!==========================Matrix A: herm(deg)*herm(deg), deg(0,Vtot)=======================!
!===========================================================================================!
     DO m=0,Vtot
        DO n=0,Vtot
           A(m,n,:) = herm(m,:)*herm(n,:)
        END DO
     END DO
!================================Evaluate Matrix Elements ==================================!
!==================Matrix Umat: PE matrix elements. U1 for partial average==================!
!===========================================================================================!
     DO j=1,Jmax         ! Loop over jmax
        DO j1=j,Jmax 
           B=A(v(1,j),v(1,j1),1)
           DO k=2,dim
              B=B*A(v(k,j),v(k,j1),k)
           END DO
           Umat(j1,j) = Umat(j1,j) + B
        END DO
     END DO
! Write out partial average and flush
     IF(MOD(i,1)==0) THEN
        Umat = Umat + U1
        U1 = 0d0
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
