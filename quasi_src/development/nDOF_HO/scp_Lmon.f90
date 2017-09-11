Module SCP_module
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
  integer :: Dim, Natoms, Dim1, Lmon
  logical withgrad, quantum
  
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
  
  subroutine inerita_tensor(Natoms,coords,IT)
    ! computes the inertia tensor
    
    implicit none
    integer :: i
    integer, intent(in) :: Natoms
    double precision, intent(in) :: coords(3,Natoms)   ! bohr
    
    double precision, intent(out)  :: IT(3,3)            ! amu*bohr**2
    
    IT=0d0
    do i=1,Natoms
       ! diagonal elements
       
       IT(1,1) = IT(1,1) + mass(i)*(coords(2,i)**2 + coords(3,i)**2)  ! I_xx
       IT(2,2) = IT(2,2) + mass(i)*(coords(3,i)**2 + coords(1,i)**2)  ! I_yy
       IT(3,3) = IT(3,3) + mass(i)*(coords(1,i)**2 + coords(2,i)**2)  ! I_zz
       
       ! off-diagonal elements, symmetric matrix
       IT(1,2) = IT(1,2) - mass(i)*coords(1,i)*coords(2,i)  ! I_xy
       IT(2,3) = IT(2,3) - mass(i)*coords(2,i)*coords(3,i)  ! I_yz
       IT(3,1) = IT(3,1) - mass(i)*coords(3,i)*coords(1,i)  ! I_zx
       
    enddo
    
    IT(2,1) = IT(1,2)
    IT(3,2) = IT(2,3)
    IT(1,3) = IT(3,1)
    
  end subroutine inerita_tensor
  
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
  
  
  subroutine Newton_Raphson_Lmon(T,q,Hess,V0,Free_Energy,omega,U,Nsobol,skip,step,freq_cutoff,n_Lmon)
    ! On input:
    ! n_Lmon - the monomer to be SCP'ed
    ! Dim - the number of degrees of freedom = 3*Natoms
    ! q(Dim) - Gaussian Center
    
    ! Output:
    
    ! Hessian(Dim,Dim)- mass-scaled Hessian averaged over the gaussian in atomic units
    ! omega(Dim) - eigenfrequencies from the mass-scaled Hessian
    
    
    implicit none
    integer :: i,k,j,Nsobol,n_Lmon
    integer(kind=8) :: skip
    double precision ::  q(Dim),q1(Dim), G(Dim1,Dim1),G1(Dim1,Dim1),U(Dim1,Dim1),rn(Dim1),&
         r(Dim1),force_r_aver(Dim1,Dim1),&
         Hess(Dim1,Dim1),omegak,d(Dim1),omega(Dim1,0:1),D1(Dim1,Dim1),&
         force_aver(Dim1),force(Dim),Ener,freq_cutoff,T,step,V0,Free_Energy
    
    
!ittm    external calc_pot_link2f90
!ittm    external calc_pot_link2f90_g
    
    G=0d0
    G1=0d0
    do k=1,Dim1
       omegak=dabs(omega(k,0))
       if(omegak<freq_cutoff) omegak=freq_cutoff
       if(omegak<freq_cutoff) then
          d(k)=0d0
          !            write(*,*) 'omega(k)=',dabs(omega(k,0))*autocm,&
          !                  'd(k) -->0 '
       else
          if(quantum) then
             if(T>1d-8) then
                d(k)=0.5/(dtanh(omegak/(2*T))*omegak)
             else
                d(k)=0.5/omegak
             endif
          else
             d(k)=T/omegak**2
          endif
          do i=1,Dim1
             do j=1,Dim1
                G(i,j)=G(i,j)+dsqrt(d(k))*U(i,k)*U(j,k) ! sqrt(D) of the mass-scaled D in atomic units
                !                if(omegak>freq_cutoff) G1(i,j)=G1(i,j)+d(k)*U(i,k)*U(j,k) ! the mass-scaled D in atomic units
                G1(i,j)=G1(i,j)+d(k)*U(i,k)*U(j,k) ! the mass-scaled D in atomic units
             enddo
          enddo
       end if
    enddo
    D1=0d0
    do k=1,Dim1
       do i=1,Dim1
          do j=1,Dim1
             if(d(k)>1d-8)&
                  D1(i,j)=D1(i,j)+U(i,k)*U(j,k)/d(k) ! 1/D of the mass-scaled D-matrix in atomic units
          enddo
       enddo
    enddo
    force_aver=0d0
    force_r_aver=0d0
    V0=0d0
    !     compute the averages
    q1=q
    do k=1,Nsobol
       call sobol_stdnormal(Dim1,skip,rn)
       call rmatvec(Dim1,G,rn,r,.false.) ! multiply by sqrt(D), r in atomic units
       do i=1,Dim1
          r(i)=r(i)/sqrt_mass(Dim1*(n_Lmon-1)+i)
          q1(Dim1*(n_Lmon-1)+i)=q(Dim1*(n_Lmon-1)+i)+r(i)
       enddo
       call water_potential(Natoms/3, q1, Ener, force)
       V0=V0+Ener
       force_aver=force_aver+force(Dim1*(n_Lmon-1)+1:Dim1*n_Lmon)
       do i=1,Dim1
          do j=1,Dim1
             force_r_aver(j,i)=force_r_aver(j,i) + force(Dim1*(n_Lmon-1)+i)*r(j)
          enddo
       enddo
    enddo
    V0=V0/Nsobol
    force_aver=force_aver/Nsobol
    do i=1,Dim1
       do j=1,Dim1
          force_r_aver(j,i)=-force_r_aver(j,i)&
               *sqrt_mass(Dim1*(n_Lmon-1)+j)/(sqrt_mass(Dim1*(n_Lmon-1)+i)*Nsobol)
       enddo
    enddo
    call rmatmat(Dim1,D1,force_r_aver,Hess) ! multiply by D^{-1} to get mass-scaled Hessian in atomic units
    !     Symmetrization:
    do i=1,Dim1
       do j=1,i-1
          Hess(i,j)=(Hess(i,j)+Hess(j,i))/2
          Hess(j,i)=Hess(i,j)
       enddo
    enddo
!    if(Lmon==3.or.Lmon==Dim) call project(Dim1,Hess,q(Dim1*(n_Lmon-1)+1:Dim1*n_Lmon),6)
    if(Lmon==3.or.Lmon==Dim) call project(Dim1,Hess,q(Dim1*(n_Lmon-1)+1:Dim1*n_Lmon),3)
    call frequencies_from_Hess(Dim1,Hess,omega(1,1),U)

    Free_Energy=V0
    do k=1,Dim1
       omegak=dabs(omega(k,1))
       if(omegak > freq_cutoff) then
          if(T>1d-8) then
             Free_Energy=Free_Energy+T*dlog(2*dsinh(omegak/(2*T)))&
                  -omegak/(4*dtanh(omegak/(2*T)))
          else
             Free_Energy=Free_Energy+omegak/4
          endif
       end if
    end do
    
    

!!$    G=0d0
!!$    do k=1,Dim1
!!$       omegak=dabs(omega(k,1))
!!$       if(omegak > freq_cutoff) then
!!$          d(k)=0d0
!!$          if(quantum) then
!!$             if(T>1d-8) then
!!$                d(k)=0.5/(dtanh(omegak/(2*T))*omegak)
!!$             else
!!$                d(k)=0.5/omegak
!!$             endif
!!$          else
!!$             d(k)=T/omegak**2
!!$          endif
!!$          do i=1,Dim1
!!$             do j=1,Dim1
!!$                G(i,j)=G(i,j)+U(i,k)*U(j,k)*d(k) ! the mass-scaled D in atomic units
!!$                !                  G(i,j)=G(i,j)+U(i,k)*U(j,k)/omegak**2 ! 1/K  of the mass-scaled Hessian in atomic units
!!$                !               if(omegak>1d-8) G(i,j)=G(i,j)+U(i,k)*U(j,k)*d(k) ! mass-scaled D in atomic units
!!$             enddo
!!$          enddo
!!$       endif
!!$    enddo
    
    !      write(*,*) 'Force:',force_aver
    if(step>0d0) then
       q(Dim1*(n_Lmon-1)+1:Dim1*n_Lmon)=q(Dim1*(n_Lmon-1)+1:Dim1*n_Lmon)+step*force_aver
    else
       do i=1,Dim1
          do j=1,Dim1
             q(Dim1*(n_Lmon-1)+i)=q(Dim1*(n_Lmon-1)+i)-&
                  step*G1(j,i)/(sqrt_mass(Dim1*(n_Lmon-1)+i)*sqrt_mass(Dim1*(n_Lmon-1)+j))*force_aver(j) ! "Newton-Raphson" using unscaled K or D
          enddo
       enddo
    endif
    return
  end subroutine Newton_Raphson_Lmon


  subroutine Moments(T,q,omega,U,Nsobol,skip,freq_cutoff)
    ! On input:
    ! n_Lmon - the monomer to be SCP'ed
    ! Dim - the number of degrees of freedom = 3*Natoms
    ! q(Dim) - Gaussian Center
    
    ! Output:
    
    ! Hessian(Dim,Dim)- mass-scaled Hessian averaged over the gaussian in atomic units
    ! omega(Dim) - eigenfrequencies from the mass-scaled Hessian
    
    
    implicit none
    integer :: i,k,j,Nsobol
    integer(kind=8) :: skip
    double precision ::  q(Dim), G(Dim1,Dim1),U(Dim1,Dim1),rn(Dim1),&
         r(Dim1),principals(3),IT_aver(3,3),IT(3,3),omegak,d(Dim1),omega(Dim1),freq_cutoff,T
    
    
    call CM(Natoms,q)
    call inerita_tensor(Natoms,q,IT)
    IT=IT*bohr**2/melectron
    call RS(3, 3, IT, principals, 0,r(1) , r(2:4), r(5:7), i)
    write(*,*) 'Principal moments for single config (amu*A**2):',principals
    write(*,*) 'Rotational constants for single config (MHz):', (505379.045961437d0/principals(i), i=1,3)
    !     Diagonalize the G-matrix G=U rho U^T
    G=0d0
    d=0
    do k=1,Dim1
       omegak=dabs(omega(k))
       if(omegak>freq_cutoff) then
          if(quantum) then
             if(T>1d-8) then
                d(k)=0.5/(dtanh(omegak/(2*T))*omegak)
             else
                d(k)=0.5/omegak
             endif
          else
             d(k)=T/omegak**2
          endif
          do i=1,Dim1
             do j=1,Dim1
                G(i,j)=G(i,j)+dsqrt(d(k))*U(i,k)*U(j,k) ! sqrt(D) of the mass-scaled D in atomic units
             enddo
          enddo
       end if
    enddo
    IT_aver=0d0
    !     compute the averages
    do k=1,Nsobol
       call sobol_stdnormal(Dim1,skip,rn)
       call rmatvec(Dim1,G,rn,r,.false.) ! multiply by sqrt(D), r in atomic units
       do i=1,Dim1
          r(i)=q(i)+r(i)/sqrt_mass(i)
       enddo
       call inerita_tensor(Natoms,r,IT)
       IT_aver=IT_aver+IT
    enddo
    IT_aver=(IT_aver/Nsobol)*bohr**2/melectron
      ! diagonalize
    call RS(3, 3, IT_aver, principals, 0,r(7) , r(1:3), r(4:6), i)
    write(*,*) 'Vibrationally averaged principal moments (amu*A**2):',principals
    write(*,*) 'Vibrationally averaged Rotational constants (MHz):', (505379.045961437d0/principals(i), i=1,3)
    return
  end subroutine Moments



  subroutine frequencies_from_Hess(Dim1,Hess,omega,U)
    implicit none
    integer :: i,IERR,Dim1
    double precision :: Hess(Dim1,Dim1),A(Dim1,Dim1),omega(Dim1),FV1(Dim1),FV2(Dim1),U(Dim1,Dim1)
    A=Hess
    call RS(Dim1,Dim1,A,omega,1,U,FV1,FV2,IERR) 
          write(*,*) 'Frequencies from the Hessian:'
          write(*,*) 
    do i=Dim1,1,-1
       omega(i)=sign(dsqrt(dabs(omega(i))),omega(i))
!       if(dabs(omega(i))>1d-8) write(*,*) omega(i)*autocm
       write(*,*) omega(i)*autocm
    enddo
    
  end subroutine frequencies_from_Hess
  
  
  
  subroutine Get_Hessian(q,H)

    implicit none
    integer :: i,j
    double precision ::  H(Dim,Dim),q(Dim),r(Dim),E,force(Dim),force0(Dim)
    double precision, parameter :: s=1d-5
    
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
   
  subroutine rtranspose(N,A)
    implicit none
    integer :: N,i,j
    double precision ::  A(N,N),b
    do i=1,N
       do j=1,i-1
          b=A(i,j)
          A(i,j)=A(j,i)
          A(j,i)=b
       enddo
    enddo
    
  end subroutine rtranspose
  
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
  
end Module SCP_module

program Self_Consistent_Phonons
  
  USE SCP_module
  
  IMPLICIT NONE
  
  real :: initial_time, final_time  
  integer(kind=8) :: skip
  integer :: Nsobol,Nsobolmin,Nsobolmax,i,j,Niter,iter,n_Lmon,j_max
  double precision, allocatable :: omega(:,:), Hess(:,:), U(:,:), H1(:,:), q(:),q0(:),force(:),omega_Lmon(:)
  double precision :: freq_cutoff,T,TK,E,step,w_max,V0,Free_Energy,principals(3)
  character (len=1) status
  character(len=50) coord_in, coord_out, frequencies_out
  
!ittm  external calc_pot_link2f90
!ittm  external calc_pot_link2f90_g
  
  CALL CPU_TIME(initial_time)
  ! Read initial configuration
  skip=1048576                      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2**20 
  read(*,*) Natoms
  Dim=3*Natoms
  read(*,*) quantum,withgrad,potential
  read(*,*) Nsobolmin,Nsobolmax
  read(*,*) Niter, step
  read(*,*) freq_cutoff
  if(step<0d0) write(*,*) 'Using the Newton-Raphson method'
  freq_cutoff=freq_cutoff*cmtoau
  read(*,*) TK
  T=TK*Ktoau   ! T in atomic units
  read(*,*) Lmon
  Dim1=9
  if(Lmon==3) then
     write(*,*) '#  Lmon-3'
  else if(Lmon==9) then
     write(*,*) '#  Lmon-9'
  else
     Lmon=Dim
     Dim1=Dim
     write(*,*) '# Full SCP'
  endif
  read(*,*) status
  read(*,*) coord_in
  read(*,*) coord_out
  read(*,*) frequencies_out
  allocate(omega(Dim,0:Niter), atom_type(Natoms), sqrt_mass(Dim),mass(Natoms), &
       q(Dim), q0(Dim), force(Dim), Hess(Dim,Dim),  H1(Dim1,Dim1), &
       U(Dim,Dim),omega_Lmon(Dim))
  open(1,file=coord_in)
  read(1,*) Natoms
  read(1,*)
  write(*,*) 'Natoms=',Natoms
  do i=1,Natoms
     read(1,*) atom_type(i), q0(3*i-2:3*i)   ! coordinates in Angstroms
     mass(i)=Atom_mass(atom_type(i))
     sqrt_mass(3*i-2:3*i)=sqrt(mass(i))
  enddo
  write(*,*) 'molecule:', atom_type
  q0=q0/bohr            ! convert coordinates to atomic units
  call water_potential(Natoms/3, q0, E, force)
  write(*,*) 'E0 =',E*autokcalmol, 'kcal/mol'
  write(*,*) 'force ='
  do i=1,Natoms
     write(*,*) force(3*i-2:3*i)
  enddo
!stop
    if (status=='n') then
       call Get_Hessian(q0,Hess)          
!    write(*,*) 'Hess ='
!    do i=1,Dim
!       write(*,*) Hess(1:Dim,i)
!    enddo
    else
       read(1,*) ((Hess(i,j),i=1,j),j=1,Dim) ! mass-scaled Hessian in atomic units
       do j=1,Dim
          do i=1,j-1
             Hess(j,i)=Hess(i,j)
          enddo
       enddo
    endif
    close(1)
    open(7,file=frequencies_out)
!    call CM_zero(nc,nh,q0,q_cm)
!    call project(Dim,Hess,q0,6)
!    write(*,*) '# Fully coupled normal modes:'
!    call frequencies_from_Hess(Dim,Hess,q,U1)  ! here q is used for frequencies
!    open(33,file='HA_freq.dat')
!    do i=1,Dim
!       write(33,*) i,q(i)*autocm
!    enddo
!    close(33)
    write(*,*) 
    write(*,*) 'T=',TK
    write(*,*) 
    write(*,*) '# SCP-Local monomer modes:'

    do n_Lmon=1,Dim/Dim1
       q=q0
       write(*,*) '# n_Lmon=',n_Lmon
       write(7,*) '# n_Lmon=',n_Lmon
       H1=Hess((n_Lmon-1)*Dim1+1:n_Lmon*Dim1,(n_Lmon-1)*Dim1+1:n_Lmon*Dim1)
!       if(Lmon==3.or.Lmon==Dim) call project(Dim1,H1,q(Dim1*(n_Lmon-1)+1:Dim1*n_Lmon),6)
       if(Lmon==Dim) call project(Dim1,H1,q(Dim1*(n_Lmon-1)+1:Dim1*n_Lmon),3)
       call frequencies_from_Hess(Dim1,H1,omega(Dim1*(n_Lmon-1)+1:Dim1*n_Lmon,0),U)
       skip=1048576 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2**20 ! 1000
       ! For each monomer n_Lmon do 9x9 SCP with all the other monomers fixed
       do iter=1,Niter
          Nsobol=(Nsobolmin*(Niter-iter)+Nsobolmax*(iter-1))/(Niter-1)
          write(*,*) 'iter=',iter,'  Nsobol=',Nsobol
          call Newton_Raphson_Lmon(T,q,H1,V0,Free_Energy,omega(Dim1*(n_Lmon-1)+1:Dim1*n_Lmon,iter-1:iter),&
               U,Nsobol,skip,step,freq_cutoff,n_Lmon)
          if(Lmon==Dim) then
             write(*,*) 'T=',TK,'  Free Energy=',Free_Energy*autokcalmol, 'kcal/mol'
             call Moments(T,q,omega,U,Nsobolmax,skip,freq_cutoff)
          endif

       enddo

       do i=1,Dim1
          do j=0,Niter
             !          if (dabs(omega(i,j) > 1d-8) write(7,*) j,omega(i,j)*autocm
             write(7,*) j,omega(Dim1*(n_Lmon-1)+i,j)*autocm
          enddo
          write(7,*)
       enddo
! Do Lmon=3 or Lmon=Dim; 
       call project(Dim1,H1,q,3)
       write(*,*) 'Using Eckart subspace...'
       call frequencies_from_Hess(Dim1,H1,omega_Lmon(Dim1*(n_Lmon-1)+1:Dim1*n_Lmon),U)
    enddo
    write(*,*)
    write(*,*) 'Lmon=',Lmon, 'Frequencies:'
    do i=1,Natoms
       w_max=-1d8
       do j=1,Dim
          if(omega(j,Niter)>w_max) then
             w_max=omega(j,Niter)
             j_max=j
          endif
       enddo
       omega(j_max,Niter)=-1d9
       write(*,*) w_max*autocm
    enddo

       Write(*,*) 'Frequencies in Eckart subspace:'
       do i=1,Natoms
          w_max=-1d8
          do j=1,Dim
             if(omega_Lmon(j)>w_max) then
                w_max=omega_Lmon(j)
                j_max=j
             endif
          enddo
          omega_Lmon(j_max)=-1d9
          write(*,*) w_max*autocm
       enddo


    if(Lmon==Dim) then
       open(1,file=coord_out)
       write(1,*) Natoms
       write(1,*)  'T=',TK,'  Free Energy=',Free_Energy*autokcalmol, 'kcal/mol'
       do i=1,Natoms
          write(1,*) atom_type(i), q(3*i-2:3*i)*bohr  ! coordinates in Angstroms
       enddo
       write(1,*) ((H1(i,j),i=1,j),j=1,Dim) 
       close(1)
    endif
    CALL CPU_TIME(final_time)
    WRITE(*,*) 'TOTAL TIME: ', final_time - initial_time
  end program Self_Consistent_Phonons
  
  
  
