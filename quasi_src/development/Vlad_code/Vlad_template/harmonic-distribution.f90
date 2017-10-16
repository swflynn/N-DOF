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
  integer :: Dim, Natoms, Dim1, Lmon   ! for a single water monomer first try Natom=3, Dim=9, Dim1=6
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
  
  
  subroutine quasi_MC(q,omega,U,Nsobol,skip,freq_cutoff)

    implicit none
    integer :: i,k,j,Nsobol
    integer(kind=8) :: skip
    double precision ::  q(Dim),q1(Dim), G(Dim1,Dim1),G1(Dim1,Dim1),U(Dim1,Dim1),y(Dim1),&
         r(Dim1),omegak,d(Dim1),omega(Dim1),Ener,freq_cutoff
       
!ittm    external calc_pot_link2f90
!ittm    external calc_pot_link2f90_g

    ! for a single water monomer first try Natom=3, Dim=9, Dim1=6   

    do k=Dim-Dim1+1,Dim
       omegak=dabs(omega(k,0))
       if(omegak<freq_cutoff) omegak=freq_cutoff
       do i=1,Dim
          U(i,k)=sqrt(0.5/omegak)/sqrt_mass(i)*U(i,k)   ! the sqrt(0.5) factor comes from using the standard normal distr exp(-y^2/2)  
          ! to sample points with the distribution exp(-y^2)
    enddo

    V0=0d0
    !     compute the averages
    do k=1,Nsobol
       call sobol_stdnormal(Dim1,skip,y)
       q1=q
       do k=Dim-Dim1+1,Dim   ! sum over Dim1 largest eigenvalues
          q1(:)=q1(:)+y(k)*U(:,k)
       enddo
       call water_potential(Natoms/3, q1, Ener, force)
       V0=V0+Ener
    enddo
    V0=V0/Nsobol
       
  end subroutine quasi_MC



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

PROGRAM nDOF_HO
  USE SCP_module
  
  IMPLICIT NONE
  
  real :: initial_time, final_time  
  integer(kind=8) :: skip
  integer :: Nsobol,i,j
  double precision, allocatable :: omega(:,:), Hess(:,:), U(:,:), H1(:,:), q(:),q0(:)
  double precision :: freq_cutoff
  character(len=50) coord_in,frequencies_out
  
!ittm  external calc_pot_link2f90
!ittm  external calc_pot_link2f90_g
  
  CALL CPU_TIME(initial_time)
  ! Read initial configuration
  skip=1000 
  read(*,*) Natoms, Dim1
  Dim=3*Natoms
  read(*,*) quantum,withgrad,potential
  read(*,*) Nsobol
  read(*,*) coord_in
  read(*,*) frequencies_out
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
  allocate(omega(Dim), atom_type(Natoms), sqrt_mass(Dim), mass(Natoms), &
       q(Dim), q0(Dim), force(Dim), Hess(Dim,Dim),  H1(Dim1,Dim1), &
       U(Dim,Dim))
  call Get_Hessian(q0,Hess)          
  open(7,file=frequencies_out)
  call frequencies_from_Hess(Dim,Hess,omega,U)
  skip=Nsobol !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2**20 ! 1000
  call  quasi_MC(q,omega,U,Nsobol,skip,freq_cutoff)
  CALL CPU_TIME(final_time)
  WRITE(*,*) 'TOTAL TIME: ', final_time - initial_time
END PROGRAM nDOF_HO
