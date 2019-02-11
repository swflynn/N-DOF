!=============================================================================80
!                 Distributed Gaussian Basis (Ground State Energy)
!==============================================================================!
!    Discussion:
!DGB analysis for 1D potential: single atom (x coordinate only, mass:=1)
!Equations formulated assuming normal mode coordinates
!Gaussian chosen with quasi-random grid 
!Gaussian Width (alpha) computed as a function of r
!==============================================================================!
!       Modified:
!   10 December 2018
!       Author:
!   Shane Flynn 
!==============================================================================!
module dgb_groundstate
implicit none
!==============================================================================!
!                            Global Paramaters
!==============================================================================!
double precision,parameter::Hmass=1d0
double precision,parameter::pi=4.*atan(1d0)
!==============================================================================!
!                            Global Variables 
!==============================================================================!
!       Discussion:
!Natoms             ==> Number of Atoms 
!Dimen              ==> System Dimensionality 
!==============================================================================!
integer::Natoms,Dimen
character(len=2),allocatable::atom_type(:)
double precision,allocatable::mass(:),sqrt_mass(:)
!==============================================================================!
contains
!==============================================================================!
function Atom_Mass(atom)
!       Discussion:
!==============================================================================!
implicit none
double precision::Atom_Mass
character(len=2)::atom
if(atom=='H'.or.atom=='h')then 
    Atom_mass=Hmass
else 
    write(*,*) 'atom ', atom, ' is not recognized'
    stop 'Check Atom_Mass Function'
endif
end function Atom_Mass
!==============================================================================!
subroutine Toy_Potential(x,energies)
!==============================================================================!
!       Discussion:
!Hard-coded Potential Energy 
!V:=1/2*(x)^2
!==============================================================================!
implicit none
double precision::x(Dimen),energies
energies=0.5*x(1)**2
end subroutine Toy_Potential
!==============================================================================!
subroutine Toy_Force(x,forces)
!==============================================================================!
!       Discussion:
!Returns the Forces associated with Toy_Potential Subroutine
!Forces are hard-coded based on Toy_Potential 
!==============================================================================!
implicit none
integer::i
double precision::x(Dimen),forces(Dimen)
forces(1)=-x(1)
end subroutine Toy_Force
!==============================================================================!
subroutine Toy_Hessian(x,Hess_Mat)
!==============================================================================!
!       Discussion:
!Numerically computed Hessian using forces from Toy_Force Subroutine
!s          ==> Perturbation Parameter 
!Hess_Mat   ==> (Dimen,Dimen); Symmetrized Mass-Scaled Hessian
!x          ==> (Dimen); XYZ Configuration at minimum
!==============================================================================!
implicit none 
integer::i,j
double precision::Hess_Mat(Dimen,Dimen),x(Dimen),r(Dimen),force(Dimen)
double precision::force0(Dimen)
double precision,parameter::s=1d-6
r=x
call Toy_Force(r,force0)
do i=1,Dimen
    r(i)=x(i)+s
    call Toy_Force(r,force)
    r(i)=x(i)
    do j=1,Dimen
        Hess_Mat(i,j)=(force0(j)-force(j))/s
    enddo
enddo
!==============================================================================!
!                   Symmetrize and Mass-Scale the Hessian
!==============================================================================!
do i=1,Dimen
    do j=1,i
        if(i.ne.j) Hess_Mat(i,j)=(Hess_Mat(i,j)+Hess_Mat(j,i))/2
        Hess_Mat(i,j)=Hess_Mat(i,j)/(sqrt_mass(i)*sqrt_mass(j))
        if(i.ne.j) Hess_Mat(j,i)=Hess_Mat(i,j)
    enddo
enddo
end subroutine Toy_Hessian
!==============================================================================!
subroutine Freq_Hess(Dimen,Hess,omega,U)
!==============================================================================!
!       Discussion:
!Compute Eigenvalues and Eigenvectors of Hessian
!Uses the LLAPACK real symmetric eigen-solver (dsygev)
!Hess   ==> (Dimen,Dimen); Hessian Matrix
!omega  ==> (Dimen); Hessian Eigenvalues 
!U      ==> (Dimen,Dimen); Hessian Eigenvectors
!       LLAPACK:
!v      ==> Compute both Eigenvalues and Eigenvectors
!u      ==> Use Upper-Triangle of matrix
!==============================================================================!
implicit none
integer::i,info,lwork,Dimen
double precision::Hess(Dimen,Dimen),omega(Dimen),U(Dimen,Dimen)
double precision,allocatable::work(:) 
lwork = max(1,3*Dimen-1)        !suggested by LAPACK Developers
allocate(work(max(1,lwork)))    !suggested by LAPACK Developers 
U=Hess  
call dsyev('v','u',Dimen,U,Dimen,omega,work,lwork,info) 
write(*,*) 'Frequencies from the Hessian:'
do i=Dimen,1,-1
    omega(i)=sign(sqrt(abs(omega(i))),omega(i))
    write(*,*) omega(i), 'normalized = 1?', sum(U(:,i)**2)
enddo
end subroutine Freq_Hess
!==============================================================================!
end module dgb_groundstate
!==============================================================================!
!==============================================================================!
program DGB_1D
use dgb_groundstate
!==============================================================================!
!==============================================================================!
!               Discussion:
!coord_in       ==> Input Geometry File-Name (xyz file)
!NG             ==> Number of Gaussians 
!Nsobol         ==> Number of Sobol Points for numerical integration
!ii,skip        ==> Integer(kind=8): necessary for using Sobol Module
!alpha          ==> (NG): Width of Gaussian
!q0             ==> (Dimen): Input Geometry coordinates 
!r              ==> (Dimen,NG): Gaussian Centers 
!r2             ==> (Dimen): ith Gaussian coordinate for integration
!force          ==> (Dimen): Forces. See "Toy_Force" 
!Hess           ==> (Dimen,Dimen): Mass-Scaled Hessian. See "Toy_Hessian"
!omega          ==> (Dimen): Frequencies (Eigenvalues). See "Freq_Hess"
!U              ==> (Dimen,Dimen): Hess Eigenvectors. See "Freq_Hess"
!z              ==> (Dimen,NG): Quasi-Random Sequence for Gaussian Centers
!z2             ==> (Dimen,Nsobol): Quasi-Random Sequence for Integration 
!Smat           ==> (NG,NG): Overlap Matrix
!Tmat           ==> (NG,NG): Kinetic Energy Matrix
!Vmat           ==> (NG,NG): Potential Energy Matrix
!Hmat           ==> (NG,NG): Hamiltonian Matrix (V+T) for eigenvalue problem
!eigenvalues    ==> (NG): Eigenvalues of the Hamiltonian Matrix
!lambda         ==> (NG): Eigenvalues of the Overlap Matrix
!r_ij           ==> center of the i,j matrix element 
!E0             ==> Energy Evaluation at the minimum configuration
!pot_ene        ==> Potential Energy evaluated in q space
!==============================================================================!
implicit none
character(len=50)::coord_in
integer::NG,Nsobol
integer::i,j,k,l,n
integer*8::ii,skip,skip2
double precision::E0,pot_ene,s_sum,alpha0
!==============================================================================!
double precision,allocatable::q0(:),force(:),r(:,:),r2(:),eigenvalues(:)
double precision,allocatable::Hess(:,:),omega(:),U(:,:),z(:,:),z2(:,:)
double precision,allocatable::Smat(:,:),alpha(:),Tmat(:,:),Vmat(:,:)
double precision,allocatable::Hmat(:,:),S1mat(:,:),lambda(:),r_ij(:)
!==============================================================================!
!                       LLAPACK dsygv variables                                !
!==============================================================================!
integer::itype,info,lwork
double precision,allocatable::work(:)
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
read(*,*) coord_in
read(*,*) NG
read(*,*) Nsobol
read(*,*) alpha0
skip=NG
skip2=Nsobol
write(*,*) 'Test 1; Successfully Read Input Data File'
!==============================================================================!
!                         Set Input Water Geometry 
!==============================================================================!
open(16,File=coord_in)
read(16,*) Natoms
read(16,*) 
Dimen=1*Natoms 
write(*,*) 'Dimensionality ==> ', Dimen
!==============================================================================!
allocate(atom_type(Natoms),mass(Natoms),sqrt_mass(Dimen),q0(Dimen),force(Dimen))
allocate(Hess(Dimen,Dimen),omega(Dimen),U(Dimen,Dimen),z(Dimen,NG))
allocate(z2(Dimen,Nsobol),r(Dimen,NG),alpha(NG),Smat(NG,NG),eigenvalues(NG))
allocate(Tmat(NG,NG),Vmat(NG,NG),Hmat(NG,NG),r2(Dimen))
allocate(S1mat(NG,NG),lambda(NG),r_ij(Dimen))
!==============================================================================!
!                         Input Configuration Energy
!==============================================================================!
do i=1,Natoms
    read(16,*) atom_type(i),q0(1)   
    mass(i)=Atom_mass(atom_type(i))
    sqrt_mass(1)=sqrt(mass(i))
enddo
close(16)
call toy_potential(q0,E0)
write(*,*) 'E0 ==> ', E0
!==============================================================================!
! 			Compute Hessian and Frequencies
!==============================================================================!
call Toy_Hessian(q0, Hess)
call Freq_Hess(Dimen,Hess,omega,U)
write(*,*) 'Test 2; Successfully Computed Hessian'
!==============================================================================!
!               Generate Gaussian Centers with quasi-random grid
!==============================================================================!
do ii=1,NG
    call sobol_stdnormal(Dimen,skip,z(:,ii))
    do i=1,Dimen
        r(:,ii)=z(:,ii)/sqrt(abs(omega(i)))
    enddo
enddo
!==============================================================================!
!                       Generate Alpha Scaling 
!==============================================================================!
do i=1,NG
    alpha(i)=alpha0*NG**(2./dimen)*exp(-1./dimen*sum(r(:,i)**2*omega(:)))
enddo
open(unit=17,file='centers.dat')
do i=1,NG
    write(17,*) alpha(i), r(1,i)
enddo
close(17)
!==============================================================================!
!                           Overlap Matrix (S)
!==============================================================================!
do i=1,NG
    do j=i,NG
        s_sum=sum(omega(:)*(r(:,i)-r(:,j))**2)               
        Smat(i,j)=sqrt(alpha(i)+alpha(j))**(-dimen)&           
                 *exp(-0.5*alpha(i)*alpha(j)/(alpha(i)+alpha(j))*s_sum)
        Smat(j,i)=Smat(i,j)
    enddo
enddo
!==============================================================================!
!                   Check to see if S is positive definite
!If this is removed, you need to allocate llapack arrays before Hamiltonian 
!==============================================================================!
S1mat=Smat
lwork=max(1,3*NG-1)
allocate(work(max(1,lwork)))
call dsyev('v','u',NG,S1mat,NG,lambda,work,Lwork,info)
write(*,*) 'Info (Overlap Matrix) ===>', info
open(unit=18,file='overlap_eigenvalues.dat')
do i=1,NG
    write(18,*) lambda(i)
enddo
close(18)
write(*,*) 'Test 3; Successfully computed Overlap Matrix'
!==============================================================================!
!                           Kinetic Matrix (T)
!==============================================================================!
do i=1,NG
    do j=i,NG
        Tmat(i,j)=Smat(i,j)*0.5*alpha(i)*alpha(j)/(alpha(i)+alpha(j))&
        *sum(omega(:)-(alpha(i)*alpha(j)*(omega(:)**2&
        *(r(:,i)-r(:,j))**2)/(alpha(i)+alpha(j))))
        Tmat(j,i)=Tmat(i,j)
    enddo
enddo
write(*,*) 'Test 4; Successfully computed Kinetic Matrix'
!==============================================================================!
!                   Generate Sequence For Evaluating Potential
! Want to generate a single sequence and then scale for each Gaussian
!==============================================================================!
do ii=1,Nsobol
    call sobol_stdnormal(Dimen,skip2,z2(:,ii))
enddo
write(*,*) 'Test 5; Successfully generated integration sequence'
!==============================================================================!
!                           Evaluate Potential 
!==============================================================================!
Vmat=0d0
do i=1,NG
   do j=i,NG
        r_ij(:)=(alpha(i)*r(:,i)+alpha(j)*r(:,j))/(alpha(i)+alpha(j))
        do l=1,Nsobol
            r2(:)=r_ij(:)+z2(:,l)/sqrt(omega(:)*(alpha(i)+alpha(j)))
            call Toy_Potential((q0+matmul(U,r2)/sqrt_mass(:)),pot_ene)
            Vmat(i,j)=Vmat(i,j)+pot_ene
        enddo
        Vmat(j,i)=Vmat(i,j)
    enddo
enddo
Vmat=Vmat*Smat/Nsobol
write(*,*) 'Test 6; Successfully Computed Potential Matrix'
!==============================================================================!
!                     Solve Generalized Eigenvalue Problem
!==============================================================================!
Hmat=Vmat
Hmat=Hmat+Tmat
itype=1
eigenvalues=0d0
!==============================================================================!
! Allocations needed if overlap is not diagonalized above
!==============================================================================!
!lwork=max(1,3*NG-1)                                                           !
!allocate(work(max(1,lwork)))                                                  !
!==============================================================================!
CALL dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info) 
write(*,*) 'Info (Hamiltonian) ==>', info
open(unit=19,file='eigenvalues.dat')
do i=1,NG
    write(19,*) eigenvalues(i), ((i-1)+.5)*omega(1), &
        ((eigenvalues(i) - ((i-1)+.5)*omega(1))/(((i-1)+.5)*omega(1)))*100
enddo
close(19)
!==============================================================================!
!                               output file                                    !
!==============================================================================!
open(90,file='simulation.dat')
write(90,*) 'Natoms ==>', Natoms
write(90,*) 'Dimensionality ==>', Dimen 
write(90,*) 'N_gauss ==>', NG
write(90,*) 'N_Sobol ==>', Nsobol
write(90,*) 'omega ==>', omega
close(90)
end program DGB_1D
