!=============================================================================80
!                       Lmon Single Water Implementation                       !
!==============================================================================!
!    Discussion:
!               QMC code will need some of the pieces as a template
!               can remove this code anytime, keep seperate from qmc project    
!    Modified:
!       14 April 2018
!    Author:
!       Shane Flynn & Vladimir Mandelshtam
!==============================================================================!
module Lmon_water_qMC
implicit none
!==============================================================================!
!                            Global Parameters                                 !
!==============================================================================!
!    integers: 
!       Dim1   => the dimensionality of the monomer, 9 for water 
!==============================================================================!
double precision, parameter :: bohr = 0.52917721092
double precision, parameter :: autocm = 2.194746313D5
double precision, parameter :: autokcalmol = 627.5096 
double precision, parameter :: melectron=1822.88839
double precision, parameter :: Hmass = 1.00782503223*melectron
double precision, parameter :: Omass = 15.99491461957*melectron
integer, parameter :: Dim1 = 9
!==============================================================================!
!                            Global Variables                                  !
!==============================================================================!
!    integers: 
!       Natoms => number of atoms read in from xyz file 
!       Dim    => the total dimensionality (xyz) of the system; 3*Natoms 
!       Jmax   => (Integer) the number of permutations for the system
!       v      => (Dim1,Jmax) all of the combinatorics
!       Vtot   => 
!       vmax   => 
!==============================================================================!
integer :: Natoms, Dim, Dim_LM, Jmax, Vtot
integer, allocatable :: v(:,:), Vmax(:)
character(len=5) :: potential
double precision, allocatable :: sqrt_mass(:),mass(:)
character(len=2), allocatable :: atom_type(:)
!===============================OMP Parallel===================================!
integer :: Nthreads, TID
!==============================================================================!
!                               Module                                         !
!==============================================================================!
CONTAINS
!==============================================================================!
!==============================================================================!
!==============================================================================!
function Atom_Mass(atom)
!==============================================================================!
!    Discussion:
!       Atom_mass computes the mass for atoms read in from file.
!       Currently only for H, O atoms, because we work with a water cluster.
!    Arguments:
!       atom      => 1 or 2 characters associated with the atom type (H,O)
!    Output:
!       Atom_mass => double percision numerical value of the atoms mass
!==============================================================================!
implicit none
double precision :: Atom_Mass
character(len=2) :: atom
if (atom == 'O'.or. atom == 'o') then  
    Atom_mass=Omass
else if (atom == 'H' .or. atom == 'h') then 
    Atom_mass=Hmass
else 
    write(*,*) 'atom ', atom, ' is not recognized'
    write(*,*) 'Check Atom_Mass Function'
    STOP 
end if
end function Atom_Mass
!==============================================================================!
!==============================================================================!
!==============================================================================!
subroutine water_potential(NO,q,energy,force)
!==============================================================================!
!    Discussion:
!       water_potential computes the energy and force for either the TIP4P or
!       MBPOL Potential Energy Surfaces. 
!    Arguments:
!       NO      => Integer; The number of water molecules in the cluster.
!       q       => DP Array: Coordinates for all of the atoms in the cluster
!                  Every water has 3 atoms, every atom has 3 coordinates; 9*NO
!       energy  => DP: the energy evaluated by the PES for given coordinates
!       force   => DP array: Forces associated with all of the coordinates
!    Output:
!       force and energy associated with the input coordinates
!    Other:
!       iso_c_binding => intrinsic Fortran module for linking C code and Fortran
!                     => The MBPOL potential is written in C and must be called
!       TIP4P_module  => Fortran Module for the TIP4P PES, see TIP4P.f90 for src
!==============================================================================!
use iso_c_binding
use TIP4P_module
implicit none
integer :: NO                              
double precision :: energy
double precision, DIMENSION(9*NO) :: q, force
if(potential =='tip4p' .or. potential == 'TIP4P') then
    call TIP4P(NO, q, energy, force)
else if(potential == 'mbpol' .or. potential == 'MBPOL') then
    call calcpotg(NO, energy, q*bohr, force)
    force = -force*bohr/autokcalmol
    energy = energy/autokcalmol
else
    write(*,*) 'cannot identify the potential'
    write(*,*) 'Check water_potential subroutine'
    STOP
end if
end subroutine water_potential
!==============================================================================!
!==============================================================================!
!==============================================================================!
subroutine Get_Hessian(q,Hess)
!==============================================================================!
!    Discussion:
!       Get_Hessian computes the Hessian Matrix for a given set of coordinates
!       This is done through computing the forces via the water potential
!       Forces are then computed for a small pertenation of the coordinates
!       The Hessian is found through the difference (numerical 2nd derivative)
!    Arguments:
!       q(Dim)               => coordinates for every atom (xyz)
!       Hess(Dim1,Dim1)      => The Hessian 
!    Output:
!       Hess(Dim1,Dim1)      => Mass Scaled Hessian (square-symmetric)
!    Other:
!       s                    => DP: finite difference parameter
!==============================================================================!
implicit none
integer :: i,j
double precision ::  Hess(Dim1,Dim1), q(Dim), r(Dim), E, force(Dim), force0(Dim)
double precision, parameter :: s=1d-6
r=q
call water_potential(Natoms/3, r, E, force0)
do  i=1,Dim1
    r(i)=q(i)+s
    call water_potential(Natoms/3, r, E, force)
    r(i)=q(i)
    do j=1,Dim1
        Hess(i,j)=(force0(j)-force(j))/s
    end do
end do
! Symmetrize and Mass-Scale the Hessian
do i=1,Dim1
    do j=1,i
        if(i.ne.j) Hess(i,j)=(Hess(i,j)+Hess(j,i))/2
        Hess(i,j)=Hess(i,j)/(sqrt_mass(i)*sqrt_mass(j))
        if(i.ne.j)  Hess(j,i)=Hess(i,j)
    end do
end do
end subroutine Get_Hessian
!==============================================================================!
!==============================================================================!
!==============================================================================!
!==============================================================================!
subroutine Frequencies_From_Hess(Dim1,Hess,omega,U)
!==============================================================================!
!    Discussion:
!       Frequencies_From_Hess uses the Hessian to compute the associated
!       Eigenvalues (frequencies) and eigenvectors (normal modes)
!       LAPACK is used to compute eigenvalues and eigenvectors
!       dsyev is the real symmetric subroutine in LAPACK
!    Arguments:
!       Dim1    => (Integer) The dimensionality of the monomer
!       Hess    => (Dim1,Dim1) The Hessian 
!       omega   => (Dim1) the Eigenvalues of the Hessian
!       U       => (Dim1,Dim1) the Eigenvectors of the Hessian
!    Output:
!       omega   => The frequencies (eigenvalues)
!       U       => The eigenvectors (normal modes)
!    Other:
!       dsyev   => v: eigenvalues and eigenvectors, u: upper triangle of matrix
!                  Dim1: matrix order, U(Dim1,Dim1): matrix, Dim1: leading dim
!                  omega: eigenvalues, work: algorithm, lwork: algorithm 
!                  info: check matrix properties. ==> See LAPACK Documentation
!                  lwork is suggested to be max(1,3*N-1) by LAPACK authors
!==============================================================================!
implicit none
integer :: i, info, lwork, Dim1
double precision :: Hess(Dim1,Dim1), omega(Dim1), U(Dim1,Dim1),x
double precision, allocatable :: work(:) 
lwork = max(1,3*Dim1-1)
allocate(work(max(1,lwork)))
! dsyev replaces original matrix with eigenvectors
U=Hess  
call dsyev('v', 'u', Dim1, U, Dim1, omega, work, lwork, info) 
write(*,*) 'Frequencies from the Hessian:'
! Switch to the descending order
do i=1,Dim1/2
   x=omega(i)
   omega(i)=omega(Dim1-i+1)
   omega(Dim1-i+1)=x
   work(1:Dim1)=U(1:Dim1,i)
   U(1:Dim1,i)=U(1:Dim1,Dim1-i+1)
   U(1:Dim1,Dim1-i+1)=work(1:Dim1)
enddo
do i=1,Dim1
    omega(i)=sign(sqrt(abs(omega(i))),omega(i))
    write(*,*) omega(i)*autocm !, 'normalized = 1?',sum(U(:,i)**2)
end do
end subroutine frequencies_from_Hess
!==============================================================================!
!==============================================================================!
!==============================================================================!
subroutine permutation()
!==============================================================================!
!    Discussion:
!       permutation computes the number of permutations associated with a given
!       number of basis functions, and excitation in each basis, and a total
!       excitation for the system. The subroutine computes the total number of
!       combinations, and then computes the sequence of excitations in each 
!       basis for each combination. A water monomer is completly described by 9 
!       basis functions, therefore Dim1>9 is not considered
!    Output:
!       Jmax    => (Integer) the total number of permutations
!       v       => (dim1,Jmax) all of the different combinatorics
!                  for 3 basis functions v(3,1) = (0 0 0)
!    Other:
!       This code can most probably be written in a more efficient manner.
!==============================================================================!
implicit none
integer :: j,vv(Dim_LM),k,v1,v2,v3,v4,v5,v6,v7,v8,v9
j=0
if(Dim_LM>9) stop 'Spatial_Dim>9, check permutation subroutine'
do v1=0,Vmax(1)
    vv(1)=v1
    if(Dim_LM>= 2) then
        do v2=0,min(Vtot-vv(1),Vmax(2))
           vv(2)=v2
           if(Dim_LM>= 3) then
               do v3=0,min(Vtot-sum(vv(1:2)),Vmax(3))
                   vv(3)=v3
                   if(Dim_LM>= 4) then
                       do v4=0,min(Vtot-sum(vv(1:3)),Vmax(4))
                           vv(4)=v4
                           if(Dim_LM>= 5) then
                               do v5=0,min(Vtot-sum(vv(1:4)),Vmax(5))
                                   vv(5)=v5
                                   if(Dim_LM>= 6) then
                                       do v6=0,min(Vtot-sum(vv(1:5)),Vmax(6))
                                           vv(6)=v6
                                           if(Dim_LM>= 7) then
                                               do v7=0,min(Vtot-sum(vv(1:6)),Vmax(7))
                                                   vv(7)=v7
                                                   if(Dim_LM>= 8) then
                                                       do v8=0,min(Vtot-sum(vv(1:7)),Vmax(8))
                                                           vv(8)=v8
                                                           if(Dim_LM== 9) then
                                                               do v9=0,min(Vtot-sum(vv(1:8)),Vmax(9))
                                                                   vv(9)=v9
                                                                   j=j+1
                                                               end do
                                                           else
                                                               j=j+1
                                                           end if
                                                       end do
                                                   else
                                                       j=j+1
                                           end if
                                       end do
                                   else
                                       j=j+1
                                   end if
                               end do
                           else
                               j=j+1
                            end if
                       end do
                   else 
                       j=j+1 
                    end if
               end do
           else
               j=j+1
           end if
        end do
    else
        j=j+1
    end if
end do
else
    j=j+1
end if
end do
Jmax = j                
write(*,*) 'Jmax = ', Jmax
ALLOCATE(v(Dim_LM,Jmax))
j=0
if(Dim_LM>9) stop 'Spatial_Dim>9'
do v1=0,Vmax(1)
vv(1)=v1
if(Dim_LM>= 2) then
  do v2=0,min(Vtot-vv(1),Vmax(2))
     vv(2)=v2
     if(Dim_LM>= 3) then
        do v3=0,min(Vtot-sum(vv(1:2)),Vmax(3))
           vv(3)=v3
           if(Dim_LM>= 4) then
              do v4=0,min(Vtot-sum(vv(1:3)),Vmax(4))
                 vv(4)=v4
                 if(Dim_LM>= 5) then
                    do v5=0,min(Vtot-sum(vv(1:4)),Vmax(5))
                       vv(5)=v5
                       if(Dim_LM>= 6) then
                          do v6=0,min(Vtot-sum(vv(1:5)),Vmax(6))
                             vv(6)=v6
                             if(Dim_LM>= 7) then
                                do v7=0,min(Vtot-sum(vv(1:6)),Vmax(7))
                                   vv(7)=v7
                                   if(Dim_LM>= 8) then
                                      do v8=0,min(Vtot-sum(vv(1:7)),Vmax(8))
                                         vv(8)=v8
                                         if(Dim_LM== 9) then
                                            do v9=0,min(Vtot-sum(vv(1:8)),Vmax(9))
                                               vv(9)=v9
                                               j=j+1
                                               v(:,j)=vv 
                                            end do
                                         else
                                            j=j+1
                                            v(:,j)=vv           
                                         end if
                                      end do
                                   else
                                      j=j+1
                                      v(:,j)=vv           
                                   end if
                                end do
                             else
                                j=j+1
                                v(:,j)=vv           
                             end if
                          end do
                       else
                          j=j+1
                          v(:,j)=vv                            
                       end if
                    end do
                 else 
                    j=j+1 
                    v(:,j)=vv
                  end if
              end do
           else
              j=j+1
              v(:,j)=vv           
           end if
        end do
     else
        j=j+1
        v(:,j)=vv           
     end if
  end do
else
  j=j+1
  v(:,j)=vv           
end if
end do
end subroutine permutation
!==============================================================================!
!==============================================================================!
!==============================================================================!
subroutine quad_root(j,m,n)
!==============================================================================!
!    Discussion:
!       quad_root uses the quadratic formula for explicity setting the thread
!       allocations for the parallel region
!    Arguments:
!    Output:
!==============================================================================!
integer :: j,m,n
m=ceiling(0.5*(sqrt(1.+8*j)-1))
n=j-(m*(m-1))/2
end subroutine quad_root

!==============================================================================!
!==============================================================================!
!==============================================================================!
end module Lmon_water_qmc
!==============================================================================!
!==============================================================================!
!==============================================================================!
program Single_monomer_in_cluster
use Lmon_water_qmc
use omp_lib
!==============================================================================!
implicit none
integer(kind=8) :: skip
integer:: Nsobol, g, i, j, j2(1), k, n, m, data_freq, Jmax2, jstart, jfinish, &
          deltaj
double precision, allocatable:: omega(:), Hess(:,:), U(:,:), q(:), q0(:), &
                                q1(:), force(:), y(:), herm(:,:), A(:,:,:), &
                                Umat(:), U1mat(:), C(:,:), eigenvalue(:)
double precision:: freq_cutoff, freq_replace, E, V0, potdif, pot, B, E0(1), &
                   E1(1), E2(1), E3(1)
character(len=50) :: coord_in
!==============================================================================!
!                         Variables to run dsyev                               !
!==============================================================================!
integer :: info, lwork 
double precision, allocatable :: work(:)
!==============================================================================!
!                         Variables to run ssobol.f                            !
!==============================================================================!
integer :: TAUS, IFLAG, max_seq, SAM
logical :: FLAG(2)
!==============================================================================!
!                              Read Input Files                                !
!==============================================================================!
read(*,*) Dim_LM             !VM: changed
read(*,*) Vtot
allocate(Vmax(Dim_LM))
read(*,*) Vmax(1:Dim_LM)
read(*,*) potential                
!read(*,*) freq_cutoff
!read(*,*) freq_replace
read(*,*) Nsobol
read(*,*) skip
read(*,*) data_freq
read(*,*) Nthreads
read(*,*) coord_in
write(*,*) 'Test 1; Successfully Read Input File'
!freq_cutoff = freq_cutoff / autocm
!freq_replace = freq_replace / autocm
!==============================================================================!
!                           Variables to run ssobol.f                          !
!==============================================================================!
SAM = 1
max_seq = 30
IFLAG = 1
!==============================================================================!
!                                 Read xyz File                                !
!==============================================================================!
OPEN(61,File=coord_in)
read(61,*) Natoms
read(61,*)
Dim= 3*Natoms ! cartesian coordinates
write(*,*) 'Test 2; Successfully Read Coordinates!'
!==============================================================================!
ALLOCATE(atom_type(Natoms), sqrt_mass(Dim), mass(Natoms), q(Dim), q0(Dim), &
        q1(Dim), force(Dim), Hess(Dim1,Dim1), omega(Dim1), U(Dim1,Dim1))   
!==============================================================================!
!                           Read Atom Type, Get Mass                           !
!            q0 contains the xyz coordinates of the input configuration        !
!==============================================================================!
! Input coordinates are in Angstroms
do i=1,Natoms
    read(61,*) atom_type(i), q0(3*i-2:3*i)   
    mass(i)=Atom_mass(atom_type(i))
    sqrt_mass(3*i-2:3*i)=SQRT(mass(i))
end do
close(61)
! Convert coordinates to Atomic Units
q0=q0/bohr            
write(*,*) 'Test 3; Successfully Read Coordinates'
!==============================================================================!
!                        Input Configuration Energy                            !
!==============================================================================!
CALL water_potential(Natoms/3, q0, E, force)
write(*,*) 'E0 =',E*autocm, 'cm-1'
!write(*,*) 'Test 5; Successfully Evaluate Input Geometry Energy'
!==============================================================================!
!                     Compute Hessian and Frequencies                          !
!==============================================================================!
CALL Get_Hessian(q0,Hess)
CALL Frequencies_From_Hess(Dim1,Hess,omega,U)
write(*,*) 'Test 6; Successfully Compute Hessian and Frequencies'
!==============================================================================!
!                    U(Dim,Dim) Contains the Normal Modes
!       Replace small frequencies with a large value to improve convergence
!==============================================================================!
do k=1,Dim_LM     !VM: changed
!    if (omega(k)<freq_cutoff) omega(k) = freq_replace
   U(:,k) = U(:,k) / sqrt_mass(:) / sqrt(omega(k))    ! VM: changed
   !    do i = 1, Dim1
   !        U(i,k) = SQRT(1/omega(k)) / sqrt_mass(i)*U(i,k)   
   !    end do
end do
!write(*,*) 'Test 7; Successfully Compute Normal Modes'
!==============================================================================!
!                       Determine Jmax and Permutations
!==============================================================================!
CALL permutation()
Jmax2=(Jmax*(Jmax+1))/2
ALLOCATE(Umat(Jmax2), U1mat(Jmax2), C(Jmax,Jmax), eigenvalue(Jmax), y(Dim_LM), &
        herm(0:Vtot,Dim_LM), A(0:Vtot,0:Vtot,Dim_LM))
write(*,*) 'Test 8; Successfully Determine Jmax and Allocate Arrays '
!==============================================================================!
!                      Initialize for Potential Matrix                         !
!==============================================================================!
lwork = max(1,3*Jmax-1)
ALLOCATE(work(max(1,lwork)))
Umat=0d0  ! initialize PE matrix
U1mat=0d0 
work = 0d0
!OPEN(UNIT=81, FILE='eigenvalues.dat')
!OPEN(UNIT=82, FILE='eigenvectors.dat')
OPEN(UNIT=83, FILE='fundamentals.dat')
!==============================================================================!
!                       Begin Loop over Sobol Points                           !
!                      See ssobol.f for documentation                          !
!==============================================================================!
OPEN(UNIT=83, FILE='fundamentals.dat')
CALL INSSOBL(FLAG,Dim_LM,Nsobol,TAUS,y,max_seq,IFLAG)
write(*,*) 'Test 9; Successfully Define Sobol Scrambling Method'
do i = 1, Nsobol
    CALL GOSSOBL(y)
!==============================================================================!
!                   If Sobol point == 0 or == 1, replace                       !
!==============================================================================!
    do j=1,Dim_LM
        if(y(j)<.00000001)then
            y(j)=.00000001
        else if(y(j)>1-.00000001)then
            y(j) = 1-.00000001
        end if
    end do
!==============================================================================!
!           Compute Hermite Polynomials (see paper for normalization)
!            Factor of 1/(sqrt(2)) come from the normal distribution
!==============================================================================!
    CALL sobol_stdnormal(Dim_LM, y)
    y = y/SQRT(2.)    
    herm(0,:) = 1.0                       
    herm(1,:) = SQRT(2.)*y(:)       
    do j = 2,Vtot      
        herm(j,:) = (SQRT(2./j)*y(:)*herm(j-1,:)) - (SQRT(((j-1d0)/j))*herm(j-2,:))
    end do
!==============================================================================!
!                     Compute Hermite Polynomial Products 
!==============================================================================!
    do m=0,Vtot
        do n=0,Vtot
            A(m,n,:) = herm(m,:)*herm(n,:)
        end do
    end do
!==============================================================================!
!                     Transform to normal coordinates                          !
!==============================================================================!
    q1=q0     
    do j=1,Dim_LM
       q1(1:Dim1) = q1(1:Dim1) + y(j)*U(:,j)
    end do
!==============================================================================!
!                   Compute potential difference (PES-HA)                      !
!==============================================================================!
    CALL water_potential(Natoms/3, q1, pot, force)  
    potdif = pot - E
    do j=1,Dim_LM
       potdif = potdif - (0.5*omega(j)*y(j)**2)
    end do
!==============================================================================!
!                   Compute PE Difference Matrix (Umat)                        !
!                    U1mat is used for Partial Averages                        !
!                       Also Begin Parallel Region                             !
!==============================================================================!
    !$ call omp_set_num_threads(Nthreads)
    deltaj=Jmax2/Nthreads
    !$omp parallel &
    !$omp shared (U1mat, A, v) &
    !$omp private (B, k, j, m, n,jstart,jfinish, g)
    !$omp do
    do g=1,Nthreads
        jstart=(g-1)*deltaj
        jfinish=jstart+deltaj
        if (g==Nthreads) jfinish=Jmax2
        call quad_root(jstart,m,n)
        do j=jstart+1,jfinish
            if(n<m) then
                n=n+1        
            else
                n=1
                m=m+1
            end if
            B=potdif
            do k=1,Dim_LM        
                B=B*A(v(k,m),v(k,n),k)
            end do
            U1mat(j) = U1mat(j) + B      
        end do
    end do
    !$omp end do
    !$omp end parallel
!==============================================================================!
!                             End Parallel Region                              !
!                           Compute Partial Averages                           !
!==============================================================================!
    IF(MOD(i,data_freq)==0) THEN
        Umat = Umat + U1mat
        U1mat = 0d0
        m=1
        n=0
        do j=1,Jmax2     
            if(n<m) then
                n=n+1
            else
                n=1
                m=m+1
            end if
            C(n,m) = Umat(j)/i
            if(m==n) then
                C(n,n)=C(n,n)+sum(omega(1:Dim_LM)*(v(:,n)+0.5)) 
            else
                C(m,n) = C(n,m) 
            end if
        end do
!        write(*,*) 'C= ', C
        info = 0
        CALL dsyev('v', 'u', Jmax, C, Jmax, eigenvalue, work, Lwork, info)
        !        write(81,*) i, eigenvalue(:)    ! write eigenvalues to file
        !        write(82,*) i, C(:,:)           !  write eigenvectors to file
!==============================================================================!
!                       Compute Fundamental Frequencies                        !
!==============================================================================!
        E0 = 0d0
        E1 = 0d0
        E2 = 0d0
        E3 = 0d0
        do m = 1,Jmax
            if(sum(v(:,m))==0) then
                j2=maxloc(abs(C(m,:)))
                E0 = eigenvalue(j2)*autocm  ! ground state
             else if(sum(v(:,m))==1) then  ! first excited state
                if(v(1,m)==1) then  ! OH excited state
                   j2=maxloc(abs(C(m,:)))
                   E1 = eigenvalue(j2)*autocm              
                else if(v(2,m)==1) then   ! OH excited state
                   j2=maxloc(abs(C(m,:)))
                   E2 = eigenvalue(j2)*autocm              
                else if(v(3,m)==1) then    ! HOH excited state
                   j2=maxloc(abs(C(m,:)))
                   E3 = eigenvalue(j2)*autocm
                end if
             end if
        end do
        write(*,*) 'funds =', i,E0, E1 - E0, E2-E0, E3-E0
        write(83,*) i,E0, E1 - E0, E2-E0, E3-E0
        FLUSH(83)
    end IF
end do ! end loop over sobol points
close(UNIT=83)
close(UNIT=84)
write(*,*) 'Final Check Successful; Hello Universe!'
!==============================================================================!
!                               output.dat                                     !
!==============================================================================!
open(90,file='simulation.dat')
write(90,*) 'Here is the output file for your calculation'
write(90,*) 'Natoms=                     ', Natoms
write(90,*) 'Dim=                        ', Dim
write(90,*) 'Dim1=                        ', Dim1
write(90,*) 'Vmax =                      ', Vmax
write(90,*) 'Vtot =                      ', Vtot
write(90,*) 'Jmax =                      ', Jmax
write(90,*) 'potential=                  ', potential
write(90,*) 'Nsobol=                     ', Nsobol
write(90,*) 'Nthreads=                     ', Nthreads
write(90,*) 'freq_cutoff=                     ', freq_cutoff*autocm
write(90,*) 'freq_replace=                     ', freq_replace*autocm
write(90,*) 'skip=                       ', skip
write(90,*) 'coord_in=                   ', coord_in
write(90,*) 'The HA ~ gave an energy of =  ', E
write(90,*) 'Fundamental Frequency Calculation!'
close(90)
end program Single_monomer_in_cluster
