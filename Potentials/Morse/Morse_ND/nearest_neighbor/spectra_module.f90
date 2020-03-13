!=================20==================40==================60==================80
!                            nD-Morse EigenSpectra
!==============================================================================!
!Compute the EigenSpectra for a given Potential
!example uses the nD morse potential, all equations written for general nD case
!Basis Functions:
!             phi(r):=(2*alpha_i/pi)^(d/4)*exp[-alpha_i(r-r_i)^2]
!Assumes a QRG; generates alpha (Gaussian Widths) based on nearest neghbor
!Needs gen_hermite_rule.f90 code (Gauss-Hermite-Quadriture) for potential eval.
!Needs LLAPACK to solve the generalized eigenvalue problem
!==============================================================================!
!       Modified:
!   13 March 2020
!       Author:
!   Shane Flynn
!==============================================================================!
!       To Do:
!Add in alternate alpha constructions to accomidate different grids
!Add in alternate basis functions for special case of d=2
!Make seperate subroutines to compute potential elements and hamiltonian
!==============================================================================!
module spectra_mod
implicit none
!==============================================================================!
contains
!==============================================================================!
subroutine read_grid(d,NG,x)
!==============================================================================!
!Reads in grid points, each line contains coordinates for 1 d-dimesional point
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG               ==>Number of grid points
!x(d,NG)          ==>All grid points coordinates (x={x^i})
!==============================================================================!
implicit none
integer::d,NG,i
double precision::x(d,NG)
open(17,File='grid.dat')
do i=1,NG
  read(17,*) x(:,i)
enddo
close(17)
end subroutine read_grid
!==============================================================================!
subroutine gaussian_widths(d,NG,x,alpha_0,alpha)
!==============================================================================!
!Generate scaling parameters for each grid point (proportional to widths)
!Uses nearest neighbor algorithm, good for uniform or QRG
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG               ==>Number of grid points
!x(d,NG)          ==>All grid points coordinates (x={x^i})
!alpha0           ==>Flat scaling parameter for optimizing widths
!alpha(NG)        ==>Gaussian width parameter
!==============================================================================!
implicit none
integer::d,NG,i,j
double precision::alpha_0,alpha(NG),r2,x(d,NG)
do i=1,NG
  alpha(i)=1d20                          !large distance for initial placeholder
  do j=1,NG
    if(j.ne.i) then
      r2=sum((x(:,i)-x(:,j))**2)                    !distance between gridpoints
      if(r2<alpha(i)) alpha(i)=r2                       !locate nearest neighbor
    endif
  enddo
  alpha(i)=alpha_0/alpha(i)
enddo
open(unit=18,file='alphas.dat')
do i=1,NG
  write(18,*) alpha(i)
enddo
close(18)
end subroutine gaussian_widths
!==============================================================================!
subroutine overlap_elements(d,x_i,x_j,alpha_i,alpha_j,S_ij)
!==============================================================================!
!Compute i-jth element of the overlap matrix
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!alpha_i          ==>i-th grid points gaussian width parameters
!Sij              ==>i-j element of the overlap matrix
!==============================================================================!
integer::d,i,j
double precision::alpha_i,alpha_j,S_ij,x_i(d),x_j(d),aij,r2
aij=alpha_i*alpha_j/(alpha_i+alpha_j)
r2=sum((x_i(:)-x_j(:))**2)
S_ij=(2*sqrt(alpha_i*alpha_j)/(alpha_i+alpha_j))**(0.5*d)*exp(-aij*r2)
end subroutine overlap_elements
!==============================================================================!
subroutine overlap_matrix(d,NG,x,alpha,Smat)
!==============================================================================!
!Compute the overlap matrix (symmetric)
!Used to check the stability of the basis; see diagonalize_overlap subroutine
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG               ==>Number of grid points
!x(d,NG)          ==>All grid points coordinates (x={x^i})
!alpha(NG)        ==>Gaussian width parameter
!Smat(NG,NG)      ==>Overlap Matrix
!==============================================================================!
integer::d,NG,i,j
double precision::alpha(NG),Smat(NG,NG),x(d,NG),aij
do i=1,NG
  do j=i,NG
    call overlap_elements(d,x(:,i),x(:,j),alpha(i),alpha(j),Smat(i,j))
    Smat(j,i)=Smat(i,j)
  enddo
enddo
end subroutine overlap_matrix
!==============================================================================!
subroutine overlap_eigenvalues(NG,Smat,eigenvalues,Lwork)
!==============================================================================!
!Compute the Eigenvalues for the overlap matrix to check numerical stability
!Must be positive definite; also check the Recriprocal Condition Number
!Lwork, work, and info defined according to LLAPACK suggestions
!==============================================================================!
!NG               ==>Number of grid points
!Smat(NG,NG)      ==>Overlap Matrix
!eigenvalues(NG)  ==>Matrix eigenvalues
!==============================================================================!
implicit none
integer::i,NG,info,Lwork
double precision::eigenvalues(NG),Smat(NG,NG),work(max(1,Lwork))
call dsyev('v','u',NG,Smat,NG,eigenvalues,work,Lwork,info)              !LLAPACK
write(*,*) 'Info (Overlap Matrix) ==> ', info
write(*,*) 'RCN ==> ', eigenvalues(1)/eigenvalues(NG)
open(unit=19,file='overlap.dat')
do i=1,NG
  write(19,*) eigenvalues(i)
enddo
close(19)
end subroutine overlap_eigenvalues
!==============================================================================!
subroutine kinetic_elements(d,x_i,x_j,alpha_i,alpha_j,T_ij)
!==============================================================================!
!Compute i-jth element of the kinetic matrix
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!alpha_i          ==>i-th grid points gaussian width parameters
!T_ij             ==>i-j element of the kinetic matrix
!==============================================================================!
implicit none
integer::d
double precision::x_i(d),x_j(d),alpha_i,alpha_j,aij,r2,T_ij
aij=alpha_i*alpha_j/(alpha_i+alpha_j)
r2=sum((x_i(:)-x_j(:))**2)
T_ij=aij*(d-2.*aij*r2)
end subroutine kinetic_elements
!==============================================================================!
function V_i(d,x_i)
!==============================================================================!
!Potential Energy
!Equations work for any d-dimesional system, define the potential accordingly
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!D_morse          ==>Parameter for Morse Potential
!omega(d)         ==>Parameter for Morse Potential
!V                ==>evaluate V(x_i)
!==============================================================================!
implicit none
integer::d
double precision::x_i(d),V_i
double precision,parameter::omega(3)=(/0.2041241,0.18371169,0.16329928/)
double precision,parameter::D_morse=12.
V_i=D_morse*sum((exp(-omega(:)*x_i(:))-1.)**2 )        !Sum of 1D Morse potentials
end function V_i
!==============================================================================!
subroutine get_hamiltonian(d,NG,x,alpha,Smat,GH_order,Hmat)
!==============================================================================!
!Compute the Hamiltonian Matrix
!Use Gauss-Hermite quadrature for the potential
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG               ==>Number of grid points
!x(d,NG)          ==>All grid points coordinates (x={x^i})
!alpha(NG)        ==>Gaussian width parameter
!Smat(NG,NG)      ==>Overlap Matrix
!GH_order         ==>Gauss-Hermite quadrature order
!Hmat(NG.NG)      ==>Hamiltonian Matrix
!z(GH_order)      ==>Quadrature points
!w(GH_order)      ==>Quadrature points weights
!==============================================================================!
implicit none
integer::d,NG,GH_order,i,j,k,ll
double precision,parameter::pi=4.*atan(1d0)
double precision::x(d,NG),alpha(NG),Smat(NG,NG),Hmat(NG,NG),z(GH_order)
double precision::w(GH_order),x_ij(d),Vij,l(d),rr(d),r2,eigenvalues(NG)
call cgqf(GH_order,6,0d0,0d0,0d0,1d0,z,w)                               !LLAPACK
w=w/sqrt(pi)
do i=1,NG
  do j=i,NG
    call overlap_elements(d,x(:,i),x(:,j),alpha(i),alpha(j),Smat(i,j))
    call kinetic_elements(d,x(:,i),x(:,j),alpha(i),alpha(j),Hmat(i,j))
!==============================================================================!
!                           Potential Energy
!==============================================================================!
    x_ij(:)=(alpha(i)*x(:,i)+alpha(j)*x(:,j))/(alpha(i)+alpha(j))
    Vij=0d0
    l(:)=1
    do ll=1,GH_order**d
      do k=1,d
        rr(k)=z(l(k))
      enddo
      rr=x_ij+rr/sqrt(alpha(i)+alpha(j))
      r2=V_i(d,rr)
      do k=1,d
        r2=r2*w(l(k))
      enddo
      Vij=Vij+r2
      do k=1,d
        l(k)=mod(l(k),float(GH_order))+1
        if(l(k).ne.1) exit
      enddo
    enddo
!==============================================================================!
!      Hamiltonian = Kinetic + Potential (all terms are scaled by Overlap)
!==============================================================================!
    Hmat(i,j)=(Hmat(i,j)+Vij)*Smat(i,j)
    Hmat(j,i)=Hmat(i,j)
  enddo
enddo
end subroutine get_hamiltonian
!==============================================================================!
subroutine hamiltonian_eigenvalues(NG,Smat,Hmat,eigenvalues,Lwork)
!==============================================================================!
!Compute the Eigenvalues for the Hamiltonian matrix
!Lwork, work, and info defined according to LLAPACK suggestions
!==============================================================================!
!NG               ==>Number of grid points
!Smat(NG,NG)      ==>Overlap Matrix
!Hmat(NG.NG)      ==>Hamiltonian Matrix
!eigenvalues(NG)  ==>Matrix eigenvalues
!==============================================================================!
implicit none
integer::NG,info,Lwork,itype,i
double precision::Smat(NG,NG),Hmat(NG,NG),eigenvalues(NG),work(max(1,lwork))
itype=1
call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
open(unit=20,file='eigenvalues.dat')
do i=1,NG
  write(20,*) eigenvalues(i)
enddo
close(20)
end subroutine hamiltonian_eigenvalues
!==============================================================================!
subroutine write_out(d,NG,alpha0,GH_order)
!==============================================================================!
!Simulation details
!==============================================================================!
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!NG               ==>Number of grid points
!alpha0           ==>Flat scaling parameter for optimizing widths
!GH_order         ==>Gauss-Hermite quadrature order
!==============================================================================!
implicit none
integer::d,NG,GH_order
double precision::alpha0
open(unit=99,file='out')
write(99,*) 'dimensionality ==> ', d
write(99,*) 'Number of  Gaussians==> ', NG
write(99,*) 'alpha0==>', alpha0
write(99,*) 'GH Order==>', GH_order
close(99)
end subroutine write_out
!==============================================================================!
end module spectra_mod
