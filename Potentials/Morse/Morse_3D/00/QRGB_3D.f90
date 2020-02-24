!=============================================================================80
!                      3D Morse Code with qLJ Basis
!==============================================================================!
!       Discussion:
!Quasi-Regular Gaussian Basis
!using optimized grid solve generalized eigenvalue problem
!compute everything here except for generating the grid
!==============================================================================!
!       Modified:
!   15 May 2019
!       Author:
!   Shane Flynn
!==============================================================================!
module QRGB_mod
implicit none
!==============================================================================!
!                            Global Variables
!==============================================================================!
!d              ==> Gaussian dsionality
!==============================================================================!
integer::d
double precision::E_cut,c_LJ,integral_P
!==============================================================================!
contains
!==============================================================================!
function V(x)
!==============================================================================!
!       Discussion:
!Hard-coded Morse Potential Energy
!==============================================================================!
implicit none
double precision::x(d),V
double precision,parameter::omega(3)=(/0.2041241,0.18371169,0.16329928/)
double precision,parameter::D_morse=12.
V=D_morse*sum((exp(-omega(:)*x(:))-1.)**2 )
end function V
!==============================================================================!
function P(x)
!==============================================================================!
!       Discussion:
!Target Distribution Function
!==============================================================================!
implicit none
double precision::x(d),P
if(V(x)<E_cut) then
   P=(E_cut-V(x))**(d/2.)/integral_P
else                                              !set equal to 0 if beyond Ecut
   P=1d-20
end if
end function P
!==============================================================================!
end module QRGB_mod
!==============================================================================!
program main
use QRGB_mod
!==============================================================================!
!       Discussion:
!==============================================================================!
!d              ==> i-th gaussian dsionality (x^i=x^i_1,x^i_2,..,x^i_d)
!NG             ==> Number of basis functions
!x              ==>(d) all atom coordinates
!==============================================================================!
implicit none
logical::unif_grid
character(len=50)::grid_in,theory_in
integer::NG,GH_order,i,j,l1,l2,l3
double precision,parameter::pi=4.*atan(1d0)
double precision::aij,r2,alpha0,Vij,RCN
double precision,allocatable,dimension(:)::alpha,eigenvalues,x_ij,z,w,rr,theory
double precision,allocatable,dimension(:,:)::x,Smat,Hmat
!==============================================================================!
!                       LLAPACK dsygv variables                                !
!==============================================================================!
integer::itype,info,lwork
double precision,allocatable,dimension(:)::work
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
read(*,*) d
read(*,*) NG
read(*,*) GH_order
read(*,*) grid_in
read(*,*) integral_P
read(*,*) theory_in
read(*,*) unif_grid
read(*,*) E_cut
read(*,*) c_LJ
read(*,*) alpha0
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,NG),x_ij(d),rr(d),alpha(NG),eigenvalues(NG),Smat(NG,NG))
allocate(Hmat(NG,NG),theory(NG),z(GH_order),w(GH_order))
write(*,*) 'Test 0; Successfully Read Input File'
!==============================================================================!
!                           Read GridPoints x(d,NG)
!==============================================================================!
open(90,File=grid_in)
do i=1,NG
    read(90,*) x(:,i)
enddo
close(90)
!==============================================================================!
!                       Generate Alphas alpha(NG)
!==============================================================================!
if(unif_grid.EQV..TRUE.)then
    alpha=alpha0
    write(*,*) 'Test 1; Uniform Grid, Alpha:=Constant'
else
    do i=1,NG
        alpha(i)=0.5*(alpha0/(c_LJ*(P(x(:,i))*NG)**(-1./d))**2)
    enddo
    write(*,*) 'Test 1; Successfully Generated Alphas for QLJ Grid'
endif
!==============================================================================!
!                          Write Alphas to File
!==============================================================================!
open(unit=91,file='alphas.dat')
do i=1,NG
    write(91,*) alpha(i)
enddo
close(91)
!==============================================================================!
!                           Overlap Matrix (S)
!==============================================================================!
do i=1,NG
    do j=i,NG
         aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
         r2=sum((x(:,i)-x(:,j))**2)
         Smat(i,j)=(2*sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d)&
             *exp(-aij*r2)
       Smat(j,i)=Smat(i,j)
    enddo
enddo
!==============================================================================!
!                   Check to see if S is positive definite
!If this is removed, you need to allocate llapack arrays before Hamiltonian
!==============================================================================!
lwork=max(1,3*NG-1)
allocate(work(max(1,lwork)))
call dsyev('v','u',NG,Smat,NG,eigenvalues,work,Lwork,info)
write(*,*) 'Info (Initial Overlap Matrix) ==> ', info
RCN = eigenvalues(1)/eigenvalues(NG)
write(*,*) 'RCN =', RCN
open(unit=92,file='overlap_eigenvalues.dat')
do i=1,NG
    write(92,*) eigenvalues(i)
enddo
close(19)
write(*,*) 'Test 2; Overlap Matrix is Positive Definite'
!==============================================================================!
!             Use Gauss Hermit quadrature to evaluate the potential matrix
! z(GH-order) w(GH-order) --- quadrature points and weights
!==============================================================================!
call cgqf(GH_order,6,0d0,0d0,0d0,1d0,z,w)
w=w/sqrt(pi)
!==============================================================================!
!                   Solve Generalized Eigenvalue Problem
!==============================================================================!
do i=1,NG
  do j=i,NG
     aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
     r2=sum((x(:,i)-x(:,j))**2)
     Smat(i,j)=(2*sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d)&
         *exp(-aij*r2)
     Smat(j,i)=Smat(i,j)
     ! kinetic energy:
     Hmat(i,j)=aij*(d-2*aij*r2)
     ! potential energy Vij matrix element
     x_ij(:)=(alpha(i)*x(:,i)+alpha(j)*x(:,j))/(alpha(i)+alpha(j))
     Vij=0d0
           do l3=1,GH_order
        do l2=1,GH_order
     do l1=1,GH_order
              rr(1)=z(l1)
              rr(2)=z(l2)
              rr(3)=z(l3)
              rr=x_ij+rr/sqrt(alpha(i)+alpha(j))
              Vij=Vij+w(l1)*w(l2)*w(l3)*V(rr)
           enddo
        enddo
     enddo
  ! kinetic + potential energy
     Hmat(i,j)=(Hmat(i,j)+Vij)*Smat(i,j)
     Hmat(j,i)=Hmat(i,j)
  enddo
enddo
itype=1
eigenvalues=0d0
call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
write(*,*) 'info ==> ', info
open(unit=21,file='eigenvalues.dat')
write(21,*) alpha0, eigenvalues(:)
close(21)
!==============================================================================!
!                              Exact Eigenvalues
!==============================================================================!
open(25,File=theory_in)
do i=1,NG
    read(25,*) theory(i)
enddo
close(25)
open(unit=73,file='abs_error.dat')
open(unit=74,file='rel_error.dat')
open(unit=75,file='alpha_abs_error.dat')
open(unit=76,file='alpha_rel_error.dat')
do i=1,NG
        write(73,*) i, abs(theory(i)-eigenvalues(i))
        write(74,*) i, (eigenvalues(i)-theory(i))/theory(i)
        write(75,*) alpha0, abs(theory(i)-eigenvalues(i))
        write(76,*) alpha0, (eigenvalues(i)-theory(i))/theory(i)
enddo
close(73)
close(74)
close(75)
close(76)
open(unit=78,file='alpha_rel_124.dat')
do i=1,124
        write(78,*) alpha0, (eigenvalues(i)-theory(i))/theory(i)
enddo
close(78)
!==============================================================================!
!                               output file                                    !
!==============================================================================!
open(90,file='simulation.dat')
write(90,*) 'dimensionality ==> ', d
write(90,*) 'NG ==> ', NG
write(90,*) 'Integral_P==>', integral_P
write(90,*) 'c LJ==>', c_LJ
write(90,*) 'Ecut ==>', E_cut
write(90,*) 'alpha0==>', alpha0
write(90,*) 'RCN Overlap Matrix==>', RCN
write(90,*) 'GH Order==>', GH_order
close(90)
write(*,*) 'Hello Universe!'
end program main
