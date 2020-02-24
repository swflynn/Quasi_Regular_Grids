!=============================================================================80
!                      3D Morse with Quasi-Monte Carlo
!==============================================================================!
!Quasi-Regular Gaussian Basis
!Compute potential using quasi Monte Carlo (not GHQ)
!Requires the sobol.f90 module (quasi-random number generator)
!NOTE: this does not seem to be working (issue with equations), see 2D qMC
!for a working code using sobol integration for the potential
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
!d              ==>i-th gaussian dsionality (x^i=x^i_1,x^i_2,..,x^i_d)
!NG             ==>Number of Gaussian Basis Functions
!integral_P     ==>Normalization factor for P(x) (from Grid Generation)
!E_cut          ==>Energy Cutoff Contour    (be consistent with Grid Generation)
!c_LJ           ==>Lennard-Jones parameter  (be consistent with Grid Generation)
!min/max        ==>Domain for 3D Morse Potential
!==============================================================================!
integer::d,NG
double precision::E_cut,c_LJ,integral_P,xmin,xmax,ymin,ymax,zmin,zmax
!==============================================================================!
contains
!==============================================================================!
function Potential(x_i)
!==============================================================================!
!Hard-coded Morse Potential Energy
!==============================================================================!
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!Potential      ==>evaluate Potential(x)
!D_morse        ==>Parameter for Morse Potential
!c_xyz          ==>(d) Parameter for Morse Potential
!==============================================================================!
implicit none
double precision,parameter::D_morse=12.
double precision,parameter::c_x=0.2041241
double precision,parameter::c_y=0.18371169
double precision,parameter::c_z=0.16329928
double precision::x_i(d),Potential
Potential=D_morse*((exp(-c_x*x_i(1))-1)**2+(exp(-c_y*x_i(2))-1)**2+&
    (exp(-c_z*x_i(3))-1)**2)
end function Potential
!==============================================================================!
function P_x(x_i)
!==============================================================================!
!Target Distribution Function
!==============================================================================!
!P_x            ==>evaluate P(x)
!x_i            ==>(d) ith particles coordinate x^i_1,..,x^i_d
!Del_par        ==>Delta parameter for P(x) distribution, :=10% of Ecut
!==============================================================================!
implicit none
double precision::Del_par,x_i(d),P_x
Del_par=0.01*E_cut
if(Potential(x_i)<E_cut) then
   P_x=(E_cut+Del_par-Potential(x_i))**(d/2.)/integral_P
else                                              !set equal to 0 if beyond Ecut
   P_x=1d-20
end if
end function P_x
!==============================================================================!
subroutine normalize_P(N_Eval)
!==============================================================================!
!Normalize the target distribution function
!For Morse: integrating over [a,b],[a,b],[a,b]
!int P(r)~Area_Square/N sum_n=1,N P(r_n)
!==============================================================================!
!norm           ==>evaluate P(r)
!r              ==>(d) coordinates for evaluating P(r)
!N_eval         ==>number of evaluations for integral approximation
!a,b            ==>bounds for square to normalize over (b-a)^2 = area
!==============================================================================!
integer::N_eval,counter
double precision::r(d),norm
integral_P=1d0                     !set equal to 1 so you can initially call P_x
norm=0d0                                              !compute the normalization
counter=0
do while(counter.lt.N_Eval)                   !only want points within the E_cut
    call random_number(r)
    r(1)=xmin+r(1)*(xmax-xmin)
    r(2)=ymin+r(2)*(ymax-ymin)
    r(3)=zmin+r(3)*(zmax-zmin)
    if(Potential(r)<E_cut)then
        norm=norm+P_x(r)
        counter=counter+1
    endif
enddo
norm=norm*(xmax-xmin)*(ymax-ymin)*(zmax-zmin)/N_eval
integral_P=norm
end subroutine normalize_P
!==============================================================================!
end module QRGB_mod
!==============================================================================!
!==============================================================================!
program main
use QRGB_mod
!==============================================================================!
!==============================================================================!
!x              ==>(d) all atom coordinates
!x0             ==>(d) store previous coordinate before trial move
!freq           ==>Interval to update mv_cutoff size
!accept         ==>number of accepted trial moves, for acceptance~50%
!counter        ==>total number of moves, for acceptance~50%
!x2 ith gaussians coordinates for matrix elements
!z integration sequence form matrix elements
!x_ij= i-jth gaussian center
!z=sequence for integration
!t_i,t_f        ==>cpu time to ~ simulation time
!==============================================================================!
implicit none
character(len=50)::grid_in,theory_in
integer::Nsobol,data_freq,n,i,j,counter,l,N_Px
integer*8::skip
double precision::time1,time2,aij,r2,alpha0
double precision,allocatable,dimension(:)::alpha,eigenvalues,x_ij,theory
double precision,allocatable,dimension(:,:)::x,Smat,z,Vmat,Hmat
!==============================================================================!
!                       LLAPACK dsygv variables                                !
!==============================================================================!
integer::itype,info,lwork
double precision,allocatable,dimension(:)::work
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
call cpu_time(time1)
read(*,*) d
read(*,*) NG
read(*,*) Nsobol
read(*,*) N_Px
read(*,*) alpha0
read(*,*) E_cut
read(*,*) c_LJ
read(*,*) grid_in
read(*,*) theory_in
read(*,*) data_freq
read(*,*) xmin
read(*,*) xmax
read(*,*) ymin
read(*,*) ymax
read(*,*) zmin
read(*,*) zmax

Nsobol=Nsobol/data_freq
Nsobol=Nsobol*data_freq
write(*,*) 'data_freq=',data_freq, '  Nsobol=',Nsobol
skip=Nsobol
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,NG),x_ij(d),alpha(NG),eigenvalues(NG),Smat(NG,NG),Vmat(NG,NG))
allocate(z(d,Nsobol),Hmat(NG,NG),theory(NG))
write(*,*) 'Test 0; Successfully Read Input File'
!==============================================================================!
!                   Determine Normalization constant for P_x
!==============================================================================!
call normalize_P(N_Px)
!==============================================================================!
!                           Input GridPoints x(d,NG)
!==============================================================================!
open(16,File=grid_in)
do n=1,NG
    read(16,*) x(:,n)
enddo
close(16)
!==============================================================================!
!                          Generate Alpha Scaling
!==============================================================================!
open(unit=17,file='alphas.dat')
do i=1,NG
    alpha(i)=alpha0/(c_LJ*(P_x(x(:,i))*NG)**(-1./d))**2
    write(17,*) alpha(i)
enddo
close(17)
write(*,*) 'Test 3; Successfully Generated Gaussian Widths'
!==============================================================================!
!                           Overlap Matrix (S)
!==============================================================================!
do i=1,NG
    do j=i,NG
         aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
         r2=sum((x(:,i)-x(:,j))**2)
         Smat(i,j)=(sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d)&
             *exp(-0.5*aij*r2)
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
open(unit=19,file='overlap_eigenvalues.dat')
do i=1,NG
    write(19,*) eigenvalues(i)
enddo
close(19)
write(*,*) 'Test 4; Overlap Matrix is Positive Definite'
!==============================================================================!
!                   Generate Sequence For Evaluating Potential
!==============================================================================!
do l=1,Nsobol
    call sobol_stdnormal(d,skip,z(:,l))
enddo
write(*,*) 'Test 5; Successfully Generated Integration Sequence'
!==============================================================================!
!                              Theory Eigenvalues Value
!==============================================================================!
open(25,File=theory_in)
do i=1,NG
    read(25,*) theory(i)
enddo
close(25)
!==============================================================================!
!                             Evaluate Potential
!==============================================================================!
Vmat=0d0
open(unit=21,file='eigenvalues.dat')
do counter=1,Nsobol/data_freq
   do i=1,NG
      do j=i,NG
         x_ij(:)=(alpha(i)*x(:,i)+alpha(j)*x(:,j))/(alpha(i)+alpha(j))
         do l=(counter-1)*data_freq+1,counter*data_freq
            Vmat(i,j)=Vmat(i,j)+Potential(x_ij(:)+z(:,l)/sqrt(alpha(i)+alpha(j)))
         enddo
      enddo
   enddo
!==============================================================================!
!                   Solve Generalized Eigenvalue Problem
!==============================================================================!
   write(*,*) 'Generalized Eigenvalue Iteration', counter*data_freq
   do i=1,NG
      do j=i,NG
         aij=alpha(i)*alpha(j)/(alpha(i)+alpha(j))
         r2=sum((x(:,i)-x(:,j))**2)
         Smat(i,j)=(sqrt(alpha(i)*alpha(j))/(alpha(i)+alpha(j)))**(0.5*d)&
             *exp(-0.5*aij*r2)
         Smat(j,i)=Smat(i,j)
         ! kinetic energy:
         Hmat(i,j)=0.5*aij*(d-aij*r2)
         ! kinetic + potential energy
         Hmat(i,j)=(Hmat(i,j)+Vmat(i,j)/(counter*data_freq))*Smat(i,j)
         Hmat(j,i)=Hmat(i,j)
      enddo
   enddo
   itype=1
   eigenvalues=0d0
   call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
   write(21,*) alpha0, eigenvalues(:)
   write(*,*) 'info ==> ', info
enddo
close(21)
close(73)
open(unit=73,file='abs_error.dat')
open(unit=74,file='rel_error.dat')
do i=1,NG
    write(74,*) i, (eigenvalues(i)-theory(i))/theory(i)
    write(73,*) i, abs(eigenvalues(i)-theory(i))
enddo
close(73)
close(74)
!==============================================================================!
!                               output file                                    !
!==============================================================================!
call cpu_time(time2)
open(99,file='simulation.dat')
write(99,*) 'particle dimensionality ==> ', d
write(99,*) 'NG ==> ', NG
write(99,*) 'P_x Normalization Iterations ==> ', N_Px
write(99,*) 'Integral_P==>', integral_P
write(99,*) 'c LJ==>', c_LJ
write(99,*) 'Ecut ==>', E_cut
write(99,*) 'alpha0==>', alpha0
write(99,*) 'Nsobol==>', Nsobol
write(99,*) 'xmin ==>', xmin
write(99,*) 'xmax ==>', xmax
write(99,*) 'ymin ==>', ymin
write(99,*) 'ymax ==>', ymax
write(99,*) 'zmin ==>', zmin
write(99,*) 'zmax ==>', zmax
write(99,*) 'Total Time ==> ', time2-time1
close(99)
write(*,*) 'Hello Universe!'
end program main
