!=================20==================40==================60==================80
!                   2D Morse Quasi-Reglar Grid Generation
!==============================================================================!
!Generate gridpoints distributed via the 2D Morse Oscillator
!This Quasi-Regular gird is optimized using a quasi-Lennard Jones Potential
!Force minimization, accept any trial moves that reduces the system's energy
!See MMC grid for metropolis optimization (acceptance criteria)
!==============================================================================!
!Start developing 2D module want to make a module for each code then remove the
!various seperate directories
!remove the moment calculations, they are not useful
!the 2D case can use more efficient basis functions, these are not generalizable
!to the nD case, see the nD code for this generalization.
!==============================================================================!
!       Modified:
!   21 Feburary 2020
!       Author:
!   Shane Flynn
!==============================================================================!
module morse_grid_mod
implicit none
!==============================================================================!
!                            Global Variables
!==============================================================================!
!d              ==> Particle Dimensionality
!Npoints        ==> Number of grid points
!integral_P     ==> Normalization for P_x
!==============================================================================!
integer::Npoints
integer,parameter::d=2
double precision,allocatable,dimension(:)::xmin,xmax
double precision,parameter::D_morse=12.
double precision,parameter::omega(2)=(/0.2041241,0.18371169/)
double precision::integral_P
!==============================================================================!
contains
!==============================================================================!
function V(x)
!==============================================================================!
!Potential Energy (sum of 1D morse potentials)
!==============================================================================!
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!V              ==>evaluate V(x)
!D_morse        ==>Parameter for Morse Potential
!omega          ==>(d) Parameter for Morse Potential
!==============================================================================!
implicit none
double precision::x(d),V
V=D_morse*sum((exp(-omega(:)*x(:))-1.)**2)
end function V
!==============================================================================!
function P(x,E_cut)
!==============================================================================!
!Target Distribution Function
!==============================================================================!
!P              ==>evaluate P(x)
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!==============================================================================!
implicit none
double precision::x(d),E_cut,P
if(V(x)<E_cut) P=(E_cut-V(x))**(d/2.)/integral_P
if(V(x)>=E_cut) P=1d-20                           !set equal to 0 if beyond Ecut
end function P
!==============================================================================!
subroutine box_size_P(N_MMC_box,E_cut)
!==============================================================================!
!Determine the box size for normalizing P; (xmin,xmax) using MMC
!==============================================================================!
!N_MMC_box      ==>Number of MMC Iterations to determine box-size
!mv_cutoff      ==>trial displacement move cutoff
!r              ==>(d) coordinates
!r_trial        ==>(d) trial coordinates
!s              ==>(d) trail displacement; random number for coordinates
!xmin           ==>(d) minimum of normalization box
!xmax           ==>(d) maximum of normalization box
!==============================================================================!
integer::N_MMC_box,i,j
double precision::E_Cut,dummy,r_trial(d),r(d),s(d)
double precision,parameter::mv_cutoff=0.1
r=0d0
Integral_P=1d0                       !set equal to 1 so you can initially call P
xmin=r
xmax=r
do i=1,N_MMC_box
!==============================================================================!
!                   Generate coordinates for Random move
!           random numbers generated (0,1), make it (-1,1) ==> s=2*s-1
!==============================================================================!
  call random_number(s)
  r_trial=r+mv_cutoff*(2*s-1)
!==============================================================================!
!                             Test Acceptance
!==============================================================================!
  call random_number(dummy)
  if(P(r_trial,E_cut)/P(r,E_cut).ge.dummy) then
    r=r_trial
    do j=1,d
      if(xmin(j)>r(j)) xmin(j)=r(j)
      if(xmax(j)<r(j)) xmax(j)=r(j)
    enddo
  endif
enddo
do i=1,d
  write(*,*) 'Box Dimensions==>', xmin(i),xmax(i)
enddo
end subroutine box_size_P
!==============================================================================!
subroutine compute_integral_P(N_1D,E_cut)
!==============================================================================!
!Use a (uniform) square grid to integrate P(r)
!Box size for the grid is determined in the box_size subroutine
!int P(r)~Area_Square/N sum_n=1,N P(r_n)
!==============================================================================!
!N_1D           ==> Number of points in a single dimension for the integral
!Ntotal         ==>Total number of evaluations for all dimensions
!xmin           ==>(d) minimum of normalization box
!xmax           ==>(d) maximum of normalization box
!Moment         ==>First Moments for the distribution
!r(d)           ==>Coordinates
!==============================================================================!
integer::N_1D
integer::i,j,index1(d),Ntotal
double precision::r(d),Moment,dummy,delx(d),E_cut
Moment=0.
Ntotal=(N_1D+1)**d
index1=0
delx(:)=(xmax(:)-xmin(:))/N_1D
do i=1,Ntotal
  do j=1,d
    if(index1(j).eq.N_1D) then
      index1(j)=0
    else
      index1(j)=index1(j)+1
      exit
    endif
  enddo
  r(:)=xmin(:)+index1(:)*delx(:)
  dummy=P(r,E_cut)
  Moment=Moment+dummy
enddo
dummy=1./N_1D**d
do j=1,d
  dummy=dummy*(xmax(j)-xmin(j))
enddo
integral_P=dummy*Moment
write(*,*) 'integral_P ==> ', integral_P
end subroutine compute_integral_P
!==============================================================================!
subroutine initial_distribution(x,Npoints,E_cut)
!==============================================================================!
!Generate initial distribution of points, accept anything within E_cut contour
!==============================================================================!
implicit none
integer::Npoints,i,j
double precision::E_Cut,x(d,Npoints),s(d)
i=1
do while(i.le.Npoints)
  call random_number(s)
  s(:)=xmin(:)+s(:)*(xmax(:)-xmin(:))
  if(V(s)<E_cut)then
    x(:,i)=s(:)
    i=i+1
  endif
enddo
open(unit=17,file='coor_ini.dat')
do i=1,Npoints
  write(17,*) x(:,i)
enddo
close(17)
end subroutine initial_distribution
!==============================================================================!
function Pair_LJ_NRG(x1,x2,Npoints,E_cut,c_LJ)
!==============================================================================!
!quasi-Lennard Jones pairwise energy between grid points
!==============================================================================!
!x1             ==>(d) ith atoms coordinates
!x2             ==>(d) jth atoms coordinates
!a/b            ==>evaluate q-LJ
!sigma          ==>Gaussian widths (c*sigma(P))
!Pair_LJ_NRG    ==>Energy of the i-j q-LJ potential
!c_LJ           ==>parameter for LJ
!==============================================================================!
implicit none
integer::Npoints
double precision::x1(d),x2(d),a,b,Pair_LJ_NRG,sigma1,sigma2,E_cut,c_LJ
a=sum((x1(:)-x2(:))**2)
sigma1=c_LJ*(P(x1,E_cut)*Npoints)**(-1./d)
sigma2=c_LJ*(P(x2,E_cut)*Npoints)**(-1./d)
b=(sigma2**2/a)
a=(sigma1**2/a)
Pair_LJ_NRG=a**(d+9)-a**(d+3)+b**(d+9)-b**(d+3)
end function Pair_LJ_NRG
!==============================================================================!
subroutine initial_pairwise_energy(U,x,Npoints,E_cut,c_LJ)
!==============================================================================!
!Compute the pairwise energies for all the initial Grid Points U[x_ij]
!==============================================================================!
implicit none
integer::Npoints,i,j
double precision::E_cut,c_LJ,U(Npoints,Npoints),x(d,Npoints),s(d)
do i=2,Npoints
  do j=1,i-1
    U(i,j)=Pair_LJ_NRG(x(:,i),x(:,j),Npoints,E_cut,c_LJ)
    U(j,i)=U(i,j)
  enddo
enddo
end subroutine initial_pairwise_energy
!!==============================================================================!





!!!!!!!!==============================================================================!
!!!!!!!function random_integer(Nmin,Nmax)
!!!!!!!!==============================================================================!
!!!!!!!!Randomly generate an integer in the range 1-Nparticles
!!!!!!!!==============================================================================!
!!!!!!!!Nmin           ==>minimum index value (1)
!!!!!!!!Nmax           ==>maximum index value (Nparticles)
!!!!!!!!random_integer ==>integer returned
!!!!!!!!a              ==>random number (0,1)
!!!!!!!!==============================================================================!
!!!!!!!implicit none
!!!!!!!integer::Nmin,Nmax,random_integer
!!!!!!!double precision::a
!!!!!!!call random_number(a)
!!!!!!!random_integer=floor(a*(Nmax-Nmin+1))+Nmin
!!!!!!!end function random_integer
!!!!!!!!==============================================================================!


!!==============================================================================!
!!                       Begin MMC to Optimize GridPoints
!!==============================================================================!
!accept=0
!counter=0
!mv_cutoff=0.01
!deltae1=0
!plt_count=0
!open(unit=18,file='mv_cut.dat')
!open(unit=19,file='delE.dat')
!do i=1,N_MC
!!==============================================================================!
!!                           Select Atom to Move
!!==============================================================================!
!    k=random_integer(1,Npoints)
!!==============================================================================!
!!                           Generate trial move
!!        random numbers generated (0,1), make it (-1,1) ==> s=2*s-1
!!==============================================================================!
!    call random_number(s)
!    x0=x(:,k)+mv_cutoff*(2*s-1)
!!==============================================================================!
!!                   Only consider point if V(trial)<Ecut
!!                 Compute Energy Change due to Trial Move
!!==============================================================================!
!    if(V(x0).lt.E_cut) then
!        counter=counter+1
!        U_move(k)=P(x0)
!        Delta_E=0d0
!        do j=1,Npoints
!            if(j.ne.k) then
!                U_move(j)=Pair_LJ_NRG(x(:,j),x0)
!                Delta_E=Delta_E+U(j,k)-U_move(j)
!            endif
!        enddo
!!==============================================================================!
!!               Accept any trial that decreases the energy
!!==============================================================================!
!        if(Delta_E.ge.0d0)then
!            U(:,k)=U_move(:)
!            U(k,:)=U_move(:)
!            accept=accept+1
!            x(:,k)=x0(:)
!            deltae1=deltae1+Delta_E
!        endif
!     endif
!!==============================================================================!
!!                          Update Cutoff Paramater
!!       acceptance ~50%, adjust trial movement displacement accordingly
!!==============================================================================!
!        if(mod(i,MMC_freq)==0)then
!            write(*,*) 'MMC Iteration', i
!            if(dble(accept)/counter<0.3)then
!                mv_cutoff=mv_cutoff*0.9
!            else
!                mv_cutoff=mv_cutoff*1.1
!            endif
!        accept=0
!        counter=0
!        call Moments(Moment,x)
!        write(*,*) 'MMC Moments:'
!        write(*,*) Moment(1:5)
!        write(*,*) 'mv cutoff', mv_cutoff
!        write(*,*) 'Deltae1==>', Deltae1
!        plt_count=plt_count+1
!        write(18,*) plt_count, mv_cutoff
!        write(19,*) plt_count, Deltae1
!        deltae1=0
!        endif
!enddo
!close(18)
!close(19)
!==============================================================================!
end module morse_grid_mod
!==============================================================================!
!==============================================================================!
program main
use morse_grid_mod
!==============================================================================!
!==============================================================================!
implicit none
integer::N_MMC_box,N_1D
double precision::E_cut,c_LJ
double precision,allocatable,dimension(:,:)::x,U
!!!!!integer::MMC_freq,accept,counter,i,j,k,plt_count
!!!!!double precision::Delta_E,deltae1,mv_cutoff,Moment(0:5)
!!!!!double precision,allocatable,dimension(:)::x0,s,U_move,xmin(:),xmax(:)
!!!!!double precision,allocatable,dimension(:,:)::x,U
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
read(*,*) N_MMC_box
read(*,*) E_cut
read(*,*) N_1D
read(*,*) Npoints
read(*,*) c_LJ
write(*,*) 'Test 0; Successfully Read Input File'
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(xmin(d),xmax(d),x(d,Npoints),U(Npoints,Npoints))
!1111allocate(U_move(Npoints))
write(*,*) 'Test 1; Successfully Allocated arrays'
!==============================================================================!
!develop call statements for subroutines to write
!==============================================================================!
call box_size_P(N_MMC_box,E_cut)
call compute_integral_P(N_1D,E_cut)
call initial_distribution(x,Npoints,E_cut)
call initial_pairwise_energy(U,x,Npoints,E_cut,c_LJ)
!call optimize_grid()
end program main
