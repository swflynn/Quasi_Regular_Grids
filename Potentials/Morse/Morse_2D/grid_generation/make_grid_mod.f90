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
!the 2D case can use more efficient basis functions, these are not generalizable
!==============================================================================!
!       Modified:
!   27 Feburary 2020
!       Author:
!   Shane Flynn
!==============================================================================!
module make_grid_mod
implicit none
!==============================================================================!
!                            Global Variables
!==============================================================================!
!d              ==> Particle Dimensionality
!Npoints        ==> Number of grid points
!integral_P     ==> Normalization for P_x
!==============================================================================!
integer::Npoints
double precision,allocatable,dimension(:)::xmin,xmax
double precision::integral_P
!==============================================================================!
contains
!==============================================================================!
function V(d,x,omega,D_morse)
!==============================================================================!
!Potential Energy (sum of 1D morse potentials)
!==============================================================================!
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!V              ==>evaluate V(x)
!D_morse        ==>Parameter for Morse Potential
!omega          ==>(d) Parameter for Morse Potential
!==============================================================================!
implicit none
integer::d
double precision::x(d),omega(d),D_morse,V
V=D_morse*sum((exp(-omega(:)*x(:))-1.)**2)
end function V
!==============================================================================!
function P(d,x,E_cut,omega,D_morse)
!==============================================================================!
!Target Distribution Function
!==============================================================================!
!P              ==>evaluate P(x)
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!==============================================================================!
implicit none
integer::d
double precision::x(d),E_cut,P,omega(d),D_morse
if(V(d,x,omega,D_morse)<E_cut) P=(E_cut-V(d,x,omega,D_morse))**(d/2.)/integral_P
if(V(d,x,omega,D_morse)>=E_cut) P=1d-20           !set equal to 0 if beyond Ecut
end function P
!==============================================================================!
subroutine box_size_P(d,N_MMC_box,E_cut,omega,D_morse)
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
integer::d,N_MMC_box,i,j
double precision::E_Cut,dummy,r_trial(d),r(d),s(d),omega(d),D_morse
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
  if(P(d,r_trial,E_cut,omega,D_morse)/P(d,r,E_cut,omega,D_morse).ge.dummy) then
    r=r_trial
    do j=1,d
      if(xmin(j)>r(j)) xmin(j)=r(j)
      if(xmax(j)<r(j)) xmax(j)=r(j)
    enddo
  endif
enddo
end subroutine box_size_P
!==============================================================================!
subroutine compute_integral_P(d,N_1D,E_cut,omega,D_morse)
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
integer::d,N_1D,i,j,index1(d),Ntotal
double precision::r(d),Moment,dummy,delx(d),E_cut,omega(d),D_morse
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
  dummy=P(d,r,E_cut,omega,D_morse)
  Moment=Moment+dummy
enddo
dummy=1./N_1D**d
do j=1,d
  dummy=dummy*(xmax(j)-xmin(j))
enddo
integral_P=dummy*Moment
end subroutine compute_integral_P
!==============================================================================!
subroutine write_grid(d,x,id,grid_name)
!==============================================================================!
!write grid to a file
!==============================================================================!
implicit none
integer::d,id,j
double precision::x(d,Npoints)
character(len=12)::grid_name
open(unit=id,file=grid_name)
do j=1,Npoints
  write(id,*) x(:,j)
enddo
close(id)
end subroutine write_grid
!==============================================================================!
subroutine initial_distribution(d,x,E_cut,omega,D_morse)
!==============================================================================!
!Generate initial distribution of points, accept anything within E_cut contour
!==============================================================================!
implicit none
integer::d,i
double precision::E_Cut,x(d,Npoints),s(d),omega(d),D_morse
i=1
do while(i.le.Npoints)
  call random_number(s)
  s(:)=xmin(:)+s(:)*(xmax(:)-xmin(:))
  if(V(d,s,omega,D_morse)<E_cut)then
    x(:,i)=s(:)
    i=i+1
  endif
enddo
call write_grid(d,x,17,'coor_ini.dat')
end subroutine initial_distribution
!==============================================================================!
function Pair_LJ_NRG(d,x1,x2,E_cut,c_LJ,omega,D_morse)
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
integer::d
double precision::x1(d),x2(d),a,b,Pair_LJ_NRG,sigma1,sigma2,E_cut,c_LJ,omega(d)
double precision::D_morse
a=sum((x1(:)-x2(:))**2)
sigma1=c_LJ*(P(d,x1,E_cut,omega,D_morse)*Npoints)**(-1./d)
sigma2=c_LJ*(P(d,x2,E_cut,omega,D_morse)*Npoints)**(-1./d)
b=(sigma2**2/a)
a=(sigma1**2/a)
Pair_LJ_NRG=a**(d+9)-a**(d+3)+b**(d+9)-b**(d+3)
end function Pair_LJ_NRG
!==============================================================================!
subroutine initial_pairwise_energy(d,x,U,E_cut,c_LJ,omega,D_morse)
!==============================================================================!
!Compute the pairwise energies for all the initial Grid Points U[x_ij]
!==============================================================================!
implicit none
integer::d,i,j
double precision::E_cut,c_LJ,U(Npoints,Npoints),x(d,Npoints),omega(d),D_morse
do i=2,Npoints
  do j=1,i-1
    U(i,j)=Pair_LJ_NRG(d,x(:,i),x(:,j),E_cut,c_LJ,omega,D_morse)
    U(j,i)=U(i,j)
  enddo
enddo
end subroutine initial_pairwise_energy
!==============================================================================!
function random_integer(Nmin,Nmax)
!==============================================================================!
!Randomly generate an integer in the range Nmin-Nmax
!==============================================================================!
!Nmin           ==>minimum index value
!Nmax           ==>maximum index value
!random_integer ==>integer returned
!a              ==>pseudo-random number (0,1)
!==============================================================================!
implicit none
integer::Nmin,Nmax,random_integer
double precision::a
call random_number(a)
random_integer=floor(a*(Nmax-Nmin+1))+Nmin
end function random_integer
!==============================================================================!
subroutine generate_grid(d,x,U,E_cut,c_LJ,omega,D_morse,N_MMC_grid,MMC_freq)
!==============================================================================!
!                           Generate QRG using MMC
!==============================================================================!
implicit none
integer::d,N_MMC_grid,MMC_freq,accept,counter,i,j,k
double precision::x(d,Npoints),U(Npoints,Npoints),E_cut,c_LJ,mv_cutoff,deltae1
double precision::Delta_E,U_move(Npoints),s(d),x0(d),omega(d),D_morse
!==============================================================================!
accept=0
counter=0
mv_cutoff=0.01
deltae1=0
do i=1,N_MMC_grid
  k=random_integer(1,Npoints)                               !Select Atom to Move
!==============================================================================!
!                           Generate trial move
!        random numbers generated (0,1), make it (-1,1) ==> s=2*s-1
!==============================================================================!
  call random_number(s)
  x0=x(:,k)+mv_cutoff*(2*s-1)
!==============================================================================!
!                   Only consider point if V(trial)<Ecut
!                 Compute Energy Change due to Trial Move
!==============================================================================!
  if(V(d,x0,omega,D_morse).lt.E_cut) then
    counter=counter+1
    U_move(k)=P(d,x0,E_cut,omega,D_morse)
    Delta_E=0d0
    do j=1,Npoints
      if(j.ne.k) then
        U_move(j)=Pair_LJ_NRG(d,x(:,j),x0,E_cut,c_LJ,omega,D_morse)
        Delta_E=Delta_E+U(j,k)-U_move(j)
      endif
    enddo
!==============================================================================!
!               Accept any trial that decreases the energy
!==============================================================================!
    if(Delta_E.ge.0d0)then
      U(:,k)=U_move(:)
      U(k,:)=U_move(:)
      accept=accept+1
      x(:,k)=x0(:)
      deltae1=deltae1+Delta_E
    endif
  endif
!==============================================================================!
! for MMC want acceptance ~50%, adjust trial movement displacement accordingly
!==============================================================================!
  if(mod(i,MMC_freq)==0)then
    if(dble(accept)/counter<0.3)then
      mv_cutoff=mv_cutoff*0.9
    else
      mv_cutoff=mv_cutoff*1.1
    endif
  accept=0
  counter=0
  deltae1=0
  endif
enddo
call write_grid(d,x,18,'grid_fin.dat')
!==============================================================================!
end subroutine generate_grid
!==============================================================================!
subroutine write_out(d,E_cut,N_1D,c_LJ,omega,D_morse,N_MMC_grid,MMC_freq)
!==============================================================================!
!write output file
!==============================================================================!
implicit none
integer::d,N_1D,N_MMC_grid, MMC_freq,i
double precision::E_cut,c_LJ,D_morse,omega(d)
open(unit=99,file='out')
write(99,*) 'd ==> ', d
write(99,*) 'D_morse ==> ',D_morse
do i=1,d
  write(99,*) 'Box Dimensions==>', xmin(i),xmax(i)
enddo
write(99,*) 'omega ==> ',omega
write(99,*) 'integral_P ==> ',integral_P
write(99,*) 'E_cut ==> ',E_cut
write(99,*) 'N_1D ==> ',N_1D
write(99,*) 'Npoints ==> ',Npoints
write(99,*) 'c_LJ ==> ',c_LJ
write(99,*) 'N_MMC_grid ==> ',N_MMC_grid
write(99,*) 'MMC_freq ==> ',MMC_freq
close(99)
end subroutine write_out
!==============================================================================!
end module make_grid_mod
