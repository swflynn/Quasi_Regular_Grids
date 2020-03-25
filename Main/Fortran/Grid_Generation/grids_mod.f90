!=================20==================40==================60==================80
!                          Grid Generation Module
!==============================================================================!
!Generate a d-dimensional grid 
!All methods use Metropolis Monte Carlo to determine the box size for evaluating 
!the distribution.
!A (uniform) square grid is then used to normalize the distribution
!The following methods are available:
!==============================================================================!
!Quasi-Regular:
!QRG is optimized using a quasi-Lennard Jones Potential
!Forced minimization; accept any trial moves that reduces the system's energy
!==============================================================================!
!Direct Product: 
!A square grid is generated such that each point is within E_{cut}
!==============================================================================!
!Pseudo-Random:
!A pseudo-random number is generated and then evaluated using either the cutoff 
!contour or the rejection method to determine if the point is kept. 
!==============================================================================!
!Quasi-Random:
!A Quasi-random number is generated and then evaluated using either the cutoff 
!contour or the rejection method to determine if the point is kept. 
!In this work the Quasi-Random Sobol Sequence is used.
!==============================================================================!
!If you added a new potential to the potentials_mod.f90 be sure to update  the
!potentials subroutine.
!==============================================================================!
!       Modified:
!   24 March 2020
!       Author:
!   Shane Flynn
!==============================================================================!
module grid_mod
implicit none
!==============================================================================!
contains
!==============================================================================!
function random_integer(Nmin,Nmax)
!==============================================================================!
!Randomly generate an integer in the range Nmin-Nmax
!==============================================================================!
!Nmin           ==>minimum index value
!Nmax           ==>maximum index value
!a              ==>uniform pseudo-random number
!random_integer ==>integer returned
!==============================================================================!
implicit none
integer::Nmin,Nmax,random_integer
double precision::a
!==============================================================================!
call random_number(a)
random_integer=floor(a*(Nmax-Nmin+1))+Nmin
end function random_integer
!==============================================================================!
subroutine potentials(potential,d,x_i,V_i)
!==============================================================================!
!Determines which Potential Energy Function to call
!See the potentials_mod.f90 for available potentials, and update this 
!subroutine accordingly.
!All equations derived for a d-dimensional system, any potential is valid.
!==============================================================================!
!potential        ==>potential name
!d                ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)           ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i              ==>Potential Energy evaluation V(x_i)
!==============================================================================!
use potential_mod
implicit none
integer::d
double precision::x_i(d),V_i
character(len=10)::potential
!==============================================================================!
if(potential=='morse_pot') then
  call morse(d,x_i,v_i)
elseif(potential=='henon_pot') then
  call henon(d,x_i,v_i)
else
  stop 'Cannot Identify Potential, Check "potentials" Subroutine'
endif
end subroutine potentials
!==============================================================================!
subroutine grid_type(grid,potential,reject,Npoints,d,x,U,V_i,E_cut,xmin,xmax,&
N_1D,N_MMC_box,c_LJ,N_MMC_grid,MMC_freq,integral_P)
!==============================================================================!
!Determines which grid generation method to call:
!Quasi-Regular, Quasi-Random, Pseudo-Random, Direct-Product
!==============================================================================!
!grid               ==>Grid type
!potential          ==>Potential name
!reject             ==>Rejection logical (pseudo-random/quasi-random grids only)
!Npoints            ==>Number of points to generate
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x(d,Npoints)       ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!U(Npoints,Npoints) ==>Pairwise energies for all i-j points
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!xmin(d)            ==>Minimum of normalization box size
!xmax(d)            ==>Maximum of normalization box size
!N_1D               ==>Number of points in 1 dimension for computing integral_P
!N_MMC_box          ==>Number of MMC Iterations to determine box-size
!c_LJ               ==>Parameter for q-LJ pseudo-potential
!N_MMC_grid         ==>Number of MMC Iterations to optimize QRG
!MMC_freq           ==>Frequency to update QRG MMC grid mv_cutoff
!integral_P         ==>Normalization constant for the distribtion P(x)
!==============================================================================!
implicit none
integer::Npoints,d,N_MMC_box,N_1D,N_MMC_grid,MMC_freq
double precision::x(d,Npoints),U(Npoints,Npoints),E_cut,c_LJ
double precision::integral_P,V_i,xmin(d)
double precision::xmax(d)
character(len=20)::grid,potential
logical::reject
!==============================================================================!
if(grid=='direct-product') then
  call Direct(potential,d,V_i,E_cut,xmin,xmax,N_1D,N_MMC_box,integral_P)
elseif(grid=='pseudo-random')then
  call P_Rand(potential,reject,Npoints,d,x,V_i,E_cut,xmin,xmax,N_1D,&
  N_MMC_box,integral_P)
  return
elseif(grid=='quasi-random')then
  call Q_Rand(potential,reject,Npoints,d,x,V_i,E_cut,xmin,xmax,N_1D,&
  N_MMC_box,integral_P)
elseif(grid=='quasi-regular') then
  call Q_Reg(potential,Npoints,d,x,U,V_i,E_cut,xmin,xmax,N_1D,N_MMC_box,c_LJ,&
  N_MMC_grid,MMC_freq,integral_P)
else
  stop 'Cannot identify the grid, Check "grid_type" Subroutine in grids_mod.f90'
endif
end subroutine grid_type
!==============================================================================!
function P_i(potential,d,x_i,V_i,E_cut,integral_P)
!==============================================================================!
!Target Distribution Function, !defined according to a semi-classical argument:
!B. Poirier, “Algebraically self-consistent quasiclassical approximation on 
!phase space,” Found. Phys. 30, 1191–1226 (2000).
!This subroutine computes P(x_i)
!==============================================================================!
!potential          ==>Potential name
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)             ==>i-th grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!integral_P         ==>Normalization constant for the distribtion P(x)
!P_i                ==>evaluate P(x_i)
!==============================================================================!
implicit none
integer::d
double precision::x_i(d),E_cut,P_i,integral_P,V_i
character(len=10)::potential
!==============================================================================!
call potentials(potential,d,x_i,V_i)
if(V_i<E_cut) P_i=(E_cut-V_i)**(d/2.)/integral_P
if(V_i>=E_cut) P_i=1d-20                      !Define distribution=0 beyond Ecut
end function P_i
!==============================================================================!
subroutine box_size_P(potential,d,V_i,E_cut,xmin,xmax,N_MMC_box,integral_P)
!==============================================================================!
!Determine the box size for normalizing P(x) using Metropolis Monte Carlo
!This subroutine Computes xmin(d) and xmax(d)
!==============================================================================!
!potential          ==>Potential name
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!xmin(d)            ==>Minimum of normalization box size
!xmax(d)            ==>Maximum of normalization box size
!N_MMC_box          ==>Number of MMC Iterations to determine box-size
!integral_P         ==>Normalization constant for the distribtion P(x)
!==============================================================================!
integer::d,N_MMC_box,i,j
double precision::E_Cut,dummy,r_trial(d),r(d),s(d),integral_P,xmin(d),xmax(d)
double precision::V_i
double precision,parameter::mv_cutoff=0.1
character(len=10)::potential
!==============================================================================!
integral_P=1d0                              !Initially set to 1 so you cancall P
r=0d0
xmin=r
xmax=r
do i=1,N_MMC_box
  call random_number(s)
  r_trial=r+mv_cutoff*(2*s-1) 
!trial move needs to be (-1,1), random numbers are (0,1) ===>s=2*s-1 Transform
  call random_number(dummy)                           !MMC acceptance criteria
  if(P_i(potential,d,r_trial,V_i,E_cut,integral_P)/&
  P_i(potential,d,r,V_i,E_cut,integral_P).ge.dummy) then
    r=r_trial
    do j=1,d
      if(xmin(j)>r(j)) xmin(j)=r(j)
      if(xmax(j)<r(j)) xmax(j)=r(j)
    enddo
  endif
enddo
end subroutine box_size_P
!==============================================================================!
subroutine compute_integral_P(potential,d,V_i,E_cut,xmin,xmax,N_1D,integral_P)
!==============================================================================!
!Use a (uniform) square grid to integrate P(x)
!Box size for the grid is determined in the box_size subroutine
!int P(x)~Area_Square/N sum_n=1,N P(x_n)
!This subroutine Computes integral_P
!==============================================================================!
!potential          ==>Potential name
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!xmin(d)            ==>Minimum of normalization box size
!xmax(d)            ==>Maximum of normalization box size
!N_1D               ==>Number of points in 1 dimension for computing integral_P
!integral_P         ==>Normalization constant for the distribtion P(x)
!Ntotal             ==>Total number of evaluations for all dimensions
!Moment             ==>First Moment for the distribution 
!==============================================================================!
integer::d,N_1D,i,j,index1(d),Ntotal
double precision::r(d),Moment,dummy,delx(d),E_cut,integral_P,xmin(d),xmax(d),V_i
character(len=10)::potential
!==============================================================================!
open(unit=55,file='direct_grid.dat')
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
  dummy=P_i(potential,d,r,V_i,E_cut,integral_P)
  Moment=Moment+dummy
  if(V_i<E_cut) write(55,*) r
enddo
dummy=1./N_1D**d
do j=1,d
  dummy=dummy*(xmax(j)-xmin(j))
enddo
integral_P=dummy*Moment
close(55)
end subroutine compute_integral_P
!==============================================================================!
subroutine Direct(potential,d,V_i,E_cut,xmin,xmax,N_1D,N_MMC_box,integral_P)
!==============================================================================!
!Generate a Direct Product Grid
!==============================================================================!
!potential          ==>Potential name
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!xmin(d)            ==>Minimum of normalization box size
!xmax(d)            ==>Maximum of normalization box size
!N_1D               ==>Number of points in 1 dimension for computing integral_P
!N_MMC_box          ==>Number of MMC Iterations to determine box-size
!integral_P         ==>Normalization constant for the distribtion P(x)
!==============================================================================!
implicit none
integer::d,N_MMC_box,N_1D
double precision::E_cut,integral_P,V_i,xmin(d),xmax(d)
character(len=20)::potential
!==============================================================================!
call box_size_P(potential,d,V_i,E_cut,xmin,xmax,N_MMC_box,integral_P)
call compute_integral_P(potential,d,V_i,E_cut,xmin,xmax,N_1D,integral_P)
end subroutine Direct
!==============================================================================!
subroutine P_Rand(potential,reject,Npoints,d,x,V_i,E_cut,xmin,xmax,N_1D,&
  N_MMC_box,integral_P)
!==============================================================================!
!Generate grid using a Pseudo-Random sequence of points, grid can be generated
!using the cutoff contour as a boundry line, or by using the rejection method.
!==============================================================================!
!potential          ==>Potential name
!reject             ==>Rejection logical (pseudo-random/quasi-random grids only)
!Npoints            ==>Number of points to generate
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x(d,Npoints)       ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!xmin(d)            ==>Minimum of normalization box size
!xmax(d)            ==>Maximum of normalization box size
!N_1D               ==>Number of points in 1 dimension for computing integral_P
!N_MMC_box          ==>Number of MMC Iterations to determine box-size
!integral_P         ==>Normalization constant for the distribtion P(x)
!unif_tol           ==>tolerance for accepting based on Ecut alone
!==============================================================================!
implicit none
integer::Npoints,d,N_MMC_box,N_1D,count_acc,i
double precision::x(d,Npoints),E_cut,integral_P,z(d),V_i,xmin(d),xmax(d),accept
double precision,parameter::unif_tol=0.00000001
character(len=20)::potential
logical::reject
!==============================================================================!
call box_size_P(potential,d,V_i,E_cut,xmin,xmax,N_MMC_box,integral_P)
call compute_integral_P(potential,d,V_i,E_cut,xmin,xmax,N_1D,integral_P)
count_acc=0
do while(count_acc.lt.Npoints)
  call pseudo_unif(d,z,xmin,xmax)
  call random_number(accept)
!==============================================================================!
  if(reject.eqv..true.)then                       !generate grid using rejection
    if(P_i(potential,d,z,V_i,E_cut,integral_P)/(E_cut).gt.accept)then
      count_acc=count_acc+1
      x(:,count_acc)=z(:)
    endif
  else                                             !generate grid based on E_cut
    if(P_i(potential,d,z,V_i,E_cut,integral_P).gt.unif_tol)then
      count_acc=count_acc+1
      x(:,count_acc)=z(:)
    endif
  endif
enddo
open(unit=18,file='pseudo_random.dat')
do i=1,Npoints
  write(18,*) x(:,i)
enddo
close(18)
end subroutine P_Rand
!==============================================================================!
subroutine pseudo_unif(d,x_unif,xmin,xmax)
!==============================================================================!
!Generate a Uniform Distribution of Pseudo-Random Points
!Scale the points from [0,1] ==> [xmin,xmax]
!==============================================================================!
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_unif(d)          ==>Uniform Pseudo-Random point 
!xmin(d)            ==>Minimum of normalization box size
!xmax(d)            ==>Maximum of normalization box size
!==============================================================================!
implicit none
integer::d
double precision::x_unif(d),xmin(d),xmax(d)
!==============================================================================!
call random_number(x_unif)
x_unif=xmin+x_unif*(xmax-xmin)        !Scale uniform distribution to span domain
end subroutine pseudo_unif
!==============================================================================!
subroutine Q_Rand(potential,reject,Npoints,d,x,V_i,E_cut,xmin,xmax,N_1D,&
  N_MMC_box,integral_P)
!==============================================================================!
!Generate grid using a Quasi-Random sequence of points, grid can be generated
!using the cutoff contour as a boundry line, or by using the rejection method.
!Uses the Sobol Sequence from sobol.f90 made available by John Burkardt
!https://people.sc.fsu.edu/~jburkardt/f_src/sobol/sobol.html under GNU LGPL
!==============================================================================!
!potential          ==>Potential name
!reject             ==>Rejection logical (pseudo-random/quasi-random grids only)
!Npoints            ==>Number of points to generate
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x(d,Npoints)       ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!xmin(d)            ==>Minimum of normalization box size
!xmax(d)            ==>Maximum of normalization box size
!N_1D               ==>Number of points in 1 dimension for computing integral_P
!N_MMC_box          ==>Number of MMC Iterations to determine box-size
!integral_P         ==>Normalization constant for the distribtion P(x)
!unif_tol           ==>tolerance for accepting based on Ecut alone
!==============================================================================!
implicit none
integer::Npoints,d,N_MMC_box,N_1D,count_acc,i
integer*8::skip                                      !*8 needed for sobol module
double precision::x(d,Npoints),E_cut,integral_P,z(d),V_i,xmin(d),xmax(d),accept
double precision,parameter::unif_tol=0.00000001
character(len=20)::potential
logical::reject
!==============================================================================!
call box_size_P(potential,d,V_i,E_cut,xmin,xmax,N_MMC_box,integral_P)
call compute_integral_P(potential,d,V_i,E_cut,xmin,xmax,N_1D,integral_P)
skip=Npoints                              !skip is for the sobol generator seed
count_acc=0
do while(count_acc.lt.Npoints)              
  z=0d0
  call sobol_unif(d,skip,z,xmin,xmax)
  call random_number(accept)
  if(reject.eqv..true.)then                       !generate grid using rejection
    if(P_i(potential,d,z,V_i,E_cut,integral_P)/(E_cut).gt.accept)then
      count_acc=count_acc+1
      x(:,count_acc)=z(:)
    endif
  else                                             !generate grid based on E_cut
    if(P_i(potential,d,z,V_i,E_cut,integral_P)/(E_cut).gt.unif_tol)then
      count_acc=count_acc+1
      x(:,count_acc)=z(:)
    endif
  endif
enddo
open(unit=18,file='quasi_random.dat')
do i=1,Npoints
  write(18,*) x(:,i)
enddo
close(18)
end subroutine Q_Rand
!==============================================================================!
subroutine sobol_unif(d,skip,x_unif,xmin,xmax)
!==============================================================================!
!Generates a uniform Quasi-Random point
!Uses the Sobol Sequence from sobol.f90 made available by John Burkardt
!https://people.sc.fsu.edu/~jburkardt/f_src/sobol/sobol.html under GNU LGPL
!==============================================================================!
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!skip               ==>Seed for the random number generator
!x_unif(d)          ==>Uniform Quasi-Random point 
!xmin(d)            ==>Minimum of normalization box size
!xmax(d)            ==>Maximum of normalization box size
!==============================================================================!
use sobol
implicit none
integer::d
INTEGER(kind = 8)::skip
double precision::x_unif(d),xmin(d),xmax(d)
!==============================================================================!
x_unif=i8_sobol(int(d, 8), skip)
x_unif=xmin+x_unif*(xmax-xmin)        !Scale uniform distribution to span domain
end subroutine sobol_unif
!==============================================================================!
subroutine Q_Reg(potential,Npoints,d,x,U,V_i,E_cut,xmin,xmax,N_1D,N_MMC_box,&
c_LJ,N_MMC_grid,MMC_freq,integral_P)
!==============================================================================!
!Generate grid using a Quasi-Regular sequence of points
!Grid is generated using MMC and forced minimization using our pseudo-potential
!==============================================================================!
!potential          ==>Potential name
!Npoints            ==>Number of points to generate
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x(d,Npoints)       ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!U(Npoints,Npoints) ==>Pairwise energies for all i-j points
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!xmin(d)            ==>Minimum of normalization box size
!xmax(d)            ==>Maximum of normalization box size
!N_1D               ==>Number of points in 1 dimension for computing integral_P
!N_MMC_box          ==>Number of MMC Iterations to determine box-size
!c_LJ               ==>Parameter for q-LJ pseudo-potential
!N_MMC_grid         ==>Number of MMC Iterations to optimize QRG
!MMC_freq           ==>Frequency to update QRG MMC grid mv_cutoff
!integral_P         ==>Normalization constant for the distribtion P(x)
!==============================================================================!
implicit none
integer::Npoints,d,N_MMC_box,N_1D,N_MMC_grid,MMC_freq,accept,counter,i,j,k
double precision::x(d,Npoints),U(Npoints,Npoints),E_cut,c_LJ,mv_cutoff,deltae1
double precision::Delta_E,U_move(Npoints),s(d),x0(d),integral_P,V_i,xmin(d)
double precision::xmax(d)
character(len=10)::potential
!==============================================================================!
call box_size_P(potential,d,V_i,E_cut,xmin,xmax,N_MMC_box,integral_P)
call compute_integral_P(potential,d,V_i,E_cut,xmin,xmax,N_1D,integral_P)
call initial_distribution(potential,Npoints,d,x,V_i,E_cut,xmin,xmax)
call initial_pair_energy(potential,Npoints,d,x,U,V_i,E_cut,c_LJ,integral_P)
accept=0
counter=0
mv_cutoff=0.01
deltae1=0
do i=1,N_MMC_grid
  k=random_integer(1,Npoints)                               !Select Atom to Move
  call random_number(s)
  x0=x(:,k)+mv_cutoff*(2*s-1)   !random number (0,1), make it (-1,1) ==> s=2*s-1
!==============================================================================!
  call potentials(potential,d,x0,V_i)
  if(V_i.lt.E_cut) then                          !Only consider if V(trial)<Ecut
    counter=counter+1
    U_move(k)=P_i(potential,d,x0,V_i,E_cut,integral_P) 
    Delta_E=0d0
    do j=1,Npoints
      if(j.ne.k) then
        U_move(j)=Pair_LJ_NRG(potential,Npoints,d,x(:,j),x0,V_i,E_cut,c_LJ,&
            integral_P)
        Delta_E=Delta_E+U(j,k)-U_move(j)        !Energy change due to trial move
      endif
    enddo
    if(Delta_E.ge.0d0)then      !Forced minimization: accept if energy decreases
      U(:,k)=U_move(:)
      U(k,:)=U_move(:)
      accept=accept+1
      x(:,k)=x0(:)
      deltae1=deltae1+Delta_E
    endif
  endif
!for MMC want acceptance ~30-50%, adjust trial movement displacement accordingly
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
open(unit=18,file='qrg.dat')
do i=1,Npoints
  write(18,*) x(:,i)
enddo
close(18)
!==============================================================================!
end subroutine Q_Reg
!==============================================================================!
subroutine initial_distribution(potential,Npoints,d,x,V_i,E_cut,xmin,xmax)
!==============================================================================!
!Generate an initial distribution to optimize with the Quasi-Regular algorithm
!Accepts any pseudo-random number generated within the cutoff contour E_cut
!==============================================================================!
!potential          ==>Potential name
!Npoints            ==>Number of points to generate
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x(d,Npoints)       ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!xmin(d)            ==>Minimum of normalization box size
!xmax(d)            ==>Maximum of normalization box size
!==============================================================================!
implicit none
integer::Npoints,d,i
double precision::E_Cut,x(d,Npoints),s(d),xmin(d),xmax(d),V_i
character(len=10)::potential
!==============================================================================!
i=1
do while(i.le.Npoints)
  call random_number(s)
  s(:)=xmin(:)+s(:)*(xmax(:)-xmin(:))
  call potentials(potential,d,s,V_i)
  if(V_i<E_cut)then
    x(:,i)=s(:)
    i=i+1
  endif
enddo
open(unit=17,file='grid_ini.dat')
do i=1,Npoints
  write(17,*) x(:,i)
enddo
close(17)
end subroutine initial_distribution
!==============================================================================!
function Pair_LJ_NRG(potential,Npoints,d,x1,x2,V_i,E_cut,c_LJ,integral_P)
!==============================================================================!
!Computes the quasi-Lennard Jones pairwise energy between grid points used in
!our QRG algorithm
!This function computes the q-LJ energy between 2 grid-points
!==============================================================================!
!potential          ==>Potential name
!Npoints            ==>Number of points to generate
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x_i(d)             ==>Grid points coordinates (x^i=x_1,x_2,...x_d)
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!c_LJ               ==>Parameter for q-LJ pseudo-potential
!integral_P         ==>Normalization constant for the distribtion P(x)
!Pair_LJ_NRG        ==>Energy of the i-j q-LJ potential
!==============================================================================!
implicit none
integer::Npoints,d
double precision::x1(d),x2(d),a,b,Pair_LJ_NRG,sigma1,sigma2,E_cut,c_LJ,V_i
double precision::integral_P
character(len=10)::potential
!==============================================================================!
a=sum((x1(:)-x2(:))**2)
sigma1=c_LJ*(P_i(potential,d,x1,V_i,E_cut,integral_P)*Npoints)**(-1./d)
sigma2=c_LJ*(P_i(potential,d,x2,V_i,E_cut,integral_P)*Npoints)**(-1./d)
b=(sigma2**2/a)
a=(sigma1**2/a)
Pair_LJ_NRG=a**(d+9)-a**(d+3)+b**(d+9)-b**(d+3)
end function Pair_LJ_NRG
!==============================================================================!
subroutine initial_pair_energy(potential,Npoints,d,x,U,V_i,E_cut,c_LJ,integral_P)
!==============================================================================!
!Compute the pairwise energies for all the initial Grid Points U[x_{ij}]
!==============================================================================!
implicit none
integer::Npoints,d,i,j
double precision::E_cut,c_LJ,U(Npoints,Npoints),x(d,Npoints),integral_P,V_i
character(len=10)::potential
!==============================================================================!
!potential          ==>Potential name
!Npoints            ==>Number of points to generate
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!x(d,Npoints)       ==>All grid points coordinates (x^i=x_1,x_2,...x_d)
!U(Npoints,Npoints) ==>Pairwise energies for all i-j points
!V_i                ==>Potential Energy evaluation V(x_i)
!E_cut              ==>Distribution cutoff contour
!c_LJ               ==>Parameter for q-LJ pseudo-potential
!integral_P         ==>Normalization constant for the distribtion P(x)
!==============================================================================!
do i=2,Npoints
  do j=1,i-1
    U(i,j)=Pair_LJ_NRG(potential,Npoints,d,x(:,i),x(:,j),V_i,E_cut,c_LJ,&
        integral_P)
    U(j,i)=U(i,j)
  enddo
enddo
end subroutine initial_pair_energy
!==============================================================================!
subroutine write_out(grid,potential,reject,Npoints,d,E_cut,xmin,xmax,N_1D,&
        N_MMC_box,c_LJ,N_MMC_grid,MMC_freq,integral_P)
!==============================================================================!
!Write output file
!==============================================================================!
!grid               ==>Grid type
!potential          ==>Potential name
!reject             ==>Rejection logical (pseudo-random/quasi-random grids only)
!Npoints            ==>Number of points to generate
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!E_cut              ==>Distribution cutoff contour
!xmin(d)            ==>Minimum of normalization box size
!xmax(d)            ==>Maximum of normalization box size
!N_1D               ==>Number of points in 1 dimension for computing integral_P
!N_MMC_box          ==>Number of MMC Iterations to determine box-size
!c_LJ               ==>Parameter for q-LJ pseudo-potential
!N_MMC_grid         ==>Number of MMC Iterations to optimize QRG
!MMC_freq           ==>Frequency to update QRG MMC grid mv_cutoff
!integral_P         ==>Normalization constant for the distribtion P(x)
!==============================================================================!
implicit none
integer::Npoints,d,N_1D,N_MMC_grid, MMC_freq,i,N_MMC_box
double precision::E_cut,c_LJ,integral_P,xmin(d),xmax(d)
character(len=20)::grid,potential
logical::reject
!==============================================================================!
open(unit=99,file='out')
write(99,*) 'd ==> ', d
write(99,*) 'grid type ==> ', grid
write(99,*) 'potential ==> ', potential
write(99,*) 'use rejection ? ==> ', reject
do i=1,d
  write(99,*) 'Box Dimensions==>', xmin(i),xmax(i)
enddo
write(99,*) 'integral_P ==> ',integral_P
write(99,*) 'E_cut ==> ',E_cut
write(99,*) 'N_1D ==> ',N_1D
write(99,*) 'Npoints ==> ',Npoints
write(99,*) 'c_LJ ==> ',c_LJ
write(99,*) 'N_MMC_box ==> ',N_MMC_box
write(99,*) 'N_MMC_grid ==> ',N_MMC_grid
write(99,*) 'MMC_freq ==> ',MMC_freq
close(99)
end subroutine write_out
!==============================================================================!
end module grid_mod
