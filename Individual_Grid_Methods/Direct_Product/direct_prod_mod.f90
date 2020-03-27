!=================20==================40==================60==================80
!                    Direct-Product Grid Generation Module                      !
!==============================================================================!
!Generate a d-dimensional (uniform) square Direct-Product Grid within E_cut
!Metropolis Monte Carlo is used to determine the box size
!If you added a new potential to the potentials_mod.f90 be sure to update the
!potentials subroutine.
!==============================================================================!
!       Modified:
!   26 March 2020
!       Author:
!   Shane Flynn
!==============================================================================!
module dp_grid_mod
implicit none
!==============================================================================!
contains
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
character(len=20)::potential
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
subroutine write_out(potential,Npoints,d,E_cut,xmin,xmax,N_1D,N_MMC_box,&
  integral_P)
!==============================================================================!
!Write output file
!==============================================================================!
!potential          ==>Potential name
!Npoints            ==>Number of points to generate
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!E_cut              ==>Distribution cutoff contour
!xmin(d)            ==>Minimum of normalization box size
!xmax(d)            ==>Maximum of normalization box size
!N_1D               ==>Number of points in 1 dimension for computing integral_P
!N_MMC_box          ==>Number of MMC Iterations to determine box-size
!integral_P         ==>Normalization constant for the distribtion P(x)
!==============================================================================!
implicit none
integer::Npoints,d,N_1D,i,N_MMC_box
double precision::E_cut,integral_P,xmin(d),xmax(d)
character(len=20)::potential
!==============================================================================!
open(unit=99,file='out')
write(99,*) 'd ==> ', d
write(99,*) 'potential ==> ', potential
do i=1,d
  write(99,*) 'Box Dimensions==>', xmin(i),xmax(i)
enddo
write(99,*) 'integral_P ==> ',integral_P
write(99,*) 'E_cut ==> ',E_cut
write(99,*) 'N_1D ==> ',N_1D
write(99,*) 'Npoints ==> ',Npoints
write(99,*) 'N_MMC_box ==> ',N_MMC_box
close(99)
end subroutine write_out
!==============================================================================!
end module dp_grid_mod
