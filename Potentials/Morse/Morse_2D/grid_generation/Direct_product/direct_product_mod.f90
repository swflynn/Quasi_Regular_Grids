!=================20==================40==================60==================80
!                  2D Morse Direct Product Grid Generation
!=============================================================================80
!Generate gridpoints distributed via the 2D Morse Oscillator
!Direct Product grid within E_cut
!==============================================================================!
!       Modified:
!   25 Feburary 2020
!       Author:
!   Shane Flynn
!==============================================================================!
module Morse2D_direct_product_mod
implicit none
!==============================================================================!
!                              Global Variables
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
!                               Begin Module
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
end subroutine box_size_P
!==============================================================================!
subroutine generate_grid(N_1D,E_cut)
!==============================================================================!
!Generate a uniform direct product grid (and integrate P(r))
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
index1=0
Npoints=0
Ntotal=(N_1D+1)**d
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
  if(V(r)<E_cut)then
     write(*,*) r
     Npoints=Npoints+1
  endif
enddo
dummy=1./N_1D**d
do j=1,d
  dummy=dummy*(xmax(j)-xmin(j))
enddo
integral_P=dummy*Moment
end subroutine generate_grid
!==============================================================================!
subroutine write_out(E_cut,N_1D,N_MMC_box)
!==============================================================================!
!write output file
!==============================================================================!
implicit none
integer::i,N_1D,N_MMC_box
double precision::E_cut
open(unit=99,file='out')
write(99,*) 'd ==> ', d
write(99,*) 'D_morse ==> ',D_morse
write(99,*) 'omega ==> ',omega
write(99,*) 'N_MMC_box ==> ',N_MMC_box
do i=1,d
  write(99,*) 'Box Dimensions==>', xmin(i),xmax(i)
enddo
write(99,*) 'integral_P ==> ',integral_P
write(99,*) 'E_cut ==> ',E_cut
write(99,*) 'N_1D ==> ',N_1D
write(99,*) 'Npoints ==> ',Npoints
close(99)
end subroutine write_out
!==============================================================================!
end module Morse2D_direct_product_mod
