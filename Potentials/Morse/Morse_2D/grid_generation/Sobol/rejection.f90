!=============================================================================80
!                  2D Morse Sobol+Rejection Grid Generation
!=============================================================================80
!Generate gridpoints distributed via the 2D Morse Oscillator
!This grid is generated using a sobol sequence and the rejection method
!Force minimization, accept any trial moves that reduces the system's energy
!See MMC grid for metropolis optimization (acceptance criteria)
!==============================================================================!
!       Modified:
!   4 June 2019
!       Author:
!   Shane Flynn
!==============================================================================!
module sob_rej_mod
implicit none
!==============================================================================!
!                              Global Variables
!==============================================================================!
!d              ==> Particle Dimensionality
!Npoints        ==> Number of grid points
!E_cut          ==> Energy Cutoff Contour
!integral_P     ==> Normalization for P_x
!==============================================================================!
integer::d,Npoints
double precision::E_cut,integral_P
!==============================================================================!
!                               Begin Module
!==============================================================================!
contains
!==============================================================================!
function V(x)
!==============================================================================!
!Potential Energy (Hard-Coded 2D Morse)
!==============================================================================!
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!V              ==>evaluate V(x)
!D_morse        ==>Parameter for Morse Potential
!omega          ==>(d) Parameter for Morse Potential
!==============================================================================!
implicit none
double precision::x(d),V
double precision,parameter::omega(2)=(/0.2041241,0.18371169/)
double precision,parameter::D_morse=12.
V=D_morse*sum((exp(-omega(:)*x(:))-1.)**2)
end function V
!==============================================================================!
function P(x)
!==============================================================================!
!Target Distribution Function
!==============================================================================!
!P              ==>evaluate P(x)
!x              ==>(d) ith particles coordinates x^i_1,..,x^i_d
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
subroutine sobol_unif(d,skip,x_unif,xmin,xmax)
!==============================================================================!
!Generate a Uniform Distribution of Sobol Points for the 2D Morse Oscillator
!Uses the sobol number generator (sobol.f90)
!==============================================================================!
!xmin,xmax      ==>(d) distribution domain
!x_unif         ==>(d) Grid Points
!skip           ==>sobol generator seed
!==============================================================================!
use sobol
implicit none
integer::i
double precision::s(d),xmin(d),xmax(d)
INTEGER(kind = 4),INTENT(IN)::d
INTEGER(kind = 8),INTENT(IN)::skip
DOUBLE PRECISION,DIMENSION(d),INTENT(OUT)::x_unif
x_unif=i8_sobol(int(d, 8), skip)
!==============================================================================!
!scale uniform distribution to span domain
!==============================================================================!
x_unif=xmin+x_unif*(xmax-xmin)
end subroutine sobol_unif
!==============================================================================!
subroutine Moments_Reg(Moment,N,xmin,xmax)
!==============================================================================!
!Compute lower moments of the distribution to verify global accuracy
!Compute moments using a regular square grid (most accurate method for 3D Morse)
!Integrate over the square [a,b],[a,b] size determined by Moments_MMC subroutine
!int P(r)~Area_Square/N sum_n=1,N P(r_n)
!==============================================================================!
!r              ==>(d) coordinates
!xmin           ==>(d) minimum of normalization box
!xmax           ==>(d) maximum of normalization box
!Moment         ==>(0:5) 5 Lowest Moments for the distribution
!==============================================================================!
integer::N,i1,i2,j
double precision::r(d),xmin(d),xmax(d),Moment(0:5),dummy
Moment=0d0
do i1=0,N
    do i2=0,N
        r(1)=xmin(1)+i1*(xmax(1)-xmin(1))/N
        r(2)=xmin(2)+i2*(xmax(2)-xmin(2))/N
        dummy=P(r)
        Moment(0)=Moment(0)+dummy
        do j=1,d
            Moment(j)=Moment(j)+dummy*r(j)
            Moment(d+j)=Moment(d+j)+dummy*r(j)**2
        enddo
        Moment(5)=Moment(5)+dummy*r(1)*r(2)
    enddo
enddo
dummy=1./N**d
do j=1,d
    dummy=dummy*(xmax(j)-xmin(j))
enddo
integral_P=dummy*Moment(0)
Moment(1:5)=Moment(1:5)/Moment(0)
end subroutine Moments_Reg
!==============================================================================!
subroutine Moments_MMC(Moment,N_MC,xmin,xmax)
!==============================================================================!
!Compute low-level moments of the distribution to verify global accuracy
!Compute moments using Metropolis Monte Carlo
!This subroutine determins the box size for normalizing P; (xmin,xmax)
!==============================================================================!
!N_MC           ==>Number of Monte Carlo Iterations
!mv_cutoff      ==>trial displacement move cutoff
!r              ==>(d) coordinates
!r_trial        ==>(d) trial coordinates
!s              ==>(d) trail displacement; random number for coordinates
!xmin           ==>(d) minimum of normalization box
!xmax           ==>(d) maximum of normalization box
!Moment         ==>(0:5) 5 Lowest Moments for the distribution
!==============================================================================!
integer::N_MC,i,j
double precision::Moment(0:5),dummy,r_trial(d),r(d),s(d),xmin(d),xmax(d)
double precision,parameter::mv_cutoff=0.1
Moment=0d0
r=0d0
xmin=r
xmax=r
do i=1,N_MC
!==============================================================================!
!                   Generate coordinates for Random move
!           random numbers generated (0,1), make it (-1,1) ==> s=2*s-1
!==============================================================================!
    call random_number(s)
    r_trial=r+mv_cutoff*(2*s-1)
!==============================================================================!
!                               Test Acceptance
!==============================================================================!
    call random_number(dummy)
    if(P(r_trial)/P(r).ge.dummy) then
        r=r_trial
        do j=1,d
            if(xmin(j)>r(j)) xmin(j)=r(j)
            if(xmax(j)<r(j)) xmax(j)=r(j)
        enddo
    endif
    do j=1,d
        Moment(j)=Moment(j)+r(j)
        Moment(d+j)=Moment(d+j)+r(j)**2
    enddo
    Moment(5)=Moment(5)+r(1)*r(2)
enddo
Moment=Moment/N_MC
end subroutine Moments_MMC
!==============================================================================!
subroutine Moments(Moment,x)
!==============================================================================!
!Compute low-level moments of the distribution to verify global accuracy
!Compute moments given a set of points
!==============================================================================!
integer::i,j
double precision::Moment(0:5),x(d,Npoints)
Moment=0d0
do i=1,Npoints
    do j=1,d
        Moment(j)=Moment(j)+x(j,i)
        Moment(d+j)=Moment(d+j)+x(j,i)**2
    enddo
    Moment(5)=Moment(5)+x(1,i)*x(2,i)
enddo
Moment(1:5)=Moment(1:5)/Npoints
end subroutine Moments
!==============================================================================!
end module sob_rej_mod
!==============================================================================!
program main
use sob_rej_mod
!==============================================================================!
!N_MC_moments   ==># of MC moves to compute moments
!N_reg          ==># grid points in each dimension for moments with regular grid
!accept         ==># of acceptances for adjusting trial move acceptance
!Moment         ==>(0:5) 5 Lowest Moments for the distribution
!xmin           ==>(d) minimum for P(x) normalization box
!xmax           ==>(d) maximum of P(x) normalization box
!x              ==>(d,Npoints) All coordinates
!==============================================================================!
implicit none
integer::i,count_acc,count_rej,N_MC_moments,N_reg
integer*8::ii,skip
double precision::accept,pxy,moment(0:5)
double precision,allocatable::z(:),x(:,:),xmin(:),xmax(:)
!==============================================================================!
!                           Read Input Data File
!==============================================================================!
read(*,*) Npoints
read(*,*) E_cut
read(*,*) n_MC_moments
read(*,*) n_reg
skip=Npoints
d=2
allocate(z(d),x(d,Npoints),xmin(d),xmax(d))
!==============================================================================!
!                              Run Monte Carlo
!==============================================================================!
integral_P=1d0                       !set equal to 1 so you can initially call P
call Moments_MMC(Moment,N_MC_moments,xmin,xmax)
write(*,*) 'Moments from Metropolis MC using N_MC=', N_MC_moments
write(*,*) Moment(1:5)
do i=1,d
    write(*,*) i,' xmin xmax ==>', xmin(i),xmax(i)
enddo
xmin=xmin-(xmax-xmin)*0.01
xmax=xmax+(xmax-xmin)*0.01
call Moments_Reg(Moment,N_reg,xmin,xmax)
write(*,*) 'Moments using regular integration:'
write(*,*) Moment(1:5)
write(*,*) 'Integral of P =', Integral_P
!==============================================================================!
count_acc=0
count_rej=0
z=0d0
x=0d0
open(unit=66,file='all_data.dat')
!==============================================================================!
!                     Generate Points Until total=Nsobol
!==============================================================================!
do while (count_acc .lt. Npoints)
    call sobol_unif(d,skip,z(:),xmin,xmax)
    call random_number(accept)
!==============================================================================!
!                    Use PDF to Accept/Reject Coordinates
!==============================================================================!
    if(P(z)/(E_cut).gt.accept) then                    !use for rejection method
!    if(P(z) .gt. 0.00000001) then            !use for uniform sobol with cutoff
        count_acc=count_acc +1
        x(:,count_acc)=z(:)
    else
        count_rej=count_rej + 1
    endif
end do
close(66)
open(unit=99,file='grid_sobol.dat')
do i=1,Npoints
    write(99,*) x(:,i)
enddo
close(99)
open(unit=67,file='sim.dat')
write(67,*) 'Npoints', Npoints
write(67,*) 'Dimensionality', d
write(67,*) 'Integral_P', integral_P
write(67,*) 'Total Rejections', count_rej
close(67)
end program main
