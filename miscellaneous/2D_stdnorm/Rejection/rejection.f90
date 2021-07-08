!=============================================================================80
!                       Gaussian (Sobol + Rejection)
!=============================================================================80
!     Discussion:
!Fortran 90 Implementation of the Rejection method for a 2D-Gaussian pdf 
!pseudo/quasi random points are generated over the domain [xmin,xmax] and then 
!rejection is applied to generate a set of points according to the pdf.
!==============================================================================!
!     Modified:
!8 Just 2021
!     Author:
!Shane Flynn
!==============================================================================!
module rej_mod
implicit none
!==============================================================================!
!                             Global Variables
!==============================================================================!
!d                  ==>Coordinate dimensionality (x^i=x_1,x_2,...x_d)
!==============================================================================!
integer::d
double precision,parameter::pie=acos(-1d0)
!==============================================================================!
contains
!==============================================================================!
subroutine pdf_eval(x,P)
!==============================================================================!
!normalized 2D Gaussian 
!==============================================================================!
!x        ==>(d) gridpoint Coordinates (uniform distibution)
!P        ==>Gaussian Gridpoint
!==============================================================================!
implicit none
double precision::x(d),p
P=(2.*pie)**(-d/2.)*exp(-0.5*sum(x(:)**2))
end subroutine pdf_eval
!==============================================================================!
subroutine unif_point(flag,skip,x_unif,xmin,xmax)
!==============================================================================!
!Generate a Uniform Distribution of Pseudo-Random or Quasi-Random Points
!Scale the points from [0,1] ==> [xmin,xmax]
!Uses the Sobol Sequence from sobol.f90 made available by John Burkardt
!https://people.sc.fsu.edu/~jburkardt/f_src/sobol/sobol.html under GNU LGPL
!==============================================================================!
!flag               ==>True==Pseudo-Random; False==Quasi-Random
!skip               ==>Seed for the quasi-random number generator
!x_unif             ==>(d) point from uniform distribtuion
!xmin               ==>(d) Minimum of domain size
!xmax               ==>(d) Maximum of domain size
!==============================================================================!
use sobol
implicit none
integer(kind=8)::skip
double precision::x_unif(d),xmin(d),xmax(d)
logical::flag
!==============================================================================!
if(flag)then
    call random_number(x_unif)
    x_unif=xmin+x_unif*(xmax-xmin)    !Scale uniform distribution to span domain
else
    x_unif=i8_sobol(int(d,8),skip)
    x_unif=xmin+x_unif*(xmax-xmin)        
endif
end subroutine unif_point
!==============================================================================!
!==============================================================================!
!==============================================================================!
end module rej_mod
!==============================================================================!
!==============================================================================!
program reject
use rej_mod
implicit none
!==============================================================================!
!                                Discussion
!==============================================================================!
!Npoints       ==>Number of Grid Points
!count_acc     ==>Number of Accepted Points
!count_rej     ==>Number of Rejections
!pxy           ==>PDF Evaluation
!accept        ==>Randomly Generated Acceptance Condition
!z             ==>(d) Quasi-Random Points
!x             ==>(d,Nsobol) Accepted Coordinates
!==============================================================================!
integer::Npoints,count_acc,count_rej,i
integer*8::ii,skip                                        !sobol module needs *8
double precision::pxy,accept,ixmin,ixmax
double precision,allocatable::z(:),x(:,:),xmin(:),xmax(:)
logical::flag
!==============================================================================!
!                               Read Input File
!==============================================================================!
read(*,*) d
read(*,*) Npoints
read(*,*) ixmin
read(*,*) ixmax
read(*,*) flag                     !true == pseudo-random, false == quasi-random
skip=Npoints
!==============================================================================!
!                           Allocations/Initialize
!==============================================================================!
allocate(z(d),x(d,Npoints),xmin(d),xmax(d))
xmin=ixmin          !in general you may want different cutoff for dif. dimension
xmax=ixmax
count_acc=0
count_rej=0
z=0d0
x=0d0
!==============================================================================!
!                       Store All Coordinates Generated
!==============================================================================!
open(unit=17,file='all_data.dat')
do while(count_acc<Npoints)
    call unif_point(flag,skip,z,xmin,xmax)
    write(17,*) z(:)
    call random_number(accept)
    call pdf_eval(z,pxy)
!==============================================================================!
!                           Accept/Reject Criteria
!==============================================================================!
    if(pxy>accept)then
        count_acc=count_acc +1
        x(:,count_acc)=z(:)
    else
        count_rej=count_rej + 1
    endif
enddo
close(17)
!==============================================================================!
!                        Write Accepted Points to File
!==============================================================================!
open(unit=18,file='grid.dat')
do i=1,Npoints
    write(18,*) x(:,i)
enddo
close(18)
!==============================================================================!
!                           Simulation Results
!==============================================================================!
write(*,*) 'Npoints ==> ', Npoints
write(*,*) 'Dimensionality', d
write(*,*) 'Total Rejections', count_rej
end program reject
