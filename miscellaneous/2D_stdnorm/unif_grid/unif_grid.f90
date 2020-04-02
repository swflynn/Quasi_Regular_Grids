!=============================================================================80
!                             Uniform Grid
!==============================================================================!
!       Discussion:
!Generate 2D uniform grid for numerical integraion.
!==============================================================================!
!       Modified:
!   30 April 2019
!       Author:
!   Shane Flynn 
!==============================================================================!
program main
implicit none
!==============================================================================!
!               Discussion:
!==============================================================================!
!NP_1D          ==>Number of points across a single dimension 
!Npoints        ==>Total Number of Points
!lower          ==>Lower Bound
!upper          ==>Upper Bound
!points         ==>Single Dimension Coordinates
!x              ==>Grid Coordinates
!==============================================================================!
integer::NP_1D,Npoints,i,j,counter
double precision::lower,upper
double precision,allocatable,dimension(:)::points
double precision,allocatable,dimension(:,:)::x
read(*,*) lower
read(*,*) upper
read(*,*) NP_1D
Npoints=NP_1D**2
!==============================================================================!
allocate(points(NP_1D),x(2,Npoints))
!==============================================================================!
!                           1D Uniform Distribution
!==============================================================================!
do i=1,NP_1D
    points(i)=lower+(i-1.)*(upper-lower)/(NP_1D-1.)
enddo
!==============================================================================!
!                            2D Combinatiorics
!==============================================================================!
counter=1
do i=1,NP_1D
    do j=1,NP_1D
        x(1,counter)=points(i)
        x(2,counter)=points(j)
        counter=counter+1
    enddo
enddo
!==============================================================================!
!                   Write Gridpoints Configuration to file 
!==============================================================================!
open(unit=17,file='grid.dat')
do i=1,Npoints
    write(17,*) x(:,i)
enddo
close(17)
!==============================================================================!
write(*,*) 'Hello Universe!'
end program main
