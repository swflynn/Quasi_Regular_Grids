!=============================================================================80
!                           Gaussian Random Grid                               !
!==============================================================================!
!       Discussion:
!Generate a d-dimensional Multivariate Gaussian Distribution of points using the
!Fortran random number generator and the Beasley-Springer-Moro transformation.
!==============================================================================!
!       Modified:                                                              
!   16 April 2019                                                              
!       Author:                                                                
!   Shane Flynn                                                                
!==============================================================================!
program main
!==============================================================================!
!               Discussion:
!==============================================================================!
!d              ==>Grid-Point Dimensionality (x^i=x^i_1,x^i_2,..,x^i_d)
!Npoints        ==>Number of Gridpoints 
!x              ==>(d,Npoints) Grid-Point Coordinates
!==============================================================================!
implicit none
integer::d,Npoints,i
double precision,allocatable,dimension(:,:)::x
!==============================================================================!
!                           Read Input Data File                               !
!==============================================================================!
read(*,*) Npoints
read(*,*) d
!==============================================================================!
!                               Allocations
!==============================================================================!
allocate(x(d,Npoints))
!==============================================================================!
!                              Generate Grid 
!==============================================================================!
open(unit=17,file='grid.dat')
do i=1,Npoints
    call rand_stdnormal(d,x(:,i))
    write(17,*) x(:,i)
enddo
close(17)
!==============================================================================!
!                               Output File                                    !
!==============================================================================!
open(99,file='simulation.dat')
write(99,*) 'particle dimensionality ==> ', d
write(99,*) 'Number of Gridpoints ==> ', Npoints
close(99)
write(*,*) 'Hello Universe!'
end program main
