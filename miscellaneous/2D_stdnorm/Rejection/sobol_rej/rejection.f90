!=============================================================================80
!                       Gaussian (Sobol + Rejection)
!=============================================================================80
!     Discussion:
!Fortran 90 Implementation of the Rejection method for a pdf as described in
!fortran numerical methods.
!A 2D Gaussian Potential, with uniform enclosing function.
!==============================================================================!
!     Modified:
!17 September 2018
!     Author:
!Shane Flynn
!==============================================================================!
module rej_mod
implicit none
!==============================================================================!
!                             Global Variables
!==============================================================================!
!Dimen     ==>Dimensionality of the system
!==============================================================================!
integer::Dimen
!==============================================================================!
contains
!==============================================================================!
subroutine pdf_eval(x,p)
!==============================================================================!
!2D Gaussian to transform uniform sobol distribution to Gaussian distribution
!Uniform sobol dist. is transformed form [0,1] to [min,max](see sobol_stdnormal)
!==============================================================================!
!x        ==>(Dimen) gridpoint Coordinates (uniform distibution)
!p        ==>Gaussian Gridpoint
!==============================================================================!
implicit none
double precision::x(Dimen),p
p=exp(-(((x(1))**2)/2+((x(2))**2)/2))
end subroutine pdf_eval
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
!Nsobol        ==>Number of Grid Points
!count_acc     ==>Number of Accepted Points
!count_rej     ==>Number of Rejections
!pxy           ==>PDF Evaluation
!accept        ==>Randomly Generated Acceptance Condition
!z             ==>(Dimen) Quasi-Random Points
!x             ==>(Dimen,Nsobol) Accepted Coordinates
!==============================================================================!
integer::Nsobol,count_acc,count_rej,i
integer*8::ii,skip                                        !sobol module needs *8
double precision::pxy,accept
double precision,allocatable::z(:),x(:,:)
!==============================================================================!
!                               Read Input File
!==============================================================================!
read(*,*) Dimen
read(*,*) Nsobol
skip=Nsobol
!==============================================================================!
!                           Allocations/Initialize
!==============================================================================!
allocate(z(Dimen),x(Dimen,Nsobol))
count_acc=0
count_rej=0
z=0d0
x=0d0
!==============================================================================!
!                       Store All Coordinates Generated
!==============================================================================!
open(unit=17,file='all_data.dat')
do while(count_acc<Nsobol)
    call sobol_stdnormal(Dimen,skip,z(:))
    write(17,*) z(:)
    call random_number(accept)
    call pdf_eval(z(:),pxy)
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
do i=1,Nsobol
    write(18,*) x(:,i)
enddo
close(18)
!==============================================================================!
!                           Simulation Results
!==============================================================================!
open(unit=99,file='simulation.dat')
write(99,*) 'Nsobol ', Nsobol
write(99,*) 'Dimensionality', dimen
write(99,*) 'Total Rejections', count_rej
close(99)
write(*,*) 'Hello Universe!'
end program reject
