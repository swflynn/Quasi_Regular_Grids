!=============================================================================80
!                       2D Morse DVR and Truncation
!=============================================================================80
!Uniform DVR Grid (with Truncation) for 2D Morse Potential
!Requires the cg.f code (RS subroutine, solves the eigenvalue problem)
!I should go back and replace this with llapack
!I should go back and re-write this code to be in my style
!=============================================================================80
!       Modified:
!   5 May 2019
!       Author:
!   Shane Flynn
!==============================================================================!
PROGRAM QM_2d
!==============================================================================!
!x,y min,max    ==>Domain for 2D_Morse
!Nx,Ny          ==>Number of Grid Points (single dimension)
!Ndvr           ==>Total Number of Grid Points (after truncation)
!==============================================================================!
implicit none
double precision,parameter::xmin=-3.343868823,xmax=18.91310042,Ecut=11.5
double precision,parameter::ymin=-3.715644418,ymax=21.01588302,pi=dacos(-1d0)
integer,parameter::Nx=38,Ny=Nx
integer::i,j,i1,i2,j1,j2,Ndvr,IERR
double precision::x(Nx),y(Ny),dx,dy,Potential,q(2),q1(2),Vpot,T1(Nx,Nx),T2(Ny,Ny)
double precision,allocatable,dimension(:)::E,FV1,Fv2
double precision,allocatable,dimension(:,:)::H,z
!==============================================================================!
!                            Uniform Spacing
!==============================================================================!
open(1,file='Energies.dat')
dx=(xmax-xmin)/(Nx-1)
dy=(ymax-ymin)/(Ny-1)
do i=1,Nx
    x(i)=xmin+(i-1)*dx
enddo
do i=1,Ny
    y(i)=ymin+(i-1)*dy
enddo
open(7,file='grid.dat')
do i=1,Nx
    do j=1,Nx
        if(i.ne.j) then
            T1(i,j)=(-1)**(i-j)/(dx*(i-j))**2
        else
            T1(i,j)=pi**2/6/dy**2
        endif
    enddo
enddo
do i=1,Ny
    do j=1,Ny
    if(i.ne.j) then
        T2(i,j)=(-1)**(i-j)/(dy*(i-j))**2
    else
        T2(i,j)=pi**2/6/dy**2
    endif
    enddo
enddo
!==============================================================================!
!                 Total Number of Grid Points within cutoff contour
!==============================================================================!
Ndvr=0
do i=1,Nx
    q(1)=x(i)
    do j=1,Ny
        q(2)=y(j)
        if(Potential(q)<Ecut) then
            Ndvr=Ndvr+1
            write(7,*) x(i),y(j)
        endif
    enddo
enddo
write(*,*) 'Nx=',Nx,' Ny=',Ny,' dx=',dx,' dy=',dy,' Ndvr=', Ndvr
allocate(H(Ndvr,Ndvr),z(Ndvr,Ndvr),E(Ndvr),Fv1(Ndvr),Fv2(Ndvr))
!==============================================================================!
!                            Compute Hamiltonain
!==============================================================================!
i=0
H=0d0
do i1=1,Nx
    q(1)=x(i1)
    do i2=1,Ny
        q(2)=y(i2)
        Vpot=Potential(q)
        if(Vpot<Ecut) then
            i=i+1
            H(i,i)=H(i,i)+Vpot
            j=0
            do j1=1,Nx
                q1(1)=x(j1)
                do j2=1,Ny
                    q1(2)=y(j2)
                    if(Potential(q1)<Ecut) then
                        j=j+1
                        if(i1==j1) then   ! y-kinetic energy
                            if(i2.ne.j2) then
                                H(i,j)=H(i,j)+(-1)**(i2-j2)/(dy*(i2-j2))**2
                            else
                                H(i,j)=H(i,j)+pi**2/6/dy**2
                            endif
                        endif
                        if(i2==j2) then   ! x-kinetic energy
                            if(i1.ne.j1) then
                                H(i,j)=H(i,j)+(-1)**(i1-j1)/(dx*(i1-j1))**2
                            else
                                H(i,j)=H(i,j)+pi**2/6/dx**2
                            endif
                        endif
                    endif
                enddo
            enddo
        endif
    enddo
enddo
!==============================================================================!
!                         Compute Vibrational Spectra
!==============================================================================!
call RS(Ndvr,Ndvr,H,E,1,z,FV1,FV2,IERR)
write(1,*) '#   Nx=',Nx,' Ny=',Ny,' Ecut=', Ecut,'  Ndvr=', Ndvr
do i=1,122
    write(1,*) E(i)
enddo
!==============================================================================!
end program QM_2d
!==============================================================================!
function Potential(x)
!==============================================================================!
!Hard-coded Morse Potential Energy
!==============================================================================!
!x              ==>(d) ith particles coordinate x^i_1,..,x^i_d
!Potential      ==>evaluate V(x)
!D              ==>Parameter for Morse Potential
!a1,a2          ==>Parameters for Morse Potential
!==============================================================================!
implicit none
double precision,parameter::D=12.,a1=0.2041241,a2=0.18371169
double precision::x(2),Potential
Potential=D*((exp(-a1*x(1))-1)**2+(exp(-a2*x(2))-1)**2)
return
end function Potential
