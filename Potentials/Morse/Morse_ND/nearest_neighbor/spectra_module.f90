!=================20==================40==================60==================80
!                            nD-Morse EigenSpectra
!==============================================================================!
!Compute the EigenSpectra for the nD Morse Potential
!All equations are for nD case, need to define a specific potential subroutine
!generate alpha (Gaussian Widths) based on nearest neghbor
!This code needs the gen_hermite_rule.f90 code (Gauss-Hermite-Quadriture)
!This code uses LLAPACK to solve the generalized eigenvalue problem
!==============================================================================!
!       Modified:
!   11 March 2020
!       Author:
!   Shane Flynn
!==============================================================================!
module spectra_mod
implicit none
!==============================================================================!
contains
!==============================================================================!
subroutine read_grid(d,NG,x)
implicit none
integer::d,NG,i
double precision::x(d,NG)
open(17,File='grid.dat')
do i=1,NG
  read(17,*) x(:,i)
enddo
close(17)
end subroutine read_grid
!==============================================================================!
subroutine gaussian_widths(d,NG,alpha0,alpha,x)
!==============================================================================!
!use nearest neighbor to generate gaussian widths (assumes QRG grid)
!==============================================================================!
implicit none
integer::d,NG,i,j
double precision::alpha0,alpha(NG),r2,x(d,NG)
do i=1,NG
  alpha(i)=1d20                          !large initial distance for placeholder
  do j=1,NG
    if(j.ne.i) then
      r2=sum((x(:,i)-x(:,j))**2)                    !distance between gridpoints
      if(r2<alpha(i)) alpha(i)=r2
    endif
  enddo
  alpha(i)=alpha0/alpha(i)
enddo
open(unit=18,file='alphas.dat')
do i=1,NG
  write(18,*) alpha(i)
enddo
close(18)
end subroutine gaussian_widths
!==============================================================================!
subroutine overlap_elements(d,alpha_i,alpha_j,Sij,xi,xj)
!==============================================================================!
!compute i-j element of the overlap matrix (analytic expressions)
!==============================================================================!
integer::d,i,j
double precision::alpha_i,alpha_j,Sij,xi(d),xj(d),aij,r2
aij=alpha_i*alpha_j/(alpha_i+alpha_j)
r2=sum((xi(:)-xj(:))**2)
Sij=(2*sqrt(alpha_i*alpha_j)/(alpha_i+alpha_j))**(0.5*d)*exp(-aij*r2)
end subroutine overlap_elements
!==============================================================================!
subroutine overlap_matrix(d,NG,alpha,Smat,x)
!==============================================================================!
!compute the overlap matix
!==============================================================================!
integer::d,NG,i,j
double precision::alpha(NG),Smat(NG,NG),x(d,NG),aij
do i=1,NG
  do j=i,NG
    call overlap_elements(d,alpha(i),alpha(j),Smat(i,j),x(:,i),x(:,j))
    Smat(j,i)=Smat(i,j)
  enddo
enddo
!write(*,*) Smat
end subroutine overlap_matrix
!==============================================================================!
subroutine diagonalize_overlap(NG,Smat,eigenvalues,lwork)
!==============================================================================!
!Check to see if S is positive definite
!==============================================================================!
implicit none
integer::NG,info,lwork,i
double precision::eigenvalues(NG),work(max(1,lwork)),Smat(NG,NG)
call dsyev('v','u',NG,Smat,NG,eigenvalues,work,Lwork,info)
write(*,*) 'Info (Initial Overlap Matrix) ==> ', info
write(*,*) 'Recriprocal Condition Number ==>', eigenvalues(1)/eigenvalues(NG)
open(unit=19,file='overlap.dat')
do i=1,NG
  write(19,*) eigenvalues(i)
enddo
close(19)
end subroutine diagonalize_overlap
!==============================================================================!
subroutine kinetic_elements(d,alpha_i,alpha_j,Tij,xi,xj)
!==============================================================================!
!compute i-j element of the overlap matrix (analytic expressions)
!==============================================================================!
implicit none
integer::d,i,j
double precision::alpha_i,alpha_j,Sij,xi(d),xj(d),aij,Tij,r2
aij=alpha_i*alpha_j/(alpha_i+alpha_j)
r2=sum((xi(:)-xj(:))**2)
Tij=aij*(d-2.*aij*r2)
!write(*,*) 'Working kinetic elements subroutine', Tij
end subroutine kinetic_elements
!==============================================================================!
function V(d,x)
!==============================================================================!
!Hard-coded Morse Potential Energy
!==============================================================================!
!x              ==>(d) ith atoms coordinates
!V              ==>evaluate V(x_i)
!D_morse        ==>Parameter for Morse Potential
!omega(d)       ==>Parameter for Morse Potential
!==============================================================================!
implicit none
integer::d
double precision::x(d),V
double precision,parameter::omega(3)=(/0.2041241,0.18371169,0.16329928/)
double precision,parameter::D_morse=12.
!write(*,*) 'V test omega ', omega(:)
!write(*,*) 'V test dmorse', D_morse
!write(*,*) 'V test x', x(d)
V=D_morse*sum((exp(-omega(:)*x(:))-1.)**2 )
!write(*,*) 'V tester', V
end function V
!==============================================================================!
subroutine get_hamiltonian(d,NG,GH_order,x,alpha,Smat,Hmat)
!==============================================================================!
!Compute the Hamiltonian and solve the generalized eigenvalue problem
!Use Gauss Hermit quadrature to evaluate the potential matrix cgqf
!z(GH-order) w(GH-order) --- quadrature points and weights
!requires llapack and gen_hermite_rule.f90
!==============================================================================!
implicit none
integer::d,NG,GH_order,i,j,k,ll
double precision,parameter::pi=4.*atan(1d0)
double precision::w(GH_order),rr(d),Smat(NG,NG),Hmat(NG,NG),r2,Vij,x_ij(d)
double precision::eigenvalues(NG),z(GH_order),l(d),alpha(NG),x(d,NG)
call cgqf(GH_order,6,0d0,0d0,0d0,1d0,z,w)
w=w/sqrt(pi)
!write(*,*) 'w test', w
!==============================================================================!
!                   Solve Generalized Eigenvalue Problem
!==============================================================================!
do i=1,NG
  do j=i,NG
    call overlap_elements(d,alpha(i),alpha(j),Smat(i,j),x(:,i),x(:,j))
    Smat(j,i)=Smat(i,j)
!    write(*,*) 'Smat test initial ', Smat(i,j)
!==============================================================================!
!                          Kinetic Energy Matrix
!==============================================================================!
    call kinetic_elements(d,alpha(i),alpha(j),Hmat(i,j),x(:,i),x(:,j))
!==============================================================================!
!                         Potential Energy Matrix
!==============================================================================!
    x_ij(:)=(alpha(i)*x(:,i)+alpha(j)*x(:,j))/(alpha(i)+alpha(j))
!    write(*,*) 'xij test', x_ij
    Vij=0d0
    l(:)=1
    do ll=1,GH_order**d
      do k=1,d
        rr(k)=z(l(k))
      enddo
!      write(*,*) 'rr test', rr
!      write(*,*) 'rr test xij', x_ij
      rr=x_ij+rr/sqrt(alpha(i)+alpha(j))
!      write(*,*) 'rr test 2 ', rr
!      write(*,*) 'rr alpha(i)', alpha(i)
!      write(*,*) 'rr alpha(j)', alpha(j)
!      write(*,*) 'rr test sqrt',  sqrt(alpha(i)+alpha(j))
!      write(*,*) 'rr test + sqrt', rr/sqrt(alpha(i)+alpha(j))
      r2=V(d,rr)
!      write(*,*) 'r2 ', r2
      do k=1,d
        r2=r2*w(l(k))
      enddo
      Vij=Vij+r2
!      write(*,*) 'Vij test', Vij
      do k=1,d
        l(k)=mod(l(k),float(GH_order))+1
        if(l(k).ne.1) exit
      enddo
    enddo
!==============================================================================!
!         Hamiltonian = Kinetic + Potential (all terms are scaled by Smat)
!==============================================================================!
!    write(*,*) 'Hmat ' ,  Hmat(i,j)
!    write(*,*) 'i, j' ,  i, j
!    write(*,*) 'Hmat test V' ,  Vij
!    write(*,*) 'Hmat test Sij' , Smat(i,j)
!    write(*,*) 'Hmat test Sji' , Smat(j,i)
    Hmat(i,j)=(Hmat(i,j)+Vij)*Smat(i,j)
    Hmat(j,i)=Hmat(i,j)
  enddo
enddo
!write(*,*) 'Hamiltonian ?', Hmat
write(*,*) 'Working hamiltonian subroutine'
end subroutine get_hamiltonian
!==============================================================================!
subroutine hamiltonian_eigenvalues(NG,Hmat,Smat,eigenvalues,work,lwork)
!==============================================================================!
implicit none
integer::NG,info,lwork,itype,i
double precision::Smat(NG,NG),Hmat(NG,NG),eigenvalues(NG),work(max(1,lwork))
itype=1
call dsygv(itype,'n','u',NG,Hmat,NG,Smat,NG,eigenvalues,work,Lwork,info)
write(*,*) 'info ==> ', info
open(unit=20,file='eigenvalues.dat')
do i=1,NG
  write(20,*) eigenvalues(i)
enddo
close(20)
write(*,*) 'Working hamiltonian eigenvalues subroutine'
end subroutine hamiltonian_eigenvalues
!==============================================================================!
subroutine write_out(d,NG,alpha0,GH_order)
!==============================================================================!
!write output file
!==============================================================================!
implicit none
integer::d,NG,GH_order
double precision::alpha0
open(unit=99,file='out')
write(99,*) 'dimensionality ==> ', d
write(99,*) 'Number of  Gaussians==> ', NG
write(99,*) 'alpha0==>', alpha0
write(99,*) 'GH Order==>', GH_order
close(99)
end subroutine write_out
!==============================================================================!
end module spectra_mod
