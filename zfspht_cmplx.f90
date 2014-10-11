subroutine zfspht_cmplx(f,p,gauwt,wfftr,Lmax,mmax,ntmax,npmax,flm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!   This is the complex-to-complex spherical transform: version 4.
!   For Linux PC with ABSOFT IMSL Library.
!   Weijia Kuang, 10/2002
!   Modified by Ryan Vilim 01/2014
! 
! --------------------------------------------------------------------------
! 
!   before this subroutine must call: 
!     GAULEG: providing assembly points in colatitude
!        and the Gaussian weights;
!     ASLEGEND: providing the values of P_l^m at the assembly
!          points;
!     DFFTRI: initializing the array WFFTR for the FFT.
! 
! --------------------------------------------------------------------------
! 
!   f(i,j):   input,  complex values of F in the physical space (at the
!         assembly points);
!   p(l,m,j): input,  the values of P_l^m at the assembly points;
!   gauwt(j): input,  the Gaussian weights in colatitude;
!   flm(l,m): output, complex spectral coefficients of F.
!   wfftr: 1--forward, 2--backward  plan for fftw
!   Lmax:  maximum degree in colatitude (dealiensing);
!   mmax:  maximum degree in longitude (dealiensing);
!   ntmax:  # of assembly points in THETA (ntmax >= 2*Lmax);
!   npmax:  # of assembly points in PHI (npmax >= 2*mmax);
!        plan_forward: plan for the fftw
!
!   N.B.: This is the complex to complex transform, I am using it to calculate
!         integrals of the form (Ylm^*)x Stuff, this will not work in the the 
!         dynamo model. The dynamo model's spherical harmonic transform expects
!         real input, while this expects complex input! - Ryan
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none

  integer Lmax,mmax,ntmax,npmax

  real (kind=8) gauwt(ntmax),p(0:Lmax,0:mmax,ntmax)
  integer (kind=8) wfftr(2)

  complex (kind=8) ft(npmax)
  complex (kind=8) f(npmax,ntmax)
  complex (kind=8) flm(0:Lmax,0:npmax-1)
  complex (kind=8) temp(0:npmax-1,ntmax)
  integer k,L,m
  real (kind=8) wt,wtfac,pi,one

  one = 1.0
  pi = 4.0*atan(one)
  wtfac = sqrt(2.0*pi)/npmax

  flm = 0.0

  do k  = 1,ntmax
    ft = f(:,k)
    call dfftw_execute_dft(wfftr(1),ft,temp(0,k))
  enddo


  do k = 1,ntmax
    wt= gauwt(k)*wtfac       
    do L = 0,Lmax
      flm(L,0)= flm(L,0)+p(L,0,k)*real(temp(0,k))*wt
    enddo
    do m = 1,mmax
      do L = m,Lmax
        flm(L,m)= flm(L,m)+p(L,m,k)*wt*temp(m,k)
        flm(L,npmax-m)= flm(L,npmax-m)+p(L,m,k)*wt*temp(npmax-m,k)
!         write(*,*) flm(L,m),flm(L,npmax-m)
      enddo
    enddo
  enddo
  
!   do L = 0,Lmax
!     write(*,*) L,0,real(flm(L,0)),aimag(flm(L,0))
!   enddo
!   
!   do m=1,mmax
!     do L=0,Lmax
!       write(*,*) L,m,real(flm(L,m)),aimag(flm(L,m))
!       write(*,*) L,-m,real(flm(L,npmax-m)),aimag(flm(L,npmax-m))
!     enddo
!   enddo
    
 return
 end



