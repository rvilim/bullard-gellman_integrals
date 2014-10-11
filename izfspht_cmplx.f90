subroutine izfspht_cmplx(flm,p,wfftr,Lmax,mmax,ntmax,npmax,f)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This is the complex-to-real inverse spherical transform: version 4.
! For Linux PC with ABSOFT IMSL Library.
! Weijia Kuang, 10/2002
!
!--------------------------------------------------------------------------
!
! The subroutine calculates
! 	f(ph_i,th_j) = \sum{l,m} f_l^m Y_l^m(th_j,ph_i)
! 	Y_l^m(th,ph) = P_l^m(th) exp(i m ph)
!
!--------------------------------------------------------------------------
!
! before this subroutine must call:
! 	GAULEG: providing assembly points in colatitude
! 		 and the Gaussian weights;
! 	ASLEGEND: providing the values of P_l^m at the assembly
! 		   points;
! 	DFFTRI: initializing the array WFFTR for the FFT.
!
!--------------------------------------------------------------------------
!
! flm(l,m): input,  complex spectral coefficients of F;
! p(l,m,j): input,  the values of P_l^m at the assembly points;
! f(i,j):	  output, real values of F in the physical space (at
! 	  	  assembly points).
! wfftr: 1--forward, 2--backward plan for fftw
! Lmax:	maximum degree in colatitude;
! mmax:	maximum degree in longitude;
! ntmax:	# of assembly points in colatitude; (ntmax >= 2*mmax+1)
! npmax:	# of assembly points in longitude; (npmax >= 2*Lmax+1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none

  integer Lmax,mmax,ntmax,npmax

  real (kind=8) p(0:Lmax,0:mmax,ntmax)
  integer (kind=8) wfftr(2)
  complex (kind=8) flm(0:Lmax,0:npmax-1),f(npmax,ntmax)
  complex (kind=8) temp(0:npmax-1)
  integer i,k,L,m

  complex (kind=8) c1,c2

  f = 0.0
  do k = 1,ntmax
    temp(0)=(0.0,0.0)
    c1 = 0.0
    c2 = 0.0
    do L = 0,Lmax
      c1 = c1+p(L,0,k)*flm(L,0)
    enddo
    temp(0) = cmplx(real(c1),0.0)
    
    do m = 1,mmax
      c1=0.0
      c2=0.0
      do L = m,Lmax
        c1= c1+p(L,m,k)*flm(L,m)
        c2= c2+p(L,m,k)*flm(L,-m)
      enddo
      temp(m)=c1
      temp(-m)=c2
    enddo

    call dfftw_execute_dft(wfftr(2),temp,f(1,k))
    if(k.eq.4) then
      write(*,*) "-------------"
      write(*,*) temp
      write(*,*) "////"
      write(*,*) f(:,k)
      write(*,*) "-------------"
    endif
  enddo
  
end subroutine izfspht_cmplx

