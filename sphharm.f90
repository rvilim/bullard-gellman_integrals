module sphericalharmonics
  use mod_geoms
  use mod_params

  implicit none
  
CONTAINS
  function sphharm(L,m)
! This function calculates Ylm in real space for a given
! L and m. It does this by multiplying our associated Legendre
! polynomial by the exponential phi factor.

    complex(kind=8) sphharm(npmax,ntmax)
    complex(kind=8) uim

    integer(kind=8) i
    integer L,m
  
    uim=cmplx(0.0,1.0*m)
  
    do i=1,ntmax
      sphharm(:,i)=exp(uim*ph(:))*aslg(L,m,i)
    enddo
  end function sphharm
  
  function dphsphharm(L,m)
! This function calculates d/dph YLm in real space. It really just calls
! the sphharm then multiplies it by i*m. Phi derivatives are
! really easy ...

     complex(kind=8) dphsphharm(npmax,ntmax)
     complex(kind=8) uim
     
     integer L,m
  
     uim=cmplx(0.0,1.0*m)

     dphsphharm=uim*sphharm(L,m)
    
  end function dphsphharm
  
  function dthsphharm(L,m)

! This function calculates the theta derivative of a spherical harmonic
! mode Ylm in real space. It does this with an identitity which relates the sin(th) * d/dth Ylm
! to two spherical harmonics modes adjacent in L (YL+1,m and YL-1,m)

! At the end we divide everything by sin(th) because we want d/dth Ylm not sinth d/dth Ylm

    complex(kind=8) dthsphharm(npmax,ntmax)
    complex(kind=8) sphharm1,sphharm2
    complex(kind=8) uim

    integer(kind=8) i
    integer L,m
  
    real(kind=8) c1,c2
    
    uim=cmplx(0.0,1.0*m)
  
    sphharm1=0.0
    sphharm2=0.0
    
! The following constants are from euqation 8 on page 147 of 
! The Quantum Theory of Angular Momentum by Varshalovich.
! They are the constants for calculating sin(th) d/dth Ylm

    c1=L*sqrt(1.0*(L-m+1)*(L+m+1)/(1.0*(2*L+1)*(2*L+3)))
    c2=(L+1)*sqrt(1.0*(L-m)*(L+m)/(1.0*(2*L-1)*(2*L+1)))
    
! Here I am calculating two spherical harmonics instead of the usual one
! this is because the sin(th) d/dth Ylm identity splits the mode into a L+1
! and an L-1 mode (look up the Varshalovich refrence). I have some error checking
! to make sure that L<Lmaxa (the maximum Legendre polynomial that we calculated) 
! and to make sure that we aren't trying to get Legendre polynomials that have
! L<0 or m>L. In all these cases the sphharm1 or 2 should be zero.

    do i=1,ntmax
      if(L<Lmaxa) then
        sphharm1=aslg(L+1,m,i)*c1
      endif
      
      if((L>0).AND.(L>m)) then
        sphharm2=aslg(L-1,m,i)*c2
      endif
      
      dthsphharm(:,i)=exp(uim*ph(:))*(sphharm1-sphharm2)/sqrt(1.0-th(i)**2)
    enddo
  end function dthsphharm
  
end module sphericalharmonics

