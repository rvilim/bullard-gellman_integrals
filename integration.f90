module integration
  use mod_geoms
  use mod_params
  use sphericalharmonics
  
  implicit none
  
CONTAINS
  function adamsgaunt(L2,m2,L3,m3,conj2,conj3) 
  
! Here conj2 and conj3 are whether the term should be
! complex conjugated or not. The constants CC, and NORM
! have been defined in mod_params for this purpose

    integer L2,m2,L3,m3,conj2,conj3
    
    complex(kind=8) adamsgaunt(npmax,ntmax)
    complex(kind=8), allocatable :: Y2(:,:)
    complex(kind=8), allocatable :: Y3(:,:)
      
    allocate(Y2(npmax,ntmax))
    allocate(Y3(npmax,ntmax))
  
    Y2=sphharm(L2,m2)
    Y3=sphharm(L3,m3)
    
    if(conj2 == CC) then
      Y2=conjg(Y2)
    endif
    
    if(conj3 == CC) then
      Y3=conjg(Y3)
    endif
        
    adamsgaunt=Y2*Y3
    
  end function adamsgaunt

  function elsasser(L2,m2,L3,m3,conj2,conj3)
  
! Here conj2 and conj3 are whether the term should be
! complex conjugated or not. The constants CC, and NORM
! have been defined in mod_params for this purpose

    integer L2,m2,L3,m3,conj2,conj3
    integer i
    
    complex(kind=8) elsasser(npmax,ntmax)
    complex(kind=8), allocatable :: dthY2(:,:)
    complex(kind=8), allocatable :: dphY2(:,:)
    complex(kind=8), allocatable :: dthY3(:,:)
    complex(kind=8), allocatable :: dphY3(:,:)
      
    allocate(dthY2(npmax,ntmax))
    allocate(dthY3(npmax,ntmax))
    
    allocate(dphY2(npmax,ntmax))
    allocate(dphY3(npmax,ntmax))
    
    dthY2=dthsphharm(L2,m2)
    dphY2=dphsphharm(L2,m2)
    
    dthY3=dthsphharm(L3,m3)
    dphY3=dphsphharm(L3,m3)
    
    if(conj2 == CC) then
      dthY2=conjg(dthY2)
      dphY2=conjg(dphY2)
    endif
    
    if(conj3 == CC) then
      dthY3=conjg(dthY3)
      dphY3=conjg(dphY3)
    endif
        
    elsasser=dphY2*dthY3-dthY2*dphY3
    
! The elsasser integral is defined as dth Y1 dphY2-dphY1 dth Y2 integrated over
! 0 to 2Pi (phi) and 0 to Pi (theta). The sin theta should _not_ be included in dtheta. The
! spherical harmonic transform implicitly includes the sin theta so I am dividing it out in the
! integrand

    do i=1,npmax
      elsasser(i,:)=elsasser(i,:)/sqrt(1.0-th(:)**2)
    enddo
    
  end function elsasser
  
  function integrate(L2,m2,L3,m3,conj2,conj3, integral) 
    
    integer L2,m2,L3,m3
    integer conj2,conj3
    complex(kind=8) a(npmax,ntmax)
    complex(kind=8) integrate(0:Lmaxa,0:npmax-1)
    integer integral
    
    character(LEN=30) :: Format
    
    Format="(6I4, 2X, ES14.7, 2X, ES14.7)"
    
    if(integral==INT_ELSASSER) then
      a=elsasser(L2,m2,L3,m3,conj2,conj3)
    elseif(integral==INT_ADAMSGAUNT) then
      a=adamsgaunt(L2,m2,L3,m3,conj2,conj3)
    else
      write(*,*) "Error: I can't do that integral"
    endif

    call zfspht_cmplx(a,aslg,gauwt,table,Lmaxa,mmaxa,ntmax,npmax,integrate)
    
    
  end function integrate

end module integration
