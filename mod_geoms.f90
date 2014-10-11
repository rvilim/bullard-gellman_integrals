module mod_geoms
  
  use mod_params
  
  implicit none
  
  integer Lmaxa, mmaxa
  integer ntmax, npmax
  
  real(kind=8), allocatable :: th(:),ph(:)
  real(kind=8), allocatable :: gauwt(:)
  real(kind=8), allocatable :: aslg(:,:,:)
  integer (kind=8) table(2) 

CONTAINS
  subroutine fftinit
    complex( kind=8),allocatable,dimension(:):: tempfftout,tempfftin
    
    allocate(tempfftin(npmax))
    allocate(tempfftout(npmax))

    call dfftw_plan_dft_1d(table(1),npmax,tempfftin, tempfftout, FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_plan_dft_1d(table(2),npmax,tempfftout, tempfftin, FFTW_BACKWARD, FFTW_ESTIMATE)

  end subroutine fftinit
  
  subroutine allocatestuff
  
    Lmaxa = 3*(Lmax+2)/2
   
    ntmax = Lmaxa+1
    npmax = 3*(mmax+1)+1

    if (mod(npmax,2) == 1) then
      npmax = npmax+1
    endif

    mmaxa = (npmax-1)/2
   
    allocate(th(ntmax))
    allocate(ph(npmax))
    allocate(gauwt(ntmax))
    allocate(aslg(0:Lmaxa,0:mmaxa,ntmax))
  end subroutine allocatestuff
  
  subroutine initstuff
      
    integer(kind=8) i
    real(kind=8) one
    
    th = 0.0
    gauwt = 0.0
    aslg = 0.0
    one=1.0    
    
    call gauleg(-one,one,th,gauwt,ntmax)

    
    do i=1,npmax
      ph(i)=2*pi*(i-1)/npmax
    enddo
    
    do i=1,ntmax
      call aslegend(aslg(0,0,i),th(i),Lmaxa,mmaxa,2)
    enddo
    
  end subroutine initstuff
    
end module mod_geoms
