module mod_params
  implicit none

  INCLUDE 'fftw3.f'
  
  integer Lmax,L2max,mmax,CC,NORM,INT_ELSASSER,INT_ADAMSGAUNT
  real(kind=8) min
  real(kind=8) pi
  
  parameter(Lmax=150)
  parameter(mmax=150)

  parameter(L2max=15)

  parameter (pi=3.14159265359)
  
  parameter(CC=-1)
  parameter(NORM=1)
  
  parameter(INT_ELSASSER=1)
  parameter(INT_ADAMSGAUNT=2)
  
  parameter(min=1E-14)
end module mod_params
