subroutine aslegend(p,z,Lmax,mmax,inorm)
! ************************************************************************
! 
!   Evaluates normalized associated Legendre polynomial P(L,m) as function of 
!
!      z = cos(theta)            
!
!   up to L=LMAX, m=MMAX using recurrence relation starting with
!   P(m,m) and then increasing L keeping m fixed.
!
!   The normalization is:
!
!   for INORM = 1,
!
!        (Y(L,m)*,Y(L',m')) = 4 pi delta_{L L'} delta_{m m'},
!
!   for INORM = 2,
!
!        (Y(L,m)*,Y(L',m')) = delta_{L L'} delta_{m m'},    
!
!   where 
!
!     Y(L,m) = P(L,m) exp^{i*m*phi},
!
!   which is incorporated into the recurrence relation:
!
!     P(L,m) = z P(L-1,m) \sqrt[(2L+1)(2L-1)/(L+m)(L-m)] - 
!        p(L-2,m) \sqrt[(2L+1)(L+m-1)(L-m-1)/(2L-3)
!        (L+m)(L-m)].
!
!   In the subroutine,
!
!     p(L,m) = P(L,m).
!
!   Routine is modified from the Numerical Recp. subroutine.  The 
!   method is stable in single and double precision to L,m = 511.
!   W.Kuang, 17th, Aug. 1994.
!
!   For the spherical transform developed by W.Kuang, orthonomal
!   spherical harmonics (i.e. INORM=2) is necessary.
!
!   This subroutine is for SUN workstations.
!
! ************************************************************************
  implicit none

  integer inorm,L,Lmax,m,mmax
  real (kind=8) fac,fden,fnum,f1,f2,pi,plm,pmm,pm1,pm2,sin2,z,one,sign
  real (kind=8) p(0:Lmax,0:mmax)

  one = 1.0
  if (Lmax.lt.0 .or. mmax.gt.Lmax .or. abs(z).gt.one) then
    write(*,*) 'bad arguments'
  endif

  if (inorm.lt.0.5 .or. inorm.gt.2.5) then
    write(*,*) 'inorm incorrect:'
    write(*,*) 'inorm = 1 for Full normalisation'
    write(*,*) 'inorm = 2 for orthonormal spherical harmonics'

     stop
  endif
        
  pm2 = one
  p(0,0) = one

  if (Lmax .eq. 0) go to 25

  pm1 = sqrt(3.0*one)*z
  p(1,0) = pm1
  do L = 2,Lmax
    f1 = sqrt(one*(2*L+1)*(2*L-1))
    f2 = (L-1)*sqrt(one*(2*L+1)/(one*(2*L-3)))
    plm = (f1*z*pm1-f2*pm2)/L
    p(L,0) = plm
    pm2 = pm1
    pm1 = plm
  enddo

  if (mmax .eq. 0) go to 25

  !       Evaluating P(L,m) for m > 0

  pmm = one
  sin2 = (one-z)*(one+z)
  fnum = one
  fden = 0.0
  sign = one

  do m = 1,mmax

!----------Evaluating P(m,m) 

    sign = -sign
    fnum = fnum+2.0
    fden = fden+2.0
    pmm = pmm*sin2*fnum/fden
    pm2 = sign*sqrt(pmm)
    p(m,m) = pm2

    if (m .eq. Lmax) goto 25
!----------Evaluating P(m+1,m)
    pm1 = z*pm2*sqrt(one*(2*m+3))
    p(m+1,m)= pm1

!----------Evaluating P(L,m) for L = m+2,...,Lmax

    if (m .lt. (Lmax-1)) then
      do L = m+2,Lmax
        f1 = sqrt(one*(2*L+1)*(2*L-1)/(one*(L+m)*(L-m)))
        f2 = sqrt(one*(2*L+1)*(L-m-1)*(L+m-1)/(one*(2*L-3)*(L+m)*(L-m)))
        plm = z*f1*pm1-f2*pm2
        p(L,m) = plm
        pm2 = pm1
        pm1 = plm
      enddo
    endif
  enddo

  25 continue

  !       Choice of normalization

  if (inorm .eq. 2) then
    pi = 4.0*atan(one)
    fac = 1.0/sqrt(4.0*pi)
    do m = 0,mmax
      do L = m,Lmax
        p(L,m) = p(L,m)*fac
      enddo
    enddo
  endif

  return
end subroutine aslegend
             
