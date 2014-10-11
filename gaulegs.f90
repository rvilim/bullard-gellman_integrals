subroutine gauleg(x1,x2,root,wt,L)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!   This subroutine calculates the assembly points in the colatitude
!   and the corresponding Gaussian weights on the points.  It is
!   modified from NUM.RECP. subroutines.  W.Kuang 15/08/94
! 
! ------------------------------------------------------------------------
! 
!   The assembly points ROOT(L) are the L zeros of the Legendre
!   polynomial P_L(x) [x = cos(th)].  They are symmetric about
!   (x1+x2)/2 and are obtained via Newton method.  The Gaussian
!   weight WT(L) at the assembly points are defined as
! 
!     WT(i) = 2/(1-x_i^2)[P'_L(x_i)]^2.
! 
! ------------------------------------------------------------------------
! 
!   For the assembly points in the colatitude TH (in stead of x), see
!   GAULEG1.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  implicit none

  integer L,i,j,m
  real (kind=8) x1,x2,eps,pi,p1,p2,p3,pp,xl,xm,z,z1,one
  real (kind=8), dimension(L) :: root,wt

  parameter (eps=1.e-15)

  one = 1.0
  m = (L+1)/2
  xm = 0.5*(x2+x1)
  xl = 0.5*(x2-x1)
  pi = 4.0*atan(one)

  do i = 1,m
    !----------initial guess of Z_i
    z = cos(pi*(i-0.25)/(L+0.5))
    !----------Employing Newton method to obtain Z_i

1    continue

    p1 = 1.0
    p2 = 0.0
    do j = 1,L
       p3= p2
       p2= p1
       p1= ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/(1.0*j)
    enddo

!----------Obtaining the derivative P_L'(z) by the values of
!----------P_L(z) (p1) and P_[L-1] (p2).

    pp = L*(z*p1-p2)/(z*z-1.0)
    z1 = z
    z = z1-p1/pp
    if (dabs(z-z1) .gt. eps) goto 1

    root(i) = xm-xl*z
    root(L+1-i) = xm+xl*z
    wt(i) = 2.0*xl/((1.0-z*z)*pp*pp)
    wt(L+1-i) = wt(i)

  enddo

  z1 = sqrt(2.0*pi)

  do i = 1,L
    wt(i)= z1*wt(i)
  enddo
 
return
end subroutine gauleg

