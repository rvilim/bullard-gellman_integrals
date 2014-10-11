! This program efficiently calculates the Adams-Gaunt and Elsasser integrals
! for a very large number (10^9) of L's. For mmore information about these integrals
! see http://www.jstor.org/stable/91568 (Bullard and Gellman, 1954)

program sphharmint
  use mod_geoms
  use integration
  use mod_write
  use mod_selection
  use mod_params

  implicit none

  complex(kind=8), allocatable :: a_ag(:,:),a_el(:,:)
  complex(kind=8), allocatable :: b_ag(:,:),b_el(:,:)
  complex(kind=8), allocatable :: c_ag(:,:),c_el(:,:)
  
  integer L1,m1,L2,m2,L3,m3
  
  call allocatestuff
  call initstuff
  call fftinit

  open(unit=22,file=trim('out'),action="write", status="replace" )

!$OMP PARALLEL PRIVATE(L3,m3,L2,m2,L1,m1,a_ag, b_ag, c_ag, a_el, b_el, c_el)

  allocate(a_ag(0:Lmaxa,0:npmax-1))
  allocate(a_el(0:Lmaxa,0:npmax-1))

  allocate(b_ag(0:Lmaxa,0:npmax-1))
  allocate(b_el(0:Lmaxa,0:npmax-1))

  allocate(c_ag(0:Lmaxa,0:npmax-1))
  allocate(c_el(0:Lmaxa,0:npmax-1))
!$OMP DO SCHEDULE(DYNAMIC, 1)
  do L3=1,Lmax
    do m3=0,L3
      do L2=1,L2max
        write(*,*) "L3=", L3, ", L2=", L2
        do m2=0,L2
          a_ag=integrate(L2,m2,L3,m3,NORM,NORM,INT_ADAMSGAUNT)
          a_el=integrate(L2,m2,L3,m3,NORM,NORM,INT_ELSASSER)

          b_ag=integrate(L2,m2,L3,m3,NORM,CC,INT_ADAMSGAUNT)
          b_el=integrate(L2,m2,L3,m3,NORM,CC,INT_ELSASSER)
         
          c_ag=integrate(L2,m2,L3,m3,CC,NORM,INT_ADAMSGAUNT)
          c_el=integrate(L2,m2,L3,m3,CC,NORM,INT_ELSASSER)
        
          do L1=0,Lmax
            if((triangle(L1,L2,L3)==1).and.(mselection(0,m2,m3)==1)) then
               call writeintegral(L1,0,L2,m2,L3,m3,a_ag(L1,0),b_ag(L1,0),c_ag(L1,0),a_el(L1,0),b_el(L1,0),c_el(L1,0))
            endif
          enddo          
          do m1=1,mmax
            do L1=0,Lmax 
              if((triangle(L1,L2,L3)==1).and.(mselection(m1,m2,m3)==1)) then
                call writeintegral(L1,m1,L2,m2,L3,m3,a_ag(L1,m1),b_ag(L1,m1),c_ag(L1,m1),a_el(L1,m1),b_el(L1,m1),c_el(L1,m1))
              endif
            enddo
          enddo
        enddo
      enddo          
    enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL
  close(22)

end program sphharmint
