module mod_selection

  implicit none
    
CONTAINS

  function triangle(L1,L2,L3) result(nonzero)
    ! This function tests the triangle inequality, L1, L2, and L3 must be able
    ! to form the sides of a triangle or the integral (Both Elsasser and Adams-
    ! Gaunt is zero. If the integral is zero (L1,L2 and L3 cannot form a triangle)
    ! nonzero is 0, otherwise it is one.

    integer, intent(in) :: L1,L2,L3
    integer :: nonzero

    nonzero=0

    if((abs(L1-L2)<L3).and.(L3<(L1+L2))) then
       nonzero=1
    endif
  end function triangle


  function mselection(m1,m2,m3) result(nonzero)
    ! This function tests that one of m1 +- m2 +- m3 =0. If this is not true then
    ! both the Adams-Gaunt and Elsasser integrals are zero, if m1+-m2+-m3=0 then 
    ! nonzero=1 on exit, otherwise nonzero=0

    integer, intent(in) :: m1,m2,m3
    integer :: nonzero

    nonzero=0

    if((m1+m2-m3==0).or.(m1-m2+m3==0).or.(m1-m2-m3==0)) then
       nonzero=1
    endif
  end function mselection

  function odd(L1,L2,L3) result(nonzero)
    ! This function tests whether L1+L2+L3 is odd, if it is the Elsasser integral is non-zero
    ! If L1+L2+L3 is odd then nonzero=1 on exit, otherwise nonzero=0

    integer, intent(in) :: L1,L2,L3
    integer :: nonzero

    nonzero=0

    if(mod(L1+L2+L3,2)==1) then
      nonzero=1 
    endif
  end function odd 

  function even(L1,L2,L3) result(nonzero)
    ! This function tests whether L1+L2+L3 is even, if it is the Elsasser integral is non-zero
    ! If L1+L2+L3 is even then nonzero=1 on exit, otherwise nonzero=0

    integer, intent(in) :: L1,L2,L3
    integer :: nonzero

    nonzero=0

    if(mod(L1+L2+L3,2)==0) then
      nonzero=1 
    endif
  end function even

  function nonidentical(L1,L2,L3,m1,m2,m3) result(nonzero)
    ! This function ensures that no two of the harmonics are identical. If any are identical
    ! the Elsasser integral is zero. If none are identical then nonzero=1.

    integer, intent(in) :: L1,L2,L3,m1,m2,m3
    integer :: nonzero

    nonzero=1

    if( ((L1==L2).and.(m1==m2)) .or. ((L1==L3).and.(m1==m3)) .or. ((L2==L3).and.(m2==m3)) ) then 
       nonzero=0
    endif
  end function nonidentical

  subroutine mselectionzero(m1,m2,m3,a_ag,b_ag,c_ag,a_el,b_el,c_el)
        integer m1,m2,m3

        complex(kind=8) :: a_ag,a_el
        complex(kind=8) :: b_ag,b_el
        complex(kind=8) :: c_ag,c_el
       
        if((-m1+m2+m3)/=0) then
           a_ag=0.0
           a_el=0.0
        endif

        if((-m1+m2-m3)/=0) then
           b_ag=0.0
           b_el=0.0
        endif

        if((-m1-m2+m3)/=0) then
           c_ag=0.0
           c_el=0.0
        endif
   end subroutine mselectionzero

end module mod_selection
