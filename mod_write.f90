module mod_write

  use mod_params
  use mod_selection
  
  implicit none
    
CONTAINS
  subroutine writeintegral(L1,m1,L2,m2,L3,m3,a_ag,b_ag,c_ag,a_el,b_el,c_el)
    
    complex(kind=8) :: a_ag,a_el
    complex(kind=8) :: b_ag,b_el
    complex(kind=8) :: c_ag,c_el
    
    integer L1,m1,L2,m2,L3,m3

100 format(6(I4.1), 4X , 6(ES17.10))
   
    if( (abs(real(a_ag))>min .OR. &
        abs(real(b_ag))>min .OR. &
        abs(real(c_ag))>min .OR. &
        abs(aimag(a_el))>min .OR. &
        abs(aimag(b_el))>min .OR. &
        abs(aimag(c_el))>min))  then
       
        if(abs(real(a_ag))<min) then
           a_ag=0.0
        endif 
        if(abs(real(b_ag))<min) then
           b_ag=0.0
        endif 
        if(abs(real(c_ag))<min) then
           c_ag=0.0
        endif 
        if(abs(real(a_el))<min) then
           a_el=0.0
        endif 
        if(abs(real(b_el))<min) then
           b_el=0.0
        endif 
        if(abs(real(c_el))<min) then
           c_el=0.0
        endif 

        call mselectionzero(m1,m2,m3,a_ag,b_ag,c_ag,a_el,b_el,c_el)
 
        if ((even(L1,L2,L3)==1).or.((odd(L1,L2,L3)==1).and.(nonidentical(L1,L2,L3,m1,m2,m3)==1))) then
           write(22,100) L1,L2,L3,m1,m2,m3,real(a_ag),real(b_ag),real(c_ag),aimag(a_el),aimag(b_el),aimag(c_el)
        endif
    endif
  end subroutine writeintegral
end module mod_write
