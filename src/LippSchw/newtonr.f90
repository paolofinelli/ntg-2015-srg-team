      subroutine newtonr (f1,f2,e1,e2,e3,delta)
!     ========================================
!
!        this subroutine searches for the zero of a function f(e).
!        first, search where f changes sign, then
!        search for the zero with the Newton-Raphson method.
!
!
      implicit none

      double precision f1,f2,e1,e2,e3,delta
      double precision v1,v2 
 
      integer iv
      
      save iv

      data iv/0/
!
!        do linear search until f changes sign
!
      if (iv.eq.0) then
      v1=sign(1.d0,f1)
      v2=sign(1.d0,f2)
 
      if (v1.eq.v2) then
      e3=e2+delta
      return
!        start newton seach
      else 
      e3=e2-f2*(e2-e1)/(f2-f1)
      iv=1
      endif
 
!        newton search
      else if (iv.eq.1) then 
 
      e3=e2-f2*(e2-e1)/(f2-f1)
 
      endif
 
      return
      end
