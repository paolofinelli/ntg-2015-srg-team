      subroutine search (f1,f2,e1,e2,e3,delta,iserch)
!     ===============================================
 
!        this subroutine searches for the zero of a function f(e).
!        first, search where f changes sign, then
!        search for the zero with the Newton-Raphson method.
      
 
      implicit double precision (a-h,o-z)
 
 
 
!        do linear search until f changes sign
 
      if (iserch.eq.0) then
      v1=sign(1.d0,f1)
      v2=sign(1.d0,f2)
      if (v1.eq.v2) then
      e3=e2-sign(delta,f2)
      return
 
!        start newton seach
      else
      e3=e2-f2*(e2-e1)/(f2-f1)
      iserch=1
      endif
 
!        newton search
      else if (iserch.eq.1) then
 
      e3=e2-f2*(e2-e1)/(f2-f1)
 
      endif
 
      return
      end

