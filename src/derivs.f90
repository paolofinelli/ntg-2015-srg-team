
      subroutine derivs(sn,vp,fvp,np,ww,p)

      implicit none

      integer:: i,j,n,np
     
      double precision:: p(np),ww(np),fvp(np,np)
      double precision:: vp(np,np),sum0,sum1

      double precision:: pm,pp,sn
      
      

      do i=1,np
      do j=1,np
 
      pm=(p(i)*p(i)-p(j)*p(j))
      pp=(p(i)*p(i)+p(j)*p(j)) 
     
      sum0=0d0
      sum1=0d0
     
      do n=1,np

      sum0=sum0+ww(n)*vp(i,n)*vp(n,j)
      sum1=sum1+ww(n)*p(n)*p(n)*vp(i,n)*vp(n,j)

      end do

      fvp(i,j)=-pm*pm*vp(i,j)+pp*sum0-2d0*sum1       
   
      end do
      end do

      return
      end subroutine derivs
    
 
