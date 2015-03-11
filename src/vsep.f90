

      subroutine vsep(p,v0,np)


      implicit none

      integer:: np,i,j

      double precision:: p(np),v0(np,np)
      double precision:: g,lsq,ff(np)     

      double precision:: hbarc

      parameter(hbarc=197.3207d0)

      
      lsq=10d0*10d0

      g=-1d0
 
      ff(:)=dexp(-p(:)*p(:)/lsq)
      
      do i=1,np
      do j=1,np

     
      v0(i,j)=g*ff(i)*ff(j)

      end do

      end do

      return
      end subroutine vsep 
        
     
