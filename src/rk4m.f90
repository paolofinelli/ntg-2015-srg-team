      subroutine rk4m(y,dydx,n,x,h,yout,ww,p)

!        from numerical recipes
!        Given values for the variable y and their derivatives dxdy known at x, use
!        the fourth-order Runge-Kutta method to andvance the solution over an interval h
!        and return the incremented variables as yout, 
!        which need not be a distinct array from y.
!        The user supplies the subroutine derivs(x,y,dydx), 
!        which returns the derivatives dydx at x.
!        n is the number of equations

!        Modified by L.H.

      implicit none

      external derivs

      integer nmax
!     parameter (nmax=50)

      integer n, i,j
      double precision  h,x,dydx(n,n),y(n,n),yout(n,n)
      double precision  h6,hh,xh,dym(n,n),dyt(n,n),yt(n,n)
      double precision  ww(n),p(n)      

      hh=h*0.5
      h6=h/6.
      xh=x+hh
      

      do i=1,n
      do j=1,n  

       yt(i,j)=y(i,j)+hh*dydx(i,j)

      end do
      end do

     call derivs(xh,yt,dyt,n,ww,p)

      do  i=1,n
      do  j=1,n

        yt(i,j)=y(i,j)+hh*dyt(i,j)

      end do
      end do

      
      call derivs(xh,yt,dym,n,ww,p)
 
      do  i=1,n
      do  j=1,n

        yt(i,j)=y(i,j)+h*dym(i,j)
        dym(i,j)=dyt(i,j)+dym(i,j)

      end do
      end do

      call derivs(x+h,yt,dyt,n,ww,p)
 
      do i=1,n
      do j=1,n
 
        yout(i,j)=y(i,j)+h6*(dydx(i,j)+dyt(i,j)+2.*dym(i,j))
      end do
      end do

      return
      end
!       Copr. 1986-92 Numerical Recipes Software 2721[V3.
