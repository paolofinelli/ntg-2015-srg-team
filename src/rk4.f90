! TODO: Ask Linda about this file. It is not used at all.


      subroutine rk4(y,dydx,n,x,h,yout,derivs)

!        from numerical recipes
!        Given values for the variable y and their derivatives dxdy known at x, use
!        the fourth-order Runge-Kutta method to andvance the solution over an interval h
!        and return the incremented variables as yout, 
!        which need not be a distinct array from y.
!        The user supplies the subroutine derivs(x,y,dydx), 
!        which returns the derivatives dydx at x.
!        n is the number of equations

      implicit none

      external derivs

      integer nmax
      parameter (nmax=50)

      integer n, i
      double precision  h,x,dydx(n),y(n),yout(n)
      double precision  h6,hh,xh,dym(nmax),dyt(nmax),yt(nmax)

      hh=h*0.5
      h6=h/6.
      xh=x+hh
      do i=1,n
        yt(i)=y(i)+hh*dydx(i)
      end do

      call derivs(xh,yt,dyt)

      do  i=1,n
        yt(i)=y(i)+hh*dyt(i)
      end do

      call derivs(xh,yt,dym)

      do  i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
      end do

      call derivs(x+h,yt,dyt)

      do i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
      end do

      return
      end
!       Copr. 1986-92 Numerical Recipes Software 2721[V3.
