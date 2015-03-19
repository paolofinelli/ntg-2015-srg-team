!	One term-separable potential with parameters g = -1 and lambda = 10

!      module potential
      double precision function v1sep(p,pprime)
      implicit none
      double precision :: p, pprime, g, lambda
      parameter (lambda = 10.d0)
      parameter (g =-1.d0)
      

      v1sep =2d0*g*dexp(-(p**2+pprime**2)/lambda**2)

      end      


!========================================================================

      double precision function v1(p1,p1prime)

!	This is the gaussian potential given by 
!	V(p,p') = Sum over i of V_i/(2 pi) exp (-(p-p')^2 (alpha_i)^2/4)
!	To use this potential for either of model 1 or 2, un comment the line 
!	containing that potential.
      implicit none

      integer :: i, j, n
      double precision :: p1, p1prime, alpha, pi,hc,m0 
      double precision :: v11, v12, v21, v22, alpha11, alpha12, alpha21, alpha22

      parameter(m0=938d0, hc=197.3207d0)
      parameter (v11 = 12.d0)
      parameter (v12 = -12.d0)
      parameter (v21 = 0.d0)
      parameter (v22 = -2.d0)
      parameter (alpha11 = 0.2d0/2)
      parameter (alpha12 = 0.8d0/2)
      parameter (alpha21 = 0.d0/4)
      parameter (alpha22 = 0.8d0/4)
!      parameter(pi = 3.14158281)
      pi=4.d0*atan(1.d0)


!	Model 1	
	v1 =1.d0/(2d0*pi)*(v11*dexp(-(p1-p1prime)**2*alpha11**2) + v12*dexp(-(p1-p1prime)**2*alpha12**2))
   
!	Model 2	
!	v1 = -(1/(2*pi))*(exp(-(p1-p1prime)**2*(alpha21**2)) + exp(-(p1-p1prime)**2*(alpha22**2)))

      end


      
