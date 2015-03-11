
      program ode_solver


!     solves ode

      implicit none 

      integer n1,np,n,i,j,nmax,npr,nm

      parameter(np=70,n1=50,nmax=2000)

      double precision:: hbarc
 
      parameter(hbarc=197.3207d0)

      double precision:: v0(np,np),vn(np,np),dvds0(np,np),x(np) 
      double precision:: sn,ds,p(np),ww(np),w(np),cc,h,hpi
 
      double precision:: vi(np,np),p1,p2,p3

 
      character*90:: filename(nmax)
  

      open(7,file='vspot.out') 


!     initialize
!     ==========

      h=1d-5
      nm=1000
      npr=100

!     potential grid
!     ==============

!      call gauleg(0.d0,1.d0,x,w,np)

!      cc=1.d0
!      hpi=0.5d0*dacos(-1d0)
      
!      do i=1,np

!      p(i)=cc*dtan(hpi*x(i))

!      ww(i)=hpi*cc/(dcos(hpi*x(i))*dcos(hpi*x(i)))*w(i)

!      end do

      p1=.5d0
      p2=1.5d0
      p3=20d0

      call trns(n1,np-n1,np,p1,p2,p3,p,ww)


!     compute potential on grid
!     =========================


      call vsep(p,v0,np)

      vi=v0

     
      sn=0

      do n=1,nm


!     compute deriv at initial s=0
!     ----------------------------

      call derivs(sn,v0,dvds0,np,ww,p)
         

!     find vn at sn=sn+h
!     -------------------


      call rk4m(v0,dvds0,np,sn,h,vn,ww,p)
 
      do i=1,np

      if(mod(n,npr).eq.0.or.n.eq.0)write(7,*) p(i),v0(i,45)
      end do

      write(7,*)
      sn=sn+h
   
       
    
!     print potential
!     ---------------

      if(mod(n,npr).eq.0.or.n.eq.1) then

     
      write(filename(n), '("vkpot",i1,".out")')n
      if(n.gt.99) write(filename(n), '("vkpot",i3,".out")')n
      if(n.gt.999) write(filename(n), '("vkpot",i4,".out")')n

      open(9,file=filename(n))


      do i=1,np

      write(9,*)(v0(i,j),j=1,np)

      end do

      close(9)        
      end if

   
      v0=vn

      end do

      open(10,file= 'k1.out')
      open(11,file='k0.out')

      do i=1,np
      write (10,1000) p(i)
      write(11,1000) p(10)
      end do     

1000 format(1x,300f18.7)


      end program ode_solver

  
