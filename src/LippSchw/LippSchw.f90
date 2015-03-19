!--------------------------------------------------------------------------
!	This program computes the solution of the Integro-differential
!	Equation for the deuteron problem. It also computes the eigenvalues
!	and eigen vectors for the same. Eventually provides the binding energy
!	of the system. The computation is for momentum space.
!	
!------------------------------------------------------------------------------		
!	Code needs lapack and lblas libraries to compute the eigenvalues and 
!	the wave function. In addition it also needs other subroutines.
!------------------------------------------------------------------------------
      program momentum

!      external v1
      implicit none

      external v1,v1sep
      integer :: i,j,ng,true, true1, k, kread, kwrite

      integer:: iemax
      double precision :: mu,hc,pi,Ea,Enew, tolnr
      double precision :: v1,norm, dE, E, dk, v1sep

      parameter (ng=60)

      double precision :: x(ng),w(ng) ! points and weights from Gauss-Legendre
      double precision :: p(ng),w1(ng) ! mapping: points and weights to do the integral
      double precision :: h(ng,ng), identity(ng,ng) ! Hamiltonian  and identity matrix
      double precision :: diff(ng), delta,cc
      double precision :: WR(ng),WI(ng),VL(ng,ng),VR(ng,ng),WORK(4*ng), dn(ng)
      integer :: INFO,isearch

      double precision xx(ng)
      double precision fxx(ng)
      integer ie, iemax

      data kread/5/, kwrite/6/
      data iemax/200000/


!	 define constants

      hc=197.32705d0
      pi=dacos(-1.d0)
      Ea =3d0

      tolnr=1.d-10
      mu=1d0!938.926d0
      cc=3.d0
      isearch=0 
!	 get Gauss points and weights
      call gauleg (-1.d0,1.d0,x,w,ng)

!	 Converting the indefinite integral to a finite integral by a transformation
      do i=1,ng
        p(i)=cc*dtan(pi/4.d0*(1+x(i))) 
        w1(i)=cc*pi/4.d0*w(i)/(dcos(pi/4.d0*(1.d0+x(i))))**2
!	Identity matrix
        identity(i,i) = 1.d0
      end do
      open (7, file ="pot.dat", status ="unknown")
 !     write(7,*) "#         p(i)      qprime = 1.0   qprime = 1.0"
      do i = 1, ng  
      write(*,10002) p(i), (v1sep(p(i),1.d0)), (v1sep(p(i),1.5d0))
      end do
      close(7)

      
	k = 0
	dk = 1d-2
!	do k = 1, 10
10	k = k+1
	
	if (k.eq.1) then 
	E = Ea
	elseif(k.eq.2) then
        e=ea+dk
        else 
        ea=e
 	e =enew
	end if
! 	build the Hamiltonian matrix
      do i=1,ng
        do j=1,ng
        h(i,j)=-w1(j)/(E+p(i)**2/mu)*v1sep(p(i),p(j))
!print*, h(i,j)
      end do
      end do
!      if(k.eq.1)print*,h(3,:)

! use dgeev.f to calculate eigenvector and eigenvalues
      call DGEEV('N', 'V', ng, h, ng, WR, WI, VL, ng, VR, ng, WORK, 4*ng, INFO)
!      if(k.eq.1) print*,h(:,20) 
   
!	 find the position of binding energy
       
!      if(k.gt.1)print*,'ddg',ea,e
      true=ng
 
      dn(true)=1d3

      do i=1,ng

        dn(i)=dabs(1d0-wr(i))
        if (dn(i).lt.dn(true)) then
        true=i!	Choose the eigen value and ignore others
        end if
!      print*,dn(i),wr(i),wr(true)
      end do
!	Print binding energy to screen
!      write(*,*) 'eigenvalue ',wr(true),'MeV'


!	Newton Raphson method for finding zeros
	diff(k) = dn(true)
!	print*, wr(true)

	if(k<2) then 
!	k = k+1
	goto 10
	end if
	
	open(10, file = "test.dat", status ="unknown")
	
!        do i=1,ng
        write(10,*)e, diff(k)
!        end do

!        do newton search
!        starting from xinitial

      fxx(1) = diff(k-1)
      fxx(2) = diff(k)
!      print*,k
!      print*,ea-e,diff(1)-diff(2)
!      ie=1
  
! 20   ie=ie+1 
!      if (ie.gt.iemax) stop!

!     write (kwrite,*) ' f(',xx(ie),') = ',fxx(ie),ie
!    
      
      call newtonr (fxx(1),fxx(2),ea,e,enew,dk)
     
      if (abs(fxx(2)).lt.tolnr) then
      print*,''  
      print*,'final result with tol= 1d-10'
      print*,'============================'
      
       write(6,*) '           E_b [MeV]','       delta'
      end if
  

      if(k.eq.2) write(6,*) '           E_b [MeV]','       delta'
      write(*,10002) e,fxx(2)
!      endif

!      Enew = xx(ie+1)


!	write(*,*) k
        if (abs(fxx(2)).lt.tolnr) stop
	if (k>iemax) stop
	go to 10



10000 format(2x,2f20.12) 
10002 format(2x,1f20.12, 3d15.7)
      close(6)
      close(10)
      stop 'output written in psi_10.dat'
      end program

!========================================================================
include 'search.f'
