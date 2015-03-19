!Note that Vs has to be evaluated at p1 and p2
SUBROUTINE srg_flow(Lambda,Cutoff,Vs,function)
   
   IMPLICIT NONE
   INTEGER :: i,j,l,k, NN, N, NNN
   PARAMETER(N=10,NN=10,NNN=10)
      
   DOUBLE PRECISION,DIMENSION(N,N) :: function,Vs
   DOUBLE PRECISION ::absc,weights
    DOUBLE PRECISION ::abscc,weightsc
   DOUBLE PRECISION:: Lambda,Vs,test!,h,g
   DIMENSION absc(NN),weights(NN)
   DIMENSION abscc(NNN),weightsc(NNN)
!   Lambda=10.d0
!   g=-1.d0
   CALL gauleg(0.d0,Cutoff,absc,weights,NN)
   CALL gauleg(0.d0,Lambda,abscc,weightsc,NNN)

DO i=1,NNN
DO j=1,NNN

test=0.d0

      DO k=1,NN

       test=test+weights(k)*((abscc(i)*abscc(i)+abscc(j)*abscc(j))-2.d0*absc(k)*absc(k))&
              &*Vs(i,k)*Vs(k,j)

      END DO
       


function(i,j)=-((abscc(i)*abscc(i)-abscc(j)*abscc(j))**2.d0)*Vs(i,j)+test
END DO
END DO
return
END PROGRAM srg_flow
