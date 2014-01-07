! Functions and subroutines stolen from Numerical Recipes

module recipes

contains

  FUNCTION bessj0_s(x)
    USE nrtype; USE nrutil, ONLY : poly
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: x
    REAL(SP) :: bessj0_s
    REAL(SP) :: ax,xx,z
    REAL(DP) :: y
    REAL(DP), DIMENSION(5) :: p = (/1.0_dp,-0.1098628627e-2_dp,&
         0.2734510407e-4_dp,-0.2073370639e-5_dp,0.2093887211e-6_dp/)
    REAL(DP), DIMENSION(5) :: q = (/-0.1562499995e-1_dp,&
         0.1430488765e-3_dp,-0.6911147651e-5_dp,0.7621095161e-6_dp,&
         -0.934945152e-7_dp/)
    REAL(DP), DIMENSION(6) :: r = (/57568490574.0_dp,-13362590354.0_dp,&
         651619640.7_dp,-11214424.18_dp,77392.33017_dp,&
         -184.9052456_dp/)
    REAL(DP), DIMENSION(6) :: s = (/57568490411.0_dp,1029532985.0_dp,&
         9494680.718_dp,59272.64853_dp,267.8532712_dp,1.0_dp/)
    if (abs(x) < 8.0) then
       y=x**2
       bessj0_s=poly(y,r)/poly(y,s)
    else
       ax=abs(x)
       z=8.0_sp/ax
       y=z**2
       xx=ax-0.785398164_sp
       bessj0_s=sqrt(0.636619772_sp/ax)*(cos(xx)*&
            poly(y,p)-z*sin(xx)*poly(y,q))
    end if
  END FUNCTION bessj0_s
           
           
           
  FUNCTION bessj0_v(x)
    USE nrtype; USE nrutil, ONLY : poly
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: x
    REAL(SP), DIMENSION(size(x)) :: bessj0_v
    REAL(SP), DIMENSION(size(x)) :: ax,xx,z
    REAL(DP), DIMENSION(size(x)) :: y
    LOGICAL(LGT), DIMENSION(size(x)) :: mask
    REAL(DP), DIMENSION(5) :: p = (/1.0_dp,-0.1098628627e-2_dp,&
         0.2734510407e-4_dp,-0.2073370639e-5_dp,0.2093887211e-6_dp/)
    REAL(DP), DIMENSION(5) :: q = (/-0.1562499995e-1_dp,&
         0.1430488765e-3_dp,-0.6911147651e-5_dp,0.7621095161e-6_dp,&
         -0.934945152e-7_dp/)
    REAL(DP), DIMENSION(6) :: r = (/57568490574.0_dp,-13362590354.0_dp,&
         651619640.7_dp,-11214424.18_dp,77392.33017_dp,&
         -184.9052456_dp/)
    REAL(DP), DIMENSION(6) :: s = (/57568490411.0_dp,1029532985.0_dp,&
         9494680.718_dp,59272.64853_dp,267.8532712_dp,1.0_dp/)
    mask = (abs(x) < 8.0)
    where (mask)
       y=x**2
       bessj0_v=poly(y,r,mask)/poly(y,s,mask)
    elsewhere
       ax=abs(x)
       z=8.0_sp/ax
       y=z**2
       xx=ax-0.785398164_sp
       bessj0_v=sqrt(0.636619772_sp/ax)*(cos(xx)*&
            poly(y,p,.not. mask)-z*sin(xx)*poly(y,q,.not. mask))
    end where
  END FUNCTION bessj0_v



  SUBROUTINE polint(xa,ya,x,y,dy)
    USE nrtype; USE nrutil, ONLY : assert_eq,iminloc,nrerror
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: xa,ya
    REAL(SP), INTENT(IN) :: x
    REAL(SP), INTENT(OUT) :: y,dy
    INTEGER(I4B) :: m,n,ns
    REAL(SP), DIMENSION(size(xa)) :: c,d,den,ho
    n=assert_eq(size(xa),size(ya),'polint')
    c=ya
    d=ya
    ho=xa-x
    ns=iminloc(abs(x-xa))
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       if (any(den(1:n-m) == 0.0)) &
            call nrerror('polint: calculation failure')
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do
  END SUBROUTINE polint


end module recipes
