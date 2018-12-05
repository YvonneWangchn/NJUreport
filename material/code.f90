module integration
    implicit none
    real(kind=8) ::  a=0.0,b=1.0
contains
    real(kind=8) function f(x)
    implicit none
    real(kind=8) x
    f=sqrt(x)*log(x)
end 
end
program work6_1
    use integration
    implicit none
    real(kind=8):: T,S,h
    external T,S
    integer n
    open(10,file="..//Result//6_1.csv")
    do n=100,60000,100
        h=(b-a)/real(n);
        WRITE(*,"('h=',F12.9,4x,'T=',F12.9,4x,'S=',F12.9)") h,T(n),S(n)
        WRITE(10,"(F12.9,',',F12.9,',',F12.9)") h,T(n),S(n) 
    enddo
    call Romberg()
end program
!T,S,C,R stands for integration(using different methods)
!using trapzoid integration
real(kind=8) function T(n)
    use integration
    implicit none
    real(kind=8):: h
    integer k,n
    h=(b-a)/real(n);T=f(b)
    do k=1,n-1
    T=T+2.0*f(a+k*h)
    enddo
    T=h*T/2.0
end
!using simpson(h)
real(kind=8) function S(n)
    use integration
    implicit none
    real(kind=8):: h
    integer k,n
    h=(b-a)/real(n);S=f(b)
    S=S+4*f(a+h/2.0)
    do k=1,n-1
    S=S+4.0*f(a+k*h+h/2.0)+2.0*f(a+k*h)
    enddo
    S=h*S/6.0
end
!using Romberg
subroutine Romberg()
    use integration
    implicit none
    real(kind=8):: h,T(10)=0.0,TT,s,q
    real(kind=8),parameter :: e=1e-4
    integer j,k,m  
    T(1)=(b-a)/2.0*f(b);k=1;m=1;h=b-a
    do while(.true.)
        TT=0.0
        do j=0,k-1
        TT=TT+f(a+(j+0.5)*h)
        enddo
        TT=(TT*h+T(1))/2.0
        s=1.0
        do j=1,m
            s=4*s;q=(s*TT-T(j))/(s-1)!linear combination
            T(j)=TT;TT=q
        enddo
        write(*,"('m = ',I2,2x,'T(m+1)= ',F14.11,' T(m)=',F14.11,' difference = ',F14.11)") m,q,T(m),q-T(m)
        !if(abs(q-T(m))<=e)    exit 
        if(abs(q+4/9)<=e) exit
        m=m+1;T(m)=q;k=2*k
        h=h/2.0
    enddo
    write(*,"('m = ',I2,2x,'Romberg method = ',F14.11)") m,q
end

