module chebyshev
    !! Chebyshev fit
    use kind_module
    implicit none
    
contains
    SUBROUTINE chebft1(a,b,c,n,func)
    !! Chebyshev fit: Given a function func, lower and upper limits
    !! of the interval [a,b], and a maximum degree n, this routine 
    !! computes the n coefficients c(k) such that func(x) approximately =
    !! SUMM_(k=1)^(k=n)[c(k)*T(k-1)(y)]-c(1)/2, where y and x are related by
    !! (5.8.10). This routine is to be used with moderately large n 
    !! (e.g., 30 or 50), the array of cs subsequently to be truncated
    !! at the smaller value m such that c(m+1) and subsequent elements 
    !! are negligible. Parameters: Maximum expected value of n, and ð. 
        implicit none
        INTEGER n,NMAX
        real(wp) a,b,c(n),func,PI
        EXTERNAL func
        PARAMETER (NMAX=50, PI=3.141592653589793d0)
        INTEGER j,k
        real(wp) bma,bpa,fac,y,f(NMAX)
        real(wp) sum
        bma=0.5d0*(b-a)
        bpa=0.5d0*(b+a)
        do k=1,n
            y=cos(PI*(k-0.5d0)/n)
            f(k)=func(y*bma+bpa)
        end do
        fac=2.d0/n
        do j=1,n
            sum=0.d0
            do k=1,n
                sum=sum+f(k)*cos((PI*(j-1))*((k-0.5d0)/n))
            end do
            c(j)=fac*sum
        end do
        return
    END
    
    FUNCTION chebev(a,b,c,m,x)
    !! Chebyshev evaluation: All arguments are input. 
    !! c(1:m) is an array of Chebyshev coefficients, the first m elements 
    !! of c output from chebft (which must have been called with
    !! the same a and b). The Chebyshev polynomial evaluated
    !! and the result is returned as the function value.
        implicit none
        INTEGER m
        real(wp) chebev,a,b,x,c(m)
        INTEGER j
        real(wp) d,dd,sv,y,y2
        if ((x-a)*(x-b).gt.0.d0) pause 'x not in range in chebev'
        d=0.d0
        dd=0.d0
        y=(2.d0*x-a-b)/(b-a)
        y2=2.d0*y
        do j=m,2,-1
            sv=d
            d=y2*d-dd+c(j)
            dd=sv
        end do
        chebev=y*d-dd+0.5d0*c(1)
        return
    END

    SUBROUTINE chder(a,b,c,cder,n)
    !! Given a,b,c(1:n), as output from routine chebft(), and given n, 
    !! the desired degree of approximation (length of c to be used), 
    !! this routine returns the array cder(1:n), the Chebyshev 
    !! coefficients of the derivative of the function whose coefficients 
    !! are c(1:n).
        implicit none
        INTEGER n
        real(wp) a,b,c(n),cder(n)
        INTEGER j
        real(wp) con
        cder(n)=0.d0
        cder(n-1)=2*(n-1)*c(n)
        do j=n-2,1,-1
            cder(j)=cder(j+2)+2*j*c(j+1)
        end do
        con=2.d0/(b-a)
        do j=1,n
            cder(j)=cder(j)*con
        end do
    return
    END    
end module chebyshev