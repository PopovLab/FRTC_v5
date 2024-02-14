module approximation
    !! polinomial approximation
    use kind_module    
    implicit none
    
contains

real(wp) function polin(k,x)
    implicit none
    integer k
    real(wp) x
    polin=1d0
    if(k.gt.1) polin=x**(k-1)
    return
end

real(wp) function polin1(k,x)
    implicit none
    integer k
    real(wp) x
    polin1=x**k
    return
end

real(wp) function polin2(k,x)
    implicit none
    integer k
    real(wp) x
    polin2=x**(k+1)
    return
end

subroutine approx(x,y,n,f,m,b)
!!     \(y(i)=y(x(i))\)  the data to be approximated.  
!!     \(n\)  number of points in the input data.  
!!     \(m\)  number of coefficients of decomposition
!!            over base functions \(f(k,x)\) :  
!!     \(y(x)=sum_1^m [b(k)*f(k,x)]\)  
!!     \(b(i)\)  found decomposition coefficients 

    implicit real*8 (a-h,o-z)
    integer,  parameter :: np=20
    real(wp), parameter :: zero=0.d0
    real(wp) a(np,np),indx(np)
    real(wp) y(n),x(n),b(*)
    integer i,j,k,m,n
    if(m.gt.np) then
        write(*,*)'index error subroutine "approx"'
        return
    end if

    do j=1,m
        do k=1,j
            a(k,j)=zero
            do i=1,n
                a(k,j)=a(k,j)+f(j,x(i))*f(k,x(i))
            end do
        end do
    end do
    do k=2,m
        do j=1,k-1
            a(k,j)=a(j,k)
        end do
    end do

    do k=1,m
     b(k)=zero
      do i=1,n
       b(k)=b(k)+y(i)*f(k,x(i))
      end do
    end do

    call ludcmp(a,m,np,indx,d)
    call lubksb(a,m,np,indx,b)
end

subroutine ludcmp(a,n,np,indx,d)
    implicit real*8 (a-h,o-z)
    integer,  parameter :: nmax=501
    real(wp), parameter :: tiny=1.d-20, zero=0.d0
    real(wp) a(np,np),indx(n),vv(nmax)
    integer i,j,k,m,n,np,imax
    d=1.d0
    do i=1,n
        aamax=zero
        do j=1,n
            if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
        end do
        if (aamax.eq.zero) pause 'singular matrix.'
            vv(i)=1.d0/aamax
    end do
    do j=1,n
    if (j.gt.1) then
    do i=1,j-1
      sum=a(i,j)
      if (i.gt.1)then
        do k=1,i-1
          sum=sum-a(i,k)*a(k,j)
        end do
        a(i,j)=sum
      endif
    end do
    endif
    aamax=zero
    do i=j,n
        sum=a(i,j)
        if (j.gt.1)then
            do k=1,j-1
                sum=sum-a(i,k)*a(k,j)
            end do
        a(i,j)=sum
        endif
        dum=vv(i)*dabs(sum)
        if (dum.ge.aamax) then
            imax=i
            aamax=dum
        endif
    end do
    if (j.ne.imax) then
        do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
        end do
        d=-d
        vv(imax)=vv(j)
    endif
    indx(j)=imax
    if (j.ne.n) then
        if (a(j,j).eq.zero) a(j,j)=tiny
        dum=1.d0/a(j,j)
        do i=j+1,n
            a(i,j)=a(i,j)*dum
        end do
    endif
    end do
    if(a(n,n).eq.zero) a(n,n)=tiny
    return
end

subroutine lubksb(a,n,np,indx,b)
    implicit real*8 (a-h,o-z)
    real(wp), parameter :: zero=0.d0
    real(wp)  a(np,np),indx(n),b(n)
    integer i,j,ii,ll,n,np 
    ii=0
    do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
            do j=ii,i-1
                sum=sum-a(i,j)*b(j)
            end do
        else if (sum.ne.zero) then
            ii=i
        endif
        b(i)=sum
    end do
    do i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
            do j=i+1,n
                sum=sum-a(i,j)*b(j)
            end do
        endif
        b(i)=sum/a(i,i)
    end do
    return
end    

    function fdf(x, c, n, df) result(p)
        !! вычисление значения полинома и его производной
        real(wp), intent(in)  :: x
        real(wp), intent(in)  :: c(n)
        integer,  intent(in)  :: n
        real(wp), intent(out) :: df
        integer               :: j
        real(wp)              :: p, dp
        p=c(n)
        dp=0.d0
        do j=n-1,1,-1
            dp=dp*x+p
            p=p*x+c(j)
        end do
        df=dp
    end

real(wp) function fdfddf(x,c,n,df,ddf)
    !! вычисление значения полинома и первой и второй производной
    real(wp), intent(in)   :: x
    real(wp), intent(in)   :: c(n)
    integer,  intent(in)   :: n 
    real(wp), intent(out)  :: df
    real(wp), intent(out)  :: ddf
    integer j
    real(wp) p, dp,ddp
    p=c(n)
    dp=0d0
    ddp=0d0
    do j=n-1,1,-1
        ddp=ddp*x+2d0*dp
        dp=dp*x+p
        p=p*x+c(j)
    end do
    fdfddf=p
    df=dp
    ddf=ddp
end


end module approximation



