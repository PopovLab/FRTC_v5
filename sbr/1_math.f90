module math_module
    use kind_module
    implicit none
    
contains
    subroutine diff(x,y,n,dy)
        implicit real*8 (a-h,o-z)
        dimension y(*),x(*),dy(*)
        integer :: k, n
        dy(1)=(y(2)-y(1))/(x(2)-x(1))
        do k=2,n-1
            dy(k)=(y(k+1)-y(k-1))/(x(k+1)-x(k-1))
        end do
            dy(n)=(y(n)-y(n-1))/(x(n)-x(n-1))
        return
    end

    subroutine integral(ibeg,iend,x,y,fout)
        implicit real*8 (a-h,o-z)
        integer :: ibeg, iend
        dimension x(*),y(*)
        integer :: i, n1, n2, ie
        fout=0.d0
        if(ibeg.eq.iend) return
        znak=1.d0
        n1=ibeg
        n2=iend
        if(n2.lt.n1) then
            znak=-1.d0
            ie=n1
            n1=n2
            n2=ie
        end if
        sum=0.d0
        do i=n1+1,n2
            dx=x(i)-x(i-1)
            dsum=y(i)+y(i-1)
            sum=sum+.5d0*dsum*dx
        end do
        fout=znak*sum
    end    

    subroutine fsmoth4(x,y,n,ys)
        use constants, only : zero
        use approximation
        implicit real*8 (a-h,o-z)
        integer, intent(in) :: n
        !external polin2
        integer, parameter :: np=10, imax=601
        integer, parameter :: m0=1, ndp=1
        ! m0,ndp - parameters of smoothing procedure
        dimension y(n),x(n),ys(n)
        dimension yy(imax),xx(imax)
        dimension coeffs(np),cffs(np)
        dimension dys(imax)
        integer :: i, j, k, id
        integer :: m, m2, nmax, jlast
        if(n.gt.imax) stop 'small imax in subroutine fsmoth4()'
        call diff(x,y,n,dys)
        do k=1,n
            ys(k)=y(k)
        end do
        m=m0
        m2=m+2
        id=m+ndp
        nmax=n-id
        xs=x(1)
        do j=1,nmax
            do i=1,id
                xx(i)=x(j+i)-xs
                yy(i)=y(j+i)-ys(j)-dys(j)*xx(i)
            end do
            call approx(xx,yy,id,polin2,m,coeffs)
            cffs(1)=ys(j)
            cffs(2)=dys(j)
            do k=1,m
                cffs(k+2)=coeffs(k)
            end do
            xs=x(j+1)
            ys(j+1)=fdf(xx(1),cffs,m2,dys(j+1))
        end do
  
        j=nmax+1
  1     continue
        jlast=j
        id=n-jlast
        m=id-1
        if(m.eq.0) then
            j=j+1
            xs=x(j)
            ys(j)=fdf(xx(2),cffs,m2,dys(j))
            return
        end if
        m2=m+2
        do i=1,id
            xx(i)=x(jlast+i)-xs
            yy(i)=y(jlast+i)-ys(jlast)-dys(jlast)*xx(i)
        end do
        call approx(xx,yy,id,polin2,m,coeffs)
        cffs(1)=ys(jlast)
        cffs(2)=dys(jlast)
        do k=1,m
            cffs(k+2)=coeffs(k)
        end do
        j=j+1
        xs=x(j)
        ys(j)=fdf(xx(1),cffs,m2,dys(j))
        go to 1
    end    
end module math_module