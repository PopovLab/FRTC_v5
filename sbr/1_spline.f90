module spline_module
    !! сплайны
    use kind_module
    implicit none
    
contains
    subroutine splne(x,y,n,y2)
        integer, parameter :: nn=3001
        real(wp), parameter :: zero=0d0
        integer n
        real(wp) x(n),y(n),y2(n),u(nn)
        integer i,k
        real(wp) p,qn,un,sig
        if(n.gt.nn) stop 'n>nn in splne!'
            y2(1)=zero
            u(1)=zero
            do i=2,n-1
                sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
                p=sig*y2(i-1)+2.d0
                y2(i)=(sig-1.d0)/p
                u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
            end do
            qn=zero
            un=zero
            y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
            do k=n-1,1,-1
                y2(k)=y2(k)*y2(k+1)+u(k)
            end do
    return
    end

    subroutine splnt(xa,ya,y2a,n,x,y,dy)
        real(wp), parameter :: zero=0d0
        integer n
        real(wp) xa(n),ya(n),y2a(n)
        integer k, klo, khi
        real(wp) x,h,a,b,aa,bb,hh,ax,bx,y,dy
        klo=1
        khi=n
        do while(khi-klo.gt.1)
            k=(khi+klo)/2
            if(xa(k).gt.x)then
                khi=k
            else
                klo=k
            endif
        end do
        h=xa(khi)-xa(klo)
        if(h.eq.zero) then
            write(*,*)'bad x input in splnt(), x=',x
            write(*,*)'klo=',klo,' kho=',khi
            stop
        end if
        a=(xa(khi)-x)/h
        b=(x-xa(klo))/h
        aa=a**2
        bb=b**2
        hh=h**2/6d0
        ax=-1d0/h
        bx=-ax
        y=a*ya(klo)+b*ya(khi)+(a*(aa-1d0)*y2a(klo)+b*(bb-1d0)*y2a(khi))*hh
        dy=ax*ya(klo)+bx*ya(khi)+ax*((3.d0*aa-1d0)*y2a(klo)-(3.d0*bb-1d0)*y2a(khi))*hh
    end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dsplnt(xa,ya,y2a,n,x,y,dy,ddy)
        real(wp), parameter :: zero=0d0
        integer n
        real(wp) xa(n),ya(n),y2a(n)
        integer k, klo, khi
        real(wp) x,h,a,b,aa,bb,hh,ax,bx,y,dy,ddy        
        klo=1
        khi=n
        do while(khi-klo.gt.1)
            k=(khi+klo)/2
            if(xa(k).gt.x)then
                khi=k
            else
                klo=k
            endif
        end do
        h=xa(khi)-xa(klo)
        if(h.eq.zero) then
            write(*,*)'bad x input in splnt(), x=',x
            write(*,*)'klo=',klo,' kho=',khi
            stop
        end if
        a=(xa(khi)-x)/h
        b=(x-xa(klo))/h
        aa=a**2
        bb=b**2
        hh=h**2/6d0
        ax=-1d0/h
        bx=-ax
        y=a*ya(klo)+b*ya(khi)+(a*(aa-1d0)*y2a(klo)+b*(bb-1d0)*y2a(khi))*hh
        dy=ax*ya(klo)+bx*ya(khi)+ax*((3.d0*aa-1d0)*y2a(klo)-(3.d0*bb-1d0)*y2a(khi))*hh
        ddy=6.d0*ax*ax*(a*y2a(klo)+b*y2a(khi))*hh
    end   
end module spline_module