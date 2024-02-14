module chang_cooper_module
    use kind_module
    implicit none
    
contains
subroutine chang_cooper_solver(alfa2, nt, h, dt, n, ybeg, yend, d1,d2,d3, y)
    ! схема Ченга-Купера для уравнения Фоккера-Планка
    implicit none
    real(wp), intent(in)  :: alfa2      
    integer, intent(in) :: nt, n
    real(wp), intent(in)  :: h, dt
    real(wp), intent(in)  :: ybeg, yend
    real(wp), intent(in)  :: d1(n+1),d2(n+1),d3(n+1)
    real(wp), intent(inout) :: y(n+2)
    integer i, it, iz
    real(wp) xx(n+1), a(n),b(n),c(n),f(n)
    real(wp) y1(n)
    !print *, 'TK abc n = ', n
    do i=1,n+1
        xx(i)=h/2.d0+h*dble(i-1) !+shift
    end do
    
    do i=1,n
        y1(i)=y(i+1)
    end do
    
    do it=1, nt
           !call ABCcoef(a,b,c,f,y1,dt,n,ybeg,yend,x,xx,h,D1)
        call chang_cooper_abcoef(alfa2, a,b,c,f,y1, dt, n, ybeg, yend, xx, h, d1)
        call tridag(a,b,c,f,y1,n)
                  
        !do i=1,n
        !    if (y1(i).lt.0.d0) then
        !        if (y1(i) > epsilon(y1(i))) then
        !            y1(i)=0.d0
        !        else
        !            write(*,*) n, i, 'y(i)=',y1(i),' lt negative epsilon=',epsilon(y1(i))
        !            !pause
        !            !stop
        !        endif
        !    endif
        !enddo
        iz = n
        do i=1,n
            if (y1(i).lt.epsilon(yend)) then
                iz = i
                !print *, epsilon(yend)
                !print *, iz, n, yend, y1(i)
                exit
            endif
        enddo
        do i= iz, n
            y1(i)=yend
        enddo

    end do

    do i=1,n
        y(i+1)=y1(i)
    end do

end subroutine

! --
subroutine chang_cooper_abcoef(alfa2, A,B,C,f,Y,dt,n,ybeg,yend,xx,h,df)
    implicit none
    real(wp), intent(in)    :: alfa2
    real(wp), intent(inout) :: a(n),b(n),c(n),f(n),y(n+2)
    integer, intent(in)   :: n
    real(wp), intent(in)    :: dt, ybeg, yend, h
    real(wp), intent(in)    :: xx(n+1)
    real(wp), intent(in)    :: df(n+1)

    integer i
    real(wp) z, r, tmp1,tmp2,tmp3

    r=dt/h
    do i=1,n
        tmp1=dlt(xx(i),h, df(i), alfa2) * B1(xx(i),alfa2)
        A(i)=-r*(C1(xx(i),df(i))/h-tmp1)

        tmp2=C1(xx(i+1),df(i+1))/h-dlt(xx(i+1),h,df(i+1), alfa2) * B1(xx(i+1),alfa2)
        tmp3=(1.d0-dlt(xx(i),h,df(i), alfa2)) * B1(xx(i),alfa2)
        B(i)=r*(tmp2+tmp3+C1(xx(i),df(i))/h)+1.d0
  
        tmp1=(1.d0-dlt(xx(i+1),h,df(i+1), alfa2)) * B1(xx(i+1),alfa2)
        C(i)=-r*(tmp1+C1(xx(i+1),df(i+1))/h)

        f(i)=Y(i)
    enddo
    f(1)=f(1)-A(1)*ybeg
    f(n)=f(n)-C(n)*yend !yend in either way=0 all the time

    contains

    function B1(xx, alfa2) result(res)
        implicit none
        real(wp) xx,alfa2,beta, res
        res = -alfa2+1.d0/(xx*xx)
    end function

      
    function C1(xx,dif) result(res)
        implicit none
        real(wp) xx, dif, res
        res = dif+1.d0/(xx*xx*xx)
    end function
         
    function dlt(xx,h,dif, alfa2) result(res)
        implicit none
        real(wp) res
        real(wp) xx, h, dif, alfa2
        real(wp) w
        w = h*B1(xx, alfa2)/C1(xx,dif)
        res = 1.d0/w-1.d0/(dexp(w)-1.d0)
    end function

end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



end module Chang_Cooper_module