module lock_module
    use kind_module
    implicit none
    
contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine linf(x,y,t,fout,klo,khi)
        !! линейная аппроксимация
        !! TODO будет удобнее если переделать в функцию
        implicit none
        real(wp), intent(in)  :: x(*),y(*), t
        real(wp), intent(out) :: fout
        integer,  intent(in)  :: klo, khi
        real(wp) :: dout
        dout=(y(khi)-y(klo))/(x(khi)-x(klo))
        fout=y(klo)+dout*(t-x(klo))
    end

    subroutine lock(xa,n,x,klo,khi,ierr)
        !! что делает?
        implicit none
        real(wp),   intent(in) :: xa(*), x
        integer,    intent(in) :: n
        integer, intent(inout) :: klo, khi, ierr
        real(wp), parameter    :: tiny=1.d-14
        integer  :: k
        real(wp) :: dx1, dx2
        klo=0
        khi=0
        dx1=x-xa(1)
        dx2=x-xa(n)
        if(dx1*dx2.ge.tiny) then
            ierr=1
            return
        end if
        ierr=0
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
        if(khi.eq.klo) ierr=1
    end

    subroutine lock2(xa,n,x,klo,khi,ierr)
        !! что делает
        use constants, only: zero
        implicit none
        real(wp),   intent(in) :: xa(*), x
        integer,    intent(in) :: n
        integer, intent(inout) :: klo, khi, ierr
        real(wp), parameter    :: tiny=1.d-7
        integer  :: k
        real(wp) :: dx1, dx2

        ierr=0
        klo=0
        khi=0
        dx1=x-xa(1)
        if (abs(dx1).lt.tiny) then
            klo=1
            khi=2
            return
        else if(dx1.lt.zero) then
            ierr=1
            return
        end if
        dx2=x-xa(n)
        if (abs(dx2).lt.tiny) then
            klo=n-1
            khi=n
            return
        else if(dx2.gt.zero) then
            ierr=2
            return
        end if
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
        if(khi.eq.klo) ierr=3
    end
end module lock_module