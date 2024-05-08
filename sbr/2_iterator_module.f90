module iterator_mod
    use kind_module   
    implicit none
    real(wp) :: vmid(100),vz1(100),vz2(100)
    integer  :: ibeg(100),iend(100)

    real(wp) :: vrj(101),dj(101),djnew(1001)
    real(wp) :: dj2(101),d2j(101)

    real(wp) :: vgrid(101,100), dfundv(101,100)
    !!common/gridv/vgrid(101,100),dfundv(101,100)
    integer  :: nvpt
    !!common/gridv/nvpt
    integer :: ipt1, ipt2, ipt


    integer  :: iterat
    real(wp) :: psum4
    !!common /vvv2/ psum4
    real(wp) ::plost,pnab
    !!common /a0a4/ plost,pnab

    real(wp) :: vlf,vrt,dflf,dfrt
    !common /a0ghp/ vlf,vrt,dflf,dfrt
        
contains

    subroutine distr(vz,j,ifound,fder)
        !use iterator_mod
        use lock_module      
        implicit none
        integer, intent(in) :: j
        integer, intent(inout) :: ifound
        real*8 vz,fder
        integer i,klo,khi,ierr,nvp
        real*8,dimension(:),allocatable:: vzj,dfdvj
        real(wp) :: dfout
        !real*8 vlf,vrt,dflf,dfrt
        !common /a0ghp/ vlf,vrt,dflf,dfrt
        !common/gridv/vgrid(101,100),dfundv(101,100),nvpt

        nvp=nvpt
        allocate(vzj(nvp),dfdvj(nvp))
        do i=1, nvp
            vzj(i)=vgrid(i,j)
            dfdvj(i)=dfundv(i,j)
        end do
        call lock2(vzj,nvp,vz,klo,khi,ierr)
        if(ierr.eq.0) then !vgrid(1,j) <= vz <= vgrid(nvpt,j)
            call linf(vzj,dfdvj,vz,dfout,klo,khi)
            ifound=klo
            vlf=vzj(klo)
            vrt=vzj(khi)
            fder=dfout
            dflf=dfdvj(klo)
            dfrt=dfdvj(khi)
        else if(ierr.eq.1) then !vz < vgrid(1,j)
            write(*,*)'exception: ierr=1 in distr()'
            pause'next key = stop'
            stop
        else if(ierr.eq.2) then !vz > vgrid(nvpt,j)
            write(*,*)'exception: ierr=2 in distr()'
            pause'next key = stop'
            stop
        else if(ierr.eq.3) then
            write(*,*)'exception in distr, klo=khi=',klo,' j=',j,' nvp=',nvp
            write(*,*)'vz=',vz,' v1=',vzj(1),' v2=',vzj(nvp)
            pause'next key = stop'
            stop
        end if
        deallocate(vzj,dfdvj)
    end   


    subroutine calculate_dfundv(ispectr)
        !! calculate dfundv что такое dfundv?
        use constants, only: zero
        use rt_parameters, only: nr
        use plasma, only: fvt, vt0, cltn
        use maxwell, only: i0, vij, dfij, dij
        !use iterator_mod, only: dfundv
        !use iterator_mod, only: ipt
        !use iterator_mod, only: vrj, dj, vgrid
        use lock_module, only: lock, linf
        implicit none
        integer, intent(in) :: ispectr
        real(wp), dimension(:), allocatable:: vvj, vdfj
        integer  :: i, j, k 
        integer  :: klo,khi,ierr
        real(wp) :: dfout
        real(wp) :: r, hr
        real(wp) :: vt, vto
        hr = 1.d0/dble(nr+1)
        allocate(vvj(i0),vdfj(i0))
        k=(3-ispectr)/2
        do j=1,nr
            r=hr*dble(j)
            vt=fvt(r)
            vto=vt/vt0
            do i=1,i0
                vvj(i)=vij(i,j)
                vdfj(i)=dfij(i,j,k) !=dfundv(i,j)*vto**2
            end do
            do i=1,ipt
                vrj(i)=vgrid(i,j)/vto   !Vpar/Vt
                call lock(vvj,i0,vrj(i),klo,khi,ierr)
                if(ierr.eq.1) then
                    write(*,*)'lock error in read distribution function'
                    write(*,*)'j=',j,'i0=',i0
                    write(*,*)'vvj(1)=',vvj(1),' vvj(i0)=',vvj(i0)
                    write(*,*)'i=',i,' vrj(i)=',vrj(i),' vmax=',cltn/vto
                    write(*,*)
                    pause'next key = stop'
                    stop
                end if
                call linf(vvj,vdfj,vrj(i),dfout,klo,khi)
                dfundv(i,j)=dfout/vto**2
                if(dfundv(i,j).gt.zero) dfundv(i,j)=zero
            end do
        end do
        deallocate(vvj,vdfj)
    end


end module iterator_mod