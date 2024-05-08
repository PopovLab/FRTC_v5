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
        
    integer, parameter :: kpt1=20, kpt3=20

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

    subroutine gridvel(v1,v2,vmax,cdel,ni1,ni2,ipt1,kpt3,vrj)
        implicit none
        real(wp), intent(in) :: v1, v2, vmax, cdel
        integer,  intent(in) :: ni1, ni2, ipt1, kpt3
        real(wp), intent(inout) :: vrj(:)
        integer kpt1, kpt2, k
        real(wp) :: v12
        kpt1=ipt1-1
        kpt2=ni1+ni2+1
        do k=1,kpt1  !0<=v<v1
            vrj(k)=dble(k-1)*v1/dble(kpt1)
        end do
        v12=v1+(v2-v1)*cdel
        do k=1,ni1+1 !v1<=v<=v12
            vrj(k+kpt1)=v1+dble(k-1)*(v12-v1)/dble(ni1)
        end do
        do k=2,ni2+1 !!v12<v<=v2
            vrj(k+kpt1+ni1)=v12+dble(k-1)*(v2-v12)/dble(ni2)
        end do     
        do k=1,kpt3  !v2<v<=vmax
            vrj(k+kpt1+kpt2)=v2+dble(k)*(vmax-v2)/dble(kpt3)
        end do
    end    

    
    subroutine recalculate_f_for_a_new_mesh(ispectr)
        !!   recalculate f' for a new mesh
        use constants, only : zero
        use rt_parameters, only : nr, ni1, ni2
        use plasma, only: vt0, fvt, cltn
        use current, only: vzmin, vzmax
        use maxwell, only: i0, vij, dfij
        use lock_module        
        !use iterator_mod
        implicit none
        integer, intent(in) :: ispectr
        
        integer i, j, k
        real(wp) :: cdel, dfout
        real(wp), dimension(:), allocatable:: vvj, vdfj
        integer :: klo,khi,ierr
        real(wp) :: r, hr, vt, vto, vmax
        real(wp) :: v1, v2, vp1, vp2

        allocate(vvj(i0),vdfj(i0))
        hr = 1.d0/dble(nr+1)
        k=(3-ispectr)/2
        do j=1,nr
            r=hr*dble(j)
            vt=fvt(r)
            vto=vt/vt0
            if(iterat.gt.0) then
                v1=dmin1(vzmin(j),vz1(j))
                v2=dmax1(vzmax(j),vz2(j))
            else
                v1=vzmin(j)
                v2=vzmax(j)
            end if
            vmax=cltn/vto
            vp1=v1/vto
            vp2=v2/vto
            call gridvel(vp1,vp2,vmax,cdel,ni1,ni2,ipt1,kpt3,vrj)
            do i=1,i0
                vvj(i)=vij(i,j)
                vdfj(i)=dfij(i,j,k) !=dfundv(i,j)*vto**2
            end do
            do i=1,ipt
                call lock(vvj,i0,vrj(i),klo,khi,ierr)
                if(ierr.eq.1) then
            !!!         if(vrj(i).gt.vvj(i0)) exit
                    write(*,*)'lock error in new v-mesh'
                    write(*,*)'j=',j,' i0=',i0
                    write(*,*)'vvj(1)=',vvj(1),' vvj(i0)=',vvj(i0)
                    write(*,*)'i=',i,' vrj(i)=',vrj(i)
                    write(*,*)
                    pause'next key = stop'
                    stop
                end if
                call linf(vvj,vdfj,vrj(i),dfout,klo,khi)
                vgrid(i,j)=vrj(i)*vto
                dfundv(i,j)=dfout/vto**2
                if(dfundv(i,j).gt.zero) dfundv(i,j)=zero
            end do
            vz1(j)=v1
            vz2(j)=v2
        end do
        deallocate(vvj,vdfj)
    end subroutine    


end module iterator_mod