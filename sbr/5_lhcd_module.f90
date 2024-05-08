module lhcd_module
    !! LHCD модуль
    use kind_module
    use spectrum_mod
    implicit none
    
    type(Spectrum) :: full_spectrum
    type(Spectrum) :: pos_spectr, neg_spectr

    !real(wp), dimension(:), allocatable:: vvj, vdfj

    integer, parameter :: kpt1=20, kpt3=20

contains
    subroutine ourlhcd2017(spectr, outpe, pe_out)      
        use constants
        use approximation
        use spline_module
        use chebyshev
        use plasma
        use rt_parameters
        use maxwell      
        use trajectory_module, only: view,  init_trajectory
        use spectrum_mod
        use manager_mod
        use dispersion_module
        use current
        use iteration_result_mod
        use iterator_mod, only: pnab, plost, psum4
        use iterator_mod, only: ipt, ipt1, ipt2, nvpt, iterat
        use iterator_mod, only: vrj, vz1, vz2, vgrid
        use iterator_mod, only: calculate_dfundv, gridvel
        use lock_module
        use math_module
        !use driver_module, only : lfree
        !use driven_current_module, only : zv1, zv2
        use decrements, only: kzero
        use source_new_mod
        
        implicit none
        type(Spectrum) spectr
        real*8 outpe,pe_out 
        dimension outpe(*)

        real(wp) :: hr, r,  pn, fnr, fnrr
        real(wp) :: vt, vto, wpq, whe, v
        real(wp) :: u, u1, e1, e2, e3, tmp
        real(wp) :: cn1, avedens
        real(wp) :: anb, fuspow, o_da
        real(wp) :: dvperp, ddens
        real(wp) :: tt, cn2, vmax, v1, v2
        real(wp) :: tdens

        real(wp) :: q_rest, q_abs, q_cond
        real(wp) :: pchg

        real(wp) :: oi
        real(wp) :: ol, oc, oa, of
        real(wp) :: zff, cnyfoc, dconst, fout
        real(wp) :: galfa(50,100), vpmin(100), vcva(100)
        real(wp) :: pdprev1(100), pdprev2(100)
        real(wp) :: source(100)

        real(wp) ::  rmx_n, rmx_t, rmx_z, rmx_ti
        ! встречает только один раз common /maxrho/ rmx_n,rmx_t,rmx_z,rmx_ti
      
        type(IterationResult) :: iteration_result
        real*8 kofpar,timecof

        integer     :: iptnew
        real(wp)    :: plaun

        integer ispectr
        integer :: i, j, k  
        integer :: klo,khi,ierr
        integer :: iww, iw0, izz

        plaun = spectr%input_power

        ispectr = spectr%direction
        !lfree=1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        hr = 1.d0/dble(nr+1)
        iw0=iw
    
        call find_volums_and_surfaces

        ppv1=zero
        ppv2=zero
        pnab=zero
        plost=zero
        psum4=zero
        anb=zero
        fuspow=zero
        o_da=zero
        !c-------------------------------------------
        !c find velocity limits and initial dfdv
        !c--------------------------------------------
        ipt1=kpt1+1
        ipt2=ni1+ni2
        ipt=ipt1+ni1+ni2+kpt3
        if(ipt.gt.101) then
            write(*,*)'ipt >101'
            pause'stop program'
            stop
        end if
        nvpt=ipt

        do j=1,nr                  ! begin 'rho' cycle
            r=hr*dble(j)
            !!!!sav2008       pn=fn(r)
            !!       pn=fn1(r,fnr)
            !!       pn=fn2(r,fnr,fnrr) !sav2008
            if(inew.eq.0) then !vardens
                pn=fn1(r,fnr)
            else
                pn=fn2(r,fnr,fnrr)
            end if
            dens(j)=pn
            vt=fvt(r)
            vto=vt/vt0
            wpq=c0**2*pn
            whe=dabs(b_tor)*c1
            v=wpq/ww**2
            u1=whe/ww
            u=u1**2
            e1=1d0-v*(1d0/xmi-1d0/u)
            e2=v/u1
            e3=v
            tmp=ft(r)/0.16d-8 !Te, keV
            cn1=dsqrt(50d0/tmp)  !sav2008
            if(itend0.gt.0) then
                eta(j)=1d0-v
                vcva(j)=cnstvc*vt*dsqrt(2d0)/valfa
                vpmin(j)=2.0d0*dsqrt(tmp/(-eta(j)))
222             continue
                dvperp=(vpmax-vpmin(j))/dble(kv-1)
                if(dvperp.le.zero) then
                    vpmax=1.3d0*vpmax
                    go to 222
                end if
                do k=1,kv
                    vperp(k,j)=vpmin(j)+dble(k-1)*dvperp
                end do
                fcoll(j)=.5d-13*dens(j)*zalfa**2*xlog/xmalfa/tmp**1.5d0
                ddens=dn1*dens(j)
                tdens=dn2*dens(j)
                tt=fti(r)**one_third    ! (ti, keV)^1/3
                source(j)=4d-12*factor*ddens*tdens*dexp(-20d0/tt)/tt**2
                anb=anb+source(j)*vk(j)
            end if
            cn2=dsqrt(dabs(e1))+e2/dsqrt(e3) !sav2008
            !vz1(j)=cleft*cltn/cn1  !Vpar/Vt0
            !vz2(j)=cright*cltn/cn2  !Vpar/Vt0
            !if(vz2(j).gt.0.9d0*cltn) vz2(j)=0.9d0*cltn
            !v1=vz1(j)/vto !Vpar/Vt(rho)
            !v2=vz2(j)/vto !Vpar/Vt(rho)
            vmax=cltn/vto
            v1=4.d0  !Vpar/Vt(rho)
            v2=10.d0 !cright*cltn/cn2 !10.d0 !Vpar/Vt(rho)
            if(v2.ge.vmax) v2=0.5d0*vmax
            if(v1.ge.v2) v1=v2-2.d0
            call gridvel(v1,v2,vmax,0.5d0,ni1,ni2,ipt1,kpt3,vrj)
            vz1(j)=v1*vto !Vpar/Vt0
            vz2(j)=v2*vto !Vpar/Vt0
            if(vz2(j).gt.0.9d0*cltn) vz2(j)=0.9d0*cltn
            do i=1,ipt
                vgrid(i,j)=vrj(i)*vto
            end do
        end do                     ! end 'rho' cycle 


        !!!!!!!!!read data !!!!!!!!!!!!       
        call calculate_dfundv(ispectr)
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(itend0.gt.0) then  ! begin alpha-source renormalisation
            call alpha_source_renormalisation(anb, fuspow, source)
        end if
        !c------------------------------------
        !c set initial values of arrays
        !c------------------------------------
        iww=0
        izz=zero
        ! 
        pdl=zero
        pdc=zero
        pda=zero
        pdfast=zero

        !! массивы для невзязки
        pdprev1=zero
        pdprev2=zero

        dql=zero
        dq1=zero
        dq2=zero
        dncount=zero
        vzmin=cltn
        vzmax=-cltn
        kzero=kv
        call init_trajectory

        if(itend0.gt.0) then
            do j=1,nr           ! begin 'rho' cycle
                do i=1,50
                    dqi0(i,j)=zero
                end do
                call alphas(dqi0,vperp,j,kv,galfa)
            end do              ! end 'rho' cycle
        end if
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!sign of driven current in right coordinate system {dro,dteta,dfi}:
        !!!!!curdir=+1.0 for current drive in positive direction "dfi"
        !!!!!curdir=-1.0 for current drive in negative direction "dfi"
        !!!!!spectrum Nz>0 is along dfi>0 and Nz<0 is along dfi<0
        !!!!!it is also OK if Npar is used instead of Nz, but for Btor>0, that is along dfi>0
        !!      curdir=-dble(ispectr)
        !!!!!!!!!!!!!!! begin iterations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        q_rest=plaun
        iterat=0

    80  continue
        call manager(iterat, iw0, ntet, spectr)

        call find_achieved_radial_points(nvpt)

        call renormalisation_power
        
        pchg = find_nevyazka(pdprev1, pdprev2)

        call calculate_total_current_and_power(ol, oc, oa, of)

        !!!!!!!!! prepare to the next iteration !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        iterat=iterat+1
        q_abs=ol+oc+oa
        q_rest=pnab/xsgs
        q_cond=zero
        if(q_abs.ne.zero) q_cond=0.5d0*q_rest/q_abs
        !!!      if(q_cond.le.pabs0.and.pchg.lt.pgiter)

        call integral(1,nspl,rh,con,avedens) 

        iteration_result = IterationResult(number = iterat, &
            spectr_direction = ispectr, P_launched = plaun, &
            P_landau = ol, P_coll = oc, P_alph = oa, &
            alphas_power = fuspow, P_fast = of, &
            P_lost = plost/xsgs, P_not_accounted = pnab/xsgs, &
            P_landau_strong_absorption = ppv1/xsgs, &
            P_landau_weak_absorption = ppv2/xsgs, &
            P_turns = psum4/xsgs, efficiency = oi/plaun, &
            avedens = avedens*1.d19, r0 = r0*1.d-2, &
            eta_eff = 1.d17*avedens*r0*oi/plaun, &
            residual = pchg)

        if(iterat.gt.5.and.q_cond.le.pabs0.and.pchg.lt.pgiter) goto 110

        if(ipri.gt.1) then
            call iteration_result%print
            call iteration_result%save(tcur)
            !pause
        end if

        if(iterat.le.niterat) then
            call recalculate_f_for_a_new_mesh(ispectr)
            call init_iteration
            goto 80
        end if
        !c------------------------------------------
        !c save results
        !c------------------------------------------
110     continue

        if(ipri.gt.0) then
            write (*,*)
            write (*,*) 'RAY-TRACING RESULTS:'
            call iteration_result%print
            call iteration_result%save(tcur)
            write (*,*) '-------------------------------------------'
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call calculate_diffusion(ispectr)

        call view(tcur,ispectr,spectr%size,ntet)  !writing trajectories into a file
        call calculate_out_power(outpe)
        pe_out=ol+oc
        
    end    
    
    subroutine alpha_source_renormalisation(anb, fuspow, source)
        use constants, only: zero, talfa, one_third
        use rt_parameters, only: nr, dra, factor
        use source_new_mod, only: rsou, sou, npta
        use plasma, only: fti, dn1, dn2, vk
        use current, only: dens
        implicit none
        real(wp), intent(inout) :: anb
        real(wp), intent(inout) :: fuspow
        real(wp), intent(inout) :: source(:)
        real(wp) :: r, hr, tt
        real(wp) :: anb0, aratio, sssour
        real(wp) :: ddens, tdens
        real(wp) :: sour(100)
        integer j
        hr = 1.d0/dble(nr+1)
        fuspow=anb*talfa*1.6022d-19
        anb0=anb
        anb=zero
        do j=1,nr
            r=hr*dble(j)
            if(r.le.dra) then
                tt=fti(zero)**one_third
            else
                tt=fti(r-dra)**one_third    ! (shifted ti, kev)^1/3
            end if
            ddens=dn1*dens(j)
            tdens=dn2*dens(j)
            sour(j)=4d-12*factor*ddens*tdens*dexp(-20d0/tt)/tt**2
            anb=anb+sour(j)*vk(j)
        end do
        aratio=anb0/anb
        rsou(1)=zero
        sou(1)=aratio*sour(1)
        do j=1,nr
            r=hr*dble(j)
            rsou(j+1)=r
            sou(j+1)=aratio*sour(j)
            if(j.eq.nr) sssour=source(j)
            source(j)=sou(j+1)
        end do
        npta=nr+2
        rsou(npta)=1.d0
        sou(npta)=aratio*sour(nr)        
    end

    subroutine calculate_diffusion(ispectr)
        !! calculate diffusion
        use constants, only: c0, pme, zero
        use rt_parameters, only: nr, inew, ni1, ni2
        use plasma, only: zefff, fn1, fn2
        use plasma, only: vt0, fvt, cltn, cnye
        use driven_current_module, only : zv1, zv2
        use current, only : dql
        use maxwell, only: i0, vij, dfij, dij
        use iterator_mod, only: ipt, ipt1
        use iterator_mod, only: vrj, dj, vgrid
        use lock_module, only: lock
        implicit none
        integer, intent(in) :: ispectr
        integer  :: i,j,k
        integer  :: klo,khi,ierr
        real(wp) :: r, hr
        real(wp) :: vt, vto, vmax, zff, cnyfoc
        real(wp) :: pn, fnr, fnrr
        real(wp) :: dconst, ddout
        real(wp) :: dijk(101,100,2), vrjnew(101,100,2)
        !встречает только один раз common/t01/dijk(101,100,2), vrjnew(101,100,2), iptnew
        !
        hr = 1.d0/dble(nr+1)
        k=(3-ispectr)/2
        do j=1,nr
            r=hr*dble(j)
            vt=fvt(r)
            vto=vt/vt0
            vmax=cltn/vto
            zff=(5d0+zefff(r))/5d0
            cnyfoc=zff*c0**4*cnye
            if(inew.eq.0) then !vardens
                pn=fn1(r,fnr)
            else
                pn=fn2(r,fnr,fnrr)
            end if
            dconst=vt0/(1.d-10*cnyfoc*pme*pn**2) !divided by 10^-10 here 
            !!!!!!!!                         and multiplied by 10^-10 in dfind()
            !!!old       dconst=vt0/(cnyfoc*pme*pn**2)
            !!!        dj(i)=dql(i,j)*dconst*vto !D_normir
            do i=1,ipt
                vrj(i)=vgrid(i,j)/vto      !Vpar/Vt
                dj(i)=dql(i,j)*dconst*vto  !D_normir
                vrjnew(i,j,k)=vrj(i)
                dijk(i,j,k)=dj(i)
            end do
            do i=1,i0
                if(vij(i,j).ge.vmax) then
                    ddout=zero
                else
                    call lock(vrj,ipt,vij(i,j),klo,khi,ierr)
                    if(ierr.eq.1) then
                        write(*,*)'lock error in output dql'
                        write(*,*)'j=',j,'ipt=',ipt
                        write(*,*)'vrj(1)=',vrj(1),' vrj(ipt)=',vrj(ipt)
                        write(*,*)'i=',i,' v=',vij(i,j),' vmax=',vmax
                        write(*,*)
                        pause'next key = stop'
                        stop
                    end if
                    !!!         call linf(vrj,dj,vij(i,j),ddout,klo,khi)
                    !!         if(ddout.le.1.d0) ddout=zero
                    ddout=dj(klo)
                end if
                dij(i,j,k)=ddout
            end do
            zv1(j,k)=vrj(ipt1)
            zv2(j,k)=vrj(ni1+ni2+ipt1)
        end do
    end 

    subroutine calculate_out_power(out_pe)
        use constants, only: zero, one
        use rt_parameters, only: nr, ismthout
        use plasma, only: rh, rh1, nspl, vk
        use current, only: pdl, pdc
        use math_module, only: fsmoth4
        use lock_module, only: lock2, linf
        implicit none
        real(wp), intent(inout) :: out_pe(*)
        real(wp) :: rxx(102)
        integer  :: i, j, nrr
        integer  :: klo,khi,ierr
        real(wp) :: hr, fout
        real(wp) :: pwe(102),wrk(102)

        hr = 1.d0/dble(nr+1)
        nrr=nr+2
        rxx(1)=zero
        rxx(nrr)=one
        do j=1,nr
            rxx(j+1)=hr*dble(j)
        end do

        if(ismthout.ne.0) then
            do j=1,nr
                pwe(j+1)=(pdl(j)+pdc(j))/vk(j)
            end do
            pwe(1)=pwe(2)
            pwe(nr+2)=zero
            do i=1,nrr
                wrk(i)=pwe(i)
            end do
            call fsmoth4(rxx,wrk,nrr,pwe)
        end if
        !
        rh(1)=rh1
        if(rh(nspl).gt.1.d0) rh(nspl)=1.d0
        do j=1,nspl
            call lock2(rxx,nrr,rh(j),klo,khi,ierr)
            if(ierr.ne.0) then
                write(*,*)'lock2 error in profiles for ASTRA'
                write(*,*)'ierr=',ierr,' j=',j,' rh(j)=',rh(j)
                write(*,*)'rxx(1)=',rxx(1),' rxx(nrr)=',rxx(nrr)
                pause
            end if
            call linf(rxx,pwe,rh(j),fout,klo,khi)
            out_pe(j)=fout
        end do
        rh(1)=zero
        !
    end 

    subroutine init_iteration
        use constants, only : zero
        use rt_parameters, only : itend0
        use current
        use iterator_mod
        use plasma, only: cltn
        implicit none
        ppv1=zero
        ppv2=zero
        psum4=zero
        pnab=zero
        plost=zero
        dql=zero
        dq1=zero
        dq2=zero
        dncount=zero
        vzmin=cltn
        vzmax=-cltn
        pdl=zero
        pdc=zero
        pda=zero
        pdfast=zero
        if(itend0.gt.0) then
              dqi0=zero
        end if
    end 



    subroutine recalculate_f_for_a_new_mesh(ispectr)
        !!   recalculate f' for a new mesh
        use constants, only : zero
        use rt_parameters, only : nr, ni1, ni2
        use plasma, only: vt0, fvt, cltn
        use current, only: vzmin, vzmax
        use maxwell, only: i0, vij, dfij
        use lock_module        
        use iterator_mod
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

    subroutine alphas(d,u,j,kmax,g)
        use decrements, only: dgdu, kzero
        use constants, only : zero, one
        implicit real*8 (a-h,o-z)
        integer, intent(in) :: j, kmax
        dimension d(50,100),u(50,100),g(50,100)
        !common /arr/ dgdu(50,100),kzero(100)
        real(wp), parameter :: tiny=1.d-30
        integer :: k, km
        km=kzero(j)
        um=u(km,j)
        if(um.ge.one) then
            do k=1,kmax
                if(u(k,j).lt.one) then
                    uk=u(k,j)
                    uk2=uk**2
                    w=dsqrt(one-uk2)
                    g(k,j)=w/uk2
                    dgdu(k,j)=-one/(w*uk)-2.d0*w/(uk*uk2)
                else
                    g(k,j)=zero
                    dgdu(k,j)=zero
                end if
            end do
            return
        end if
  
        do k=1,km
            uk=u(k,j)
            uk2=uk**2
            w=dsqrt(one-uk2)
            g(k,j)=w/uk2
            dgdu(k,j)=-one/(w*uk)-2.d0*w/(uk*uk2)
        end do
  
        do k=km+1,kmax
            du=u(k,j)-u(k-1,j)
            if(u(k,j).lt.one) then
                beta=u(k,j)*dsqrt(one-u(k,j)**2)
            else
                beta=zero
            end if
            alfa=u(k,j)**3
            g(k,j)=(d(k,j)*g(k-1,j)+beta*du)/(d(k,j)+alfa*du)
            if(d(k,j).ne.zero) then
                dgdu(k,j)=(beta-alfa*g(k,j))/d(k,j)
            else
                dgdu(k,j)=(g(k,j)-g(k-1,j))/du
            end if
        end do
        do k=1,kmax
            if(g(k,j).lt.tiny) g(k,j)=zero
            if(dabs(dgdu(k,j)).lt.tiny) dgdu(k,j)=zero
        end do
        return
    end    
end module lhcd_module