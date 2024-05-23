module trajectory_module
    use kind_module
    use trajectory_data
    implicit none

    !integer, parameter :: mpnt = 100000


    !integer nrefj(mpnt)
    !!common/refl/nrefj(mpnt)

    !integer mbeg(mpnt),mend(mpnt),mbad(mpnt)

    integer, parameter :: max_num_trajectories = 30000
    type(Trajectory), target ::  trajectories(max_num_trajectories)
contains

subroutine init_trajectory
    use constants
    use driver_module
    implicit none
    !nrefj = 0
    
    !dland = zero
    !dcoll = zero
    !perpn = zero 
    !dalf  = zero
    !vel = zero
    !jrad = zero
    !iww = zero
    !tetai = zero
    !xnpar = zero
    !izz = zero

    !mbeg = zero
    !mend = zero
    !mbad = zero

end subroutine 

subroutine view(tview, ispectr,nnz,ntet) !sav2008
!!!writing trajectories into a file
    use constants
    use approximation
    use plasma
    use decrements, only: pdec1,pdec2,pdec3,pdecv,pdecal,dfdv
    use decrements, only: zatukh
    use rt_parameters, only :  nr, itend0, kv, nmaxm    
    use iterator_mod, only : dflf, dfrt, distr
    use driver_module !, only: jrad, iww, izz, length
    use trajectory_data
    implicit none
    
    real(wp), intent(in) :: tview

    integer, intent(in) :: ispectr, nnz, ntet  !sav#

    type(Trajectory) :: traj
    type(TrajectoryPoint) ::tp
    !common /bcef/ ynz,ynpopq
    !common /vth/ vthc(length),poloidn(length)
    real(wp) vthcg,npoli
    !common /a0ghp/ vlf,vrt,dflf,dfrt
    
    integer i, n, itr, ntraj
    integer jrc,nturn,ib,ie,jr,ifast,idir,iv
    integer jznak,jdlt,mn,mm,jchek,itet,inz
    integer, parameter :: unit_bias = 10
    integer, parameter :: m=7
    real(wp), parameter :: pleft=1.d-10 !m may be chaged together with name(m)

    real(wp) :: htet, h, xr, xdl, xdlp, xly, xlyp, xgm, xgmp, th
    real(wp) :: x, xx, z, zz, pl, pc, pa
    real(wp) :: pdec1z, pdec3z, pintld, pintal
    real(wp) :: cotet, sitet
    real(wp) :: v, refr, dek3, parn, argum
    real(wp) :: df, powpr, powd, powal, pil, pic
    real(wp) :: powcol, pia, pt, denom, powdamped, domin, fff
    real(wp) :: pow 
    character(14) folder
    character(40) ver_fn
    character(40) fname

    print *, 'view_time=',tview
    !print *, name(m)
    if (ispectr>0) then
        folder = "lhcd/traj/pos/"
    else
        folder = "lhcd/traj/neg/"
    endif

    write(ver_fn,'(A, "v2")') folder
    print *, ver_fn
    open(1,file=ver_fn)
    write (1,*) 'version 2'
    close(1)

    write(fname,'(A, f9.7,".dat")') folder, tview
    print *, fname



    htet=zero
    h=1d0/dble(nr+1)
    if(ntet.ne.1) htet=(tet2-tet1)/(ntet-1)


    open(1,file=fname)

    ntraj=0 !sav2008
    do itr=1,nnz*ntet !sav2008
        pow=1.d0
        pl=zero
        pc=zero
        pa=zero
        pdec1=zero
        pdec1z=zero
        pdec3=zero
        pdec3z=zero
        pdecv=zero
        pintld=zero
        pintal=zero
        jrc=nr+1
        jznak=-1
        nturn=1

        traj = trajectories(itr)
        call traj%write_info(1)

        if(traj%mbad.eq.0) then 
            ntraj=ntraj+1
            write(1,3) !write header 
            do i=1, traj%size
                tp = traj%points(i)
                v  = tp%vel
                jr = tp%jrad
                refr  = tp%perpn
                npoli = tp%poloidn
                ifast = tp%iww
                vthcg = tp%vthc
                idir  = tp%izz
                dek3  = zero
                th = tp%tetai
                parn = tp%xnpar
                if(itend0.gt.0) then
                    argum=clt/(refr*valfa)
                    dek3=zatukh(argum,abs(jr),vperp,kv)
                end if

                call distr(v,abs(jr),iv,df)
                if(jr.lt.0) then    !case of turn
                    jr=-jr
                    !variant          pintld=-dland(i)*df
                    !!          pintld=-dland(i)*(dflf+dfrt)/2d0
                    pintld = dabs(tp%dland*(dflf+dfrt)/2d0)
                    pdec2  = dexp(-2d0*tp%dcoll)
                    pintal = dabs(tp%dalf*dek3)
                else
                    pdec2 = tp%dcoll
                    pdecv = tp%dland
                    !!          pdec1=-pdecv*df
                    pdec1 = dabs(pdecv*df)
                    pdec3 = dabs(tp%dalf*dek3)
                    pintld = (pdec1+pdec1z)/2d0*h
                    pintal = (pdec3+pdec3z)/2d0*h
                    pdec1z = pdec1
                    pdec3z = pdec3
                end if
                powpr=pow
                powd=pow*dexp(-2d0*pintld)
                powcol=powd*pdec2
                powal=powcol*dexp(-2d0*pintal)
                pow=powal
                pil=pintld
                pic=.5d0*dabs(dlog(pdec2))
                pia=pintal
                pt=1.d0-pow  !total absorbed power
                denom=pil+pic+pia
                powdamped=1.d0-dexp(-2.d0*denom)
                domin=powpr*powdamped
                if(denom.ne.zero) then
                    fff=domin/denom
                    pl=pl+dabs(pil*fff)  !el. Landau absorbed power
                    pc=pc+dabs(pic*fff)  !el. collisions absorbed power
                    pa=pa+dabs(pia*fff)  !alpha Landau absorbed power
                end if
                xr= tp%rho !h*dble(jr)
                cotet=dcos(th)
                sitet=dsin(th)
                xdl=fdf(xr,cdl,ncoef,xdlp)
                xly=fdf(xr,cly,ncoef,xlyp)
                xgm=fdf(xr,cgm,ncoef,xgmp)
                xx=-xdl+xr*cotet-xgm*sitet**2
                zz=xr*xly*sitet
                x=(r0+rm*xx)/1d2
                z=(z0+rm*zz)/1d2 
                jdlt=jr-jrc
                jrc=jr
                if(jdlt*jznak.lt.0.and.nturn.lt.m-1) then
                    nturn=nturn+1
                    jznak=-jznak
                end if

                        !   R, Z, rho, theta, N_par, N_pol, P_tot,P_land, P_coll, vth,  'slow=1','idir=1','ipow', driver, N_traj 
                write(1, 7) x, z, xr,  th,    parn,  npoli, pt,   pl,     pc,     vthcg, ifast,   idir, tp%ipow,  tp%driver,  itr
                
                if(pt.ge.1d0-pleft) go to 11 !maximal absorbed power along a ray
            end do
11          continue
            write (1,*)
        end if
    end do


    close(1)


1     format(2x,'N_traj',3x,'mbad',6x,'theta',9x,'Npar',9x,'rho_start')
2     format('R_pass',4x,'Ptot',6x,'Pland',6x,'Pcoll',8x,'Pa',7x,'dPtot',6x,'dPland',5x,'dPcoll',6x,'dPa')
3     format(5x,'R',10x,'Z',11x,'rho',8x,'theta',7x,'N_par',7x,'N_pol',6x,'P_tot',7x,'P_land',6x,'P_coll',6x,'vth',4x,'slow=1',4x,'out=1',4x,'ipow',4x,'driver',4x'N_traj',6x)
4     format(i3,5x,8(f6.3,5x))
5     format(6(e13.6,3x))
6     format(2(i6,2x),4(e13.6,1x))
7     format(10(e11.4,1x),i5,2x,i5,2x,i3,2x,i5,2x,i5)
8     format('after radial pass=',i3,2x,' P_tot=',f6.3,2x,' P_land=',f5.3,2x,' P_coll=',f6.3,2x,' P_a=',f6.3)
9     format('Total passes:           P_tot=',f6.3,2x,' P_land=',f5.3,2x,' P_coll=',f6.3,2x,' P_a=',f6.3)
20    format('written time slice (seconds) =',f9.3)
    end

    !subroutine traj(xm0, tet0, xbeg, nmax, nb1, nb2, pabs) 
    subroutine tracing(traj, nmax, nb1, nb2, pow, pabs) 
        use constants, only : tiny1
        use rt_parameters, only: eps, rrange, hdrob, nr, ipri, iw
        use dispersion_module, only: izn, yn3
        use dispersion_module, only: extd4, disp2
        use dispersion_module, only: find_all_roots, find_all_roots_simple
        use driver_module, only: im4, hrad, irs, iabsorp, iznzz, iwzz, irszz, rzz
        use driver_module, only: tetzz, xmzz
        use driver_module, only: driver2, driver4
        use dispersion_equation, only: ynz
        use trajectory_data

        implicit none
        class(Trajectory), intent(inout) :: traj
        real(wp), intent(inout)    :: pow
        real(wp), intent(in)    :: pabs
        integer,  intent(inout) :: nmax        
        integer,  intent(inout) :: nb1, nb2        
        integer :: num_roots
        !integer,  intent(in)    :: nomth, nomnz


        real(wp) :: xm0
        real(wp) :: tet0
        real(wp) :: xbeg
        integer :: nrefl
        integer :: irep
        integer :: irf, irf1
        integer :: ib2
        integer :: irs0        
        integer :: inak_saved 
        real(wp), parameter :: pgdop=0.02d0
        real(wp), parameter :: hmin=0.d-7 !sav2008, old hmin=1.d-7
        real(wp) :: eps0
        real(wp) :: rrange0, hdrob0, tet, xm, hr
        real(wp) :: xsav, xend,hsav, h1
        real(wp) :: ystart(2),yy(4)

        real(wp) :: xnr
        real(wp) :: ynz0, x1, x2, rexi, tetnew
        real(wp) :: xmnew, rnew, xnrnew
        real(wp) :: xnr_root(4)
        real(wp) :: pg1, pg2, pg3, pg4, pg

        ! copy initial parameters for a trajectory
        xm0  = traj%xmzap
        tet0 = traj%tetzap
        xbeg = traj%rzap
        yn3  = traj%yn3zap
        irs  = traj%irszap
        iw   = traj%iwzap
        izn  = traj%iznzap

        eps0=eps
        rrange0=rrange
        hdrob0=hdrob
        nrefl=0
        im4=0
        nb1=0
        nb2=0
        irep=0
        tet=tet0
        xm=xm0
        hr=1.d0/dble(nr+1) !sav2008
        hrad=hr
        !---------------------------------------
        ! find saving point and define
        ! parameters, depending on direction
        !---------------------------------------
  
  10    irf1=idnint(xbeg/hr)
        if (dabs(irf1*hr-xbeg).lt.tiny1)  then
            xsav=hr*irf1
        else
            irf=int(xbeg/hr)
            if (irs.eq.1)  xsav=hr*irf
            if (irs.eq.-1) xsav=hr*(irf+1)
        end if
        xend= 0.5d0 - 0.5d0*irs + tiny1*irs
        if (ipri.gt.2) write (*,*) 'xbeg-xend',xbeg,xend
        hsav = -hr*irs
        h1 = hsav
        !---------------------------------------
        ! solve eqs. starting from xbeg
        !---------------------------------------
        ystart(1) = tet
        ystart(2) = xm
        call driver2(ystart,xbeg,xend,xsav,hmin,h1, pow, pabs)
        tet = ystart(1)
        xm = ystart(2)
        ib2 = 0

        !---------------------------------------
        ! absorption
        !---------------------------------------
        if(iabsorp.ne.0) then
            if(ipri.gt.2) write (*,*)'in traj() iabsorp=',iabsorp
            nmax=nrefl
            !return
            goto 90
        end if
        if (xend.eq.xbeg) nb1=nb1+1
        !sav2008 20    continue
  
        !--------------------------------------------------------
        !  pass turning point
        !--------------------------------------------------------
        irs0=irs
        num_roots= find_all_roots(xend, xm, tet, xnr_root)
        !num_roots= find_all_roots_simple(xend, xm, tet, xnr_root)

        xnr = xnr_root(1)
        if (num_roots == 0) then
            print *,'no roots'
            pause
        endif

        ynz0 = ynz
  40    yy(1)=tet 
        yy(2)=xm
        yy(3)=xend
        yy(4)=xnr
        x1=0d0
        x2=1d+10
        rexi=xend
        inak_saved = current_trajectory%size !inak
        call driver4(yy,x1,x2,rexi,hmin, extd4)
        if(iabsorp.eq.-1) goto 90 !return !failed to turn

        tetnew = yy(1)
        xmnew  = yy(2)
        rnew   = yy(3)
        xnrnew = yy(4)

        if(ipri.gt.2) write (*,*) 'from r=',rexi,'to r=',rnew
 
        !---------------------------------------
        ! find mode
        !---------------------------------------

        !call disp2_iroot3(rnew, xmnew, tetnew, xnr_root)
        !print *, xnrnew
        num_roots= find_all_roots(rnew, xmnew, tetnew, xnr_root)
        !num_roots= find_all_roots_simple(rnew, xmnew, tetnew, xnr_root)

        pg1 = abs(xnrnew-xnr_root(1))
        pg2 = abs(xnrnew-xnr_root(2))
        pg3 = abs(xnrnew-xnr_root(3))
        pg4 = abs(xnrnew-xnr_root(4))

        pg = dmin1(pg1,pg2,pg3,pg4)
        if (dabs(pg/xnrnew).gt.pgdop) then
            !---------------------------------------------
            ! bad accuracy, continue with 4 equations
            !--------------------------------------------
            ib2=ib2+1
            nb2=nb2+1
            if (ib2.gt.4) then
                if (ipri.gt.1) write (*,*) 'error: cant leave 4 eqs'
                iabsorp=-1
                print *,'exit ib2.gt.4'
                !return
                goto 90
            end if
            eps=eps/5d0
            rrange=rrange*2d0
            hdrob=hdrob*2d0
            call current_trajectory%reset(inak_saved) ! костыль - восстанавливаю занчение счетчика точек траектории
            goto 40
        end if
        !-------------------------------------
        !          change wave type
        !-------------------------------------
        if (pg.ne.pg1) then
            if (pg.eq.pg2) izn=-izn
            if (pg.eq.pg3) iw=-iw
            if (pg.eq.pg4) iw=-iw
            if (pg.eq.pg4) izn=-izn
        end if
        if (irs0.ne.irs) nrefl=nrefl+1
        xbeg=rnew
        tet=tetnew
        xm=xmnew
        im4=1
        eps=eps0
        rrange=rrange0
        hdrob=hdrob0
        if(nrefl.lt.nmax) goto 10
        rzz=xbeg
        tetzz=tet
        xmzz=xm
        iznzz=izn
        iwzz=iw
        irszz=irs

        !---------------------------------------
        ! remember end point of trajectory
        !---------------------------------------
90      traj%rzap   = rzz
        traj%tetzap = tetzz
        traj%xmzap  = xmzz
        traj%yn3zap = yn3
        traj%iznzap = iznzz
        traj%iwzap  = iwzz
        traj%irszap = irszz
    end  


    

end module trajectory_module