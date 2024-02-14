        subroutine drivencurrent(outj,sigmaj)
! ******************************************************************
!   outj(i)  = LH driven current density, MA/m^2
!   dndt(i)  = d^2Jr1/dt^2/E, MA/m^2/sec^2/(V/m), ~runaway d(el.density)/dt/E
!   djdt(i)  = dJr2/dt, time drivative of runaway current Jr2, MA/m^2/sec
!   outjrun(i)  = LH driven runaway current density, MA/m^2
!   outnerun(i) = runaway electron density/10^19 m^-3
! ******************************************************************
            use utils
          implicit none
          include 'for/parameter.inc'
          include 'for/const.inc'
          include 'for/status.inc'
          real*8 outj(NRD),sigmaj(NRD)
          type(Timer) my_timer
! RTOR    major radius
! ROC     effective minor radius
! UPL(*)  toroidal loop voltage
! NRD     501 - Maximum size of the radial grid
! NA1     Edge grid point number: ROC=RHO(NA1)
! GP2     2*Pi
          call my_timer%start('lhcd/drivencurrent_time.dat', time)
          call drivencurrent95(outj, sigmaj, UPL, NRD, NA1, TIME, TAU, 
     & ROC, RTOR, GP2)
             
            
          call my_timer%stop_and_save
          print *, 'driven current elapsed time', my_timer%elapsed_time
          !pause
        end    

        double precision function aiint(arr,yr)
        ! aiint:	integral {0,r} of current density
        ! 	aiint=integral {0,r} (arr/ipol**2)dv*ipol/(gp2*ro)
        ! 	only radial dependent array may be a parameter of the function
        ! 	examples:
        !    out_iint(cu)	!radial profile of toroidal current
        !    out_iint(cd,ro)    !toroidal driven current inside {0,ro}
        !    out_iint(cub)      !total toroidal current =iint(cu,roc); (=ipl)
        !			(pereverzev 23-oct-99)
! NA        = NA1-1
! RHO(*)    main magnetic surface label
! IPOL(*)   normalized poloidal current
! G33(*)
! HRO       Radial grid step in the variable ro        
! HROA      Edge step of the radial grid: RHO(NA1)-RHO(NA)
              implicit none
              integer j,jk
              double precision arr(*),yr,dr,ya
              include 'for/parameter.inc'
              include 'for/const.inc'
              include 'for/status.inc'
              aiint=0.d0
              if(yr.le.0.d0) return
              jk=yr/hro+1.d0-1.d-4
              if(jk.gt.na) jk=na
              ya=0.d0
              do j=1,jk
               aiint=aiint+ya
               ya=arr(j)*rho(j)/(g33(j)*ipol(j)**3)
               dr=yr-jk*hro+hro
              end do
              if(jk.ge.na) then
               dr=min(dr,0.5d0*(hro+hroa))
              endif
              aiint=gp2*ipol(jk)*(hro*aiint+dr*ya)
              end
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        