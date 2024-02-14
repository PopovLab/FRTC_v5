      subroutine lhcd2017(outpe)
cc******************************************************************
cc   outj(i)  = LH driven current density, MA/m^2
cc   outpe(i) =LH power density (Landau+coll.) deposited into electrons, MW/m^3
cc   outpec(i) = LH power density (collisions) deposited into electrons, MW/m^3
cc   outpef(i) = LH power dens. dep. into el. by fast wave, MW/m^3
cc   dndt(i)  = d^2Jr1/dt^2/E, MA/m^2/sec^2/(V/m), ~runaway d(el.density)/dt/E
cc   djdt(i)  = dJr2/dt, time drivative of runaway current Jr2, MA/m^2/sec
cc   outpa(i)  = LH power density deposited into alphas, MW/m^3
cc   outda(i)  = (Na/Ne) relative density of alpha-particles
cc******************************************************************
      use approximation
      use utils
      use plasma
      use rt_parameters
      use maxwell  
      use spectrum_mod    
      use lhcd_module  
      implicit none
      integer i
      real*8 p_in,pe_p,pe_m,c_p,c_m
      real*8 vint
      include 'for/parameter.inc'
      include 'for/const.inc'
      include 'for/status.inc'
      real*8 outpe(NRD)
      real*8,dimension(:),allocatable:: outpep,outpem
      type(Spectrum) spectr
      type(Timer) my_timer
cc*********************************************************************
cc    Co-ordinates used in ray-tracing:
cc         (x-x0)/rm=r*cos(teta)-delta-gamma*sin^2(teta) !sav2008 -gamma
cc         (z-z0)/rm=ell*r*sin(teta)
cc    Definitions:
cc    (x0,z0) - magnetic axis position, centimeters
cc    rm      - minor radius in mid-plane, cenrimeters
cc    r(rho_ASTRA),delta(r),gamma(r),ell(r) - dimensionless functions
cc    rho_ASTRA=sqrt(Phi_tor/GP/BTOR)
cc    Interval for r:  0.<= r <=1.
cc*********************************************************************
      call my_timer%start('lhcd/lhcd_time.dat', time)
      print *, 'start lhcd2017 time=', time
      print *, 'NA1 NB1', NA1, NB1
      tcur=time
      outpe=zero
      p_in=dble(QLH)    ! input LH power, MW

      if(p_in.eq.zero) then 
            dij(:,:,:)=zero
            return
      end if

      call read_parameters('lhcd/ray_tracing.dat')
      call init_plasma(NA1,ABC,BTOR,RTOR,UPDWN,GP2,
     & AMETR,RHO,SHIF,ELON,TRIA,MU,NE,TE,TI,ZEF,UPL)

      call write_plasma(time)
      
      call write_lcms

      full_spectrum = read_spectrum('lhcd/spectrum.dat')
      full_spectrum%input_power = p_in

      pos_spectr = full_spectrum%get_positive_part()
      neg_spectr = full_spectrum%get_negative_part()
      
      call full_spectrum%normalization
      call full_spectrum%write('full_spectrum')
!!!!!!!!!!!!! starting ray-tracing !!!!!!!!!!!!!!!!!!!!!
      allocate(outpep(ngrid),outpem(ngrid))

!!positive spectrum:
      print *, 'positive spectrum'
      pe_p=zero
      outpep=zero
      if(pos_spectr%input_power > zero) then
            print *, 'spectrum_type',spectrum_type
            select case (spectrum_type)
            case (0)
                  spectr = make_spline_approximation(pos_spectr)
            case (1)
                  spectr = pos_spectr 
                  call spectr%normalization
            case (2)
                  spectr = pos_spectr
                  call spectr%normalization
            case (3)
                  print *, '2D spectrum'
                  stop                  
            end select
            call spectr%write('spectrum_pos')
            call ourlhcd2017(spectr, outpep,pe_p)
      else
            dij(:,:,1)=zero
      end if      

      if(pe_p.ne.zero) then
            c_p=vint(outpep,roc)
            if(c_p.ne.zero) then
                  do i=1,ngrid
                        outpep(i)=pe_p*outpep(i)/c_p
                  end do
            end if
      end if

!!negative spectrum:
       print *, 'negative spectrum'
       pe_m=zero
       outpem=zero       
       if(neg_spectr%input_power > zero) then        
            select case (spectrum_type)
            case (0)
                  spectr = make_spline_approximation(neg_spectr)
            case (1)
                  spectr = neg_spectr 
                  call spectr%normalization
            case (2)
                  spectr = neg_spectr
                  call spectr%normalization
            case (3)
                  print *, '2D spectrum'
                  stop
            end select            
            call spectr%write('spectrum_neg')            
            call ourlhcd2017(spectr, outpem,pe_m)              
       else
            dij(:,:,2)=zero
       endif     
       
       if(pe_m.ne.zero) then
            c_m=vint(outpem,roc)
            if(c_m.ne.zero) then
                  do i=1,ngrid
                        outpem(i)=pe_m*outpem(i)/c_m
                  end do
            end if
      end if

      do i=1,ngrid
       outpe(i)=outpep(i)+outpem(i)
      end do

      deallocate(outpep,outpem)
      call my_timer%stop_and_save
      print *, 'lhcd elapsed time', my_timer%elapsed_time
      !pause

      end
