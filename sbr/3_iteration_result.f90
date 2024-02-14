module iteration_result_mod
    use kind_module
    implicit none
    type IterationResult

        integer :: number
        !! iteration number 'iteration=',iterat
        integer :: spectr_direction
        !! 'ispectr=',ispectr
        real(wp) :: P_launched
        !!P_launched, MW=',plaun
        real(wp) :: P_landau
        !!'P_landau, MW=',ol
        real(wp) :: P_coll
        !! 'P_coll, MW=',oc
        real(wp) :: P_alph
        !!'P_alph, MW=',oa 
        real(wp) :: alphas_power
        !!'Alphas power, MW=',fuspow
        real(wp) :: P_fast
        !write(*,*) 'P_fast (landau+coll), MW=',of
        real(wp) :: P_lost
        !write(*,*) 'P_lost, MW=',plost/xsgs
        real(wp) :: P_not_accounted
        !write(*,*) 'P_not accounted, MW=',pnab/xsgs
        real(wp) :: P_landau_strong_absorption
        !write(*,*) 'P_landau (strong absorption), MW=',ppv1/xsgs
        real(wp) :: P_landau_weak_absorption
        !write(*,*) 'P_landau (weak absorption), MW=',ppv2/xsgs
        real(wp) :: P_turns
        !write(*,*) 'P_turns, MW=', psum4/xsgs
        real(wp) :: efficiency
        !write(*,*) 'efficiency, I(MA)/P(MW)=',oi/plaun !sav2008
        !call integral(1,nspl,rh,con,avedens) !sav2010
        real(wp) :: avedens
        real(wp) :: r0
        !write (*,*) '<Ne>, m^-3=',avedens*1.d19,' R, m=',r0*1.d-2
        real(wp) :: eta_eff
        !eta_eff=1.d17*avedens*r0*oi/plaun
        !write (*,*) 'eta_eff=<Ne>*R*I/P, A/(W*m^2)=',eta_eff !sav2010
        real(wp) :: residual 
        !! невязка 'nevyazka=', pchg        
    
    contains
        procedure :: print => iteration_result_print
        procedure :: save => iteration_result_save
    
    end type IterationResult
contains
    subroutine iteration_result_print(this)
        class(IterationResult), intent(in) :: this
        print *, ' ---------'
        print *, 'ITERATION:'
        print *, 'iteration=', this%number
        print *, 'ispectr=', this%spectr_direction
        print *, 'P_launched, MW=', this%P_launched
        print *, 'P_landau, MW=', this%P_landau
        print *, 'P_coll, MW=', this%P_coll
        print *, 'P_alph, MW=', this%P_alph 
        print *, 'Alphas power, MW=', this%alphas_power
        print *, 'P_fast (landau+coll), MW=', this%P_fast
        print *, 'P_lost, MW=', this%P_lost
        print *, 'P_not accounted, MW=', this%P_not_accounted
        print *, 'P_landau (strong absorption), MW=', this%P_landau_strong_absorption
        print *, 'P_landau (weak absorption), MW=', this%P_landau_weak_absorption
        print *, 'P_turns, MW=', this%P_turns
        print *, 'efficiency, I(MA)/P(MW)=', this%efficiency 
        !call integral(1,nspl,rh,con,avedens) !sav2010
        print *, '<Ne>, m^-3=',this%avedens,' R, m=',this%r0
        !eta_eff=1.d17*avedens*r0*oi/plaun
        print *, 'eta_eff=<Ne>*R*I/P, A/(W*m^2)=',this%eta_eff
        print *, 'nevyazka=', this%residual
    end subroutine iteration_result_print

    subroutine iteration_result_save(this, time_stamp)
        !! save Iteration Result to file
        class(IterationResult), intent(in) :: this
        real(wp), intent(in) :: time_stamp
        character(120) fname
        integer, parameter :: iu = 20
        write(fname,'("lhcd/rt-result/", f9.7,".dat")')  time_stamp
        print *, fname

        open(iu, file=fname,position="append")
        if (this%number == 1 ) then
            write (20, '(18A21)') 'iteration', 'direction','P_launched', &
            'P_landau', 'P_coll', 'P_alph', 'alphas_power', 'P_fast', 'P_lost', &
            'P_not_accounted', 'P_landau_strong_abs', 'P_landau_weak_abs', &
            'P_turns', 'efficiency', 'avedens', 'r0',  'eta_eff','residual'

        endif
        write (20, '(2(I21), 16(ES21.14))') this%number, this%spectr_direction, &
                    this%P_launched, this%P_landau, this%P_coll, this%P_alph, &
                    this%alphas_power, this%P_fast, this%P_lost, this%P_not_accounted, &
                    this%P_landau_strong_absorption, this%P_landau_weak_absorption, &
                    this%P_turns, this%efficiency, this%avedens, this%r0, &
                    this%eta_eff, this%residual
        close(iu)

    end subroutine     
end module iteration_result_mod