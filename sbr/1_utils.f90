module utils
    use kind_module
    implicit none
    type Timer
        real(wp) :: start_time
        real(wp) :: end_time
        real(wp) :: elapsed_time
        real(wp) :: plasma_time
        character(120) :: file_name
    contains
        procedure :: start => Time_start
        procedure :: stop  => Time_stop
        procedure :: stop_and_save  => Time_stop_and_save
    end type Timer

    contains

    subroutine Time_start(this, fn, pt)
        implicit none
        class(Timer),  intent(inout) :: this
        character(*), intent(in)   :: fn
        real(wp), intent(in)         :: pt
        this%plasma_time = pt
        this%start_time = sys_time()
        this%file_name = fn
    end subroutine

    subroutine Time_stop(this)
        implicit none
        class(Timer),  intent(inout) :: this
        this%end_time = sys_time()
        this%elapsed_time = this%end_time - this%start_time
    end subroutine

    subroutine Time_stop_and_save(this)
        implicit none
        class(Timer),  intent(inout) :: this
        integer io
        call this%stop
		open(newunit=io, file=this%file_name, position="append")
		write(io,'(2f22.14)') this%plasma_time, this%elapsed_time
    	close(io)
    end subroutine

    function sys_time()
    ! ** return system time
    implicit none
        real(wp) sys_time
        integer count, count_rate, count_max
        call system_clock(count, count_rate, count_max)
        sys_time = count*1.0/count_rate
    return
    end

end  module utils