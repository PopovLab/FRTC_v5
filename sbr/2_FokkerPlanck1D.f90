module FokkerPlanck1D_mod ! the module name defines the namespace
    !! модуль содержит функции для решения одномерного уравнения Фоккер-Планка
    use, intrinsic :: iso_fortran_env, only: sp=>real32, dp=>real64
    use savelyev_solver_module
    implicit none
    type FokkerPlanck1D 
        !! solver of FP eq
        !integer          :: direction = 0
        !- direction
        real(dp)         :: enorm = 0
        !! электрическое поле
        real(dp)         :: v_lim = 0
        !! верхняя граница скорости электронов
        real(dp), allocatable         :: v(:)
        !! сетка скоростей
        real(dp), allocatable         :: f(:)
        !! распределение
        integer         :: i0 = 0 
        !! size of distribution grid
        real(dp)        :: alfa2  = 0 
        !! поле со знаком
        integer         :: n = 0 
        !! size of local grid        
        real(dp)        :: h  = 0 
        !! step of local grid 
        real(dp), allocatable ::  d1(:), d2(:), d3(:)
        !! диффузия
    contains
        procedure :: print => FokkerPlanck1D_print
        procedure :: solve_time_step => FokkerPlanck1D_solve_time_step
        procedure :: init_zero_diffusion => FokkerPlanck1D_init_zero_diffusion
        procedure :: init_diffusion => FokkerPlanck1D_init_diffusion
    end type FokkerPlanck1D   

    interface FokkerPlanck1D
        module procedure :: FokkerPlanck1D_constructor
    end interface FokkerPlanck1D

    contains

    function FokkerPlanck1D_constructor(e, v_lim, v, f) result(this)
        !! конструктор для FokkerPlanck1D
        implicit none
        type(FokkerPlanck1D) :: this
        real(dp), value :: e, v_lim, v(:), f(:)
        integer  :: n
        real(dp) :: h
        real(dp), parameter :: h0 = 0.1d0
        !this%inst_field1 = cmplx(0.,0.) 
        this%enorm     = abs(e)
        this%v_lim = v_lim
        this%v = v
        this%f = f
        this%i0 = size(v)
        this%alfa2 = e
        n = v_lim/h0-1
        h = v_lim/dble(n+1)
        if (h.gt.h0) then
            n = n+1
            h = v_lim/dble(n+1)
        end if
        this%n = n
        this%h = h
    end function FokkerPlanck1D_constructor

    subroutine FokkerPlanck1D_print(this)
        class(FokkerPlanck1D), intent(in) :: this

        print *, 'e = ', this%enorm, 'i0 =', this%i0
    end subroutine FokkerPlanck1D_print

    subroutine FokkerPlanck1D_init_zero_diffusion(this)
      implicit none
      class(FokkerPlanck1D), intent(inout) :: this
      integer :: n
      n = this%n
      allocate(this%d1(n+1),this%d2(n+1),this%d3(n+1))
      this%d1(:)=0d0
      this%d2(:)=0d0
      this%d3(:)=0d0
    end subroutine FokkerPlanck1D_init_zero_diffusion


    subroutine FokkerPlanck1D_init_diffusion(this, dif)
        !! инициализация диффузии для схемы савельева
        use lock_module
        implicit none
        class(FokkerPlanck1D), intent(inout) :: this
        integer :: n
        real(dp), dimension(:), intent(in) ::  dif
        real(dp), dimension(:), allocatable :: xx
        real(dp) h
        integer :: i0
        integer i, klo, khi, ierr, klo1, khi1
        integer klo2, klo3, khi2, khi3, ierr1, ierr2, ierr3
        n = this%n
        h = this%h
        allocate(this%d1(n+1),this%d2(n+1),this%d3(n+1))    
        i0 = this%i0
    
        allocate(xx(n+1))
        do i=1,n+1
            xx(i)=h/2.d0+h*dble(i-1) !+shift
        end do
    
        do i=1,n+1
            call lock(this%v, i0, xx(i), klo1, khi1, ierr1)
            call lock(this%v, i0, xx(i)-h/2d0, klo2, khi2, ierr2)
            call lock(this%v, i0, xx(i)+h/2d0, klo3, khi3, ierr3)
            if(ierr1.eq.1) then
                write(*,*)'lock error in finction d2(x)'
                write(*,*)'j=', 123,' v=', xx(i)
                write(*,*)'klo1=', klo1, 'khi1=', khi1, 'i=', i
                write(*,*)'vj(1)=', this%v(1),' vj(i0)=', this%v(i0)
                pause
                stop
            end if
            if(ierr2.eq.1) then
                write(*,*)'lock error in finction d2(x)'
                write(*,*)'j=', 123, ' v=', xx(i)
                write(*,*)'klo2=', klo2, 'khi2=', khi2, 'i=',i
                write(*,*)'vj(1)=', this%v(1), ' vj(i0)=', this%v(i0)
                pause
                stop
            end if
            if(ierr3.eq.1) then
                write(*,*)'lock error in finction d2(x)'
                write(*,*)'j=', 123, ' v=', xx(i)
                write(*,*)'klo3=', klo3, 'khi3=', khi3, 'i=',i
                write(*,*)'vj(1)=', this%v(1), ' vj(i0)=', this%v(i0)
                pause
                stop
            end if
            this%d1(i) = dif(klo1)
            this%d2(i) = dif(klo2)
            this%d3(i) = dif(klo3)	
        end do
      end subroutine FokkerPlanck1D_init_diffusion

    subroutine FokkerPlanck1D_solve_time_step(this, dt, nt)
        use lock_module
        implicit none
        class(FokkerPlanck1D), intent(inout) :: this
        
        integer, intent(in) :: nt
        real(dp), intent(in)  :: dt


        !real*8, intent (inout), optional :: dfj0(:)

        real(dp), parameter :: zero=0.d0
        real(dp)  y(this%n+2),x(this%n+2)
        real(dp), dimension(:), allocatable:: fj, dfj,  givi
        integer i, ii, it, ibeg, klo, khi, ierr, klo1, khi1
        real(dp) shift, ybeg, yend, tend, dff

         !!!!!! grid !!!!!!!!!
        !!  shift=h*0.1d0 !0.01d0
        do i=1, this%n+2
            x(i)=this%h*dble(i-1) !+shift
        end do
  
        do i=1, this%n+1
            call lock(this%v,this%i0,x(i+1),klo,khi,ierr)
            if(ierr.eq.1) then
                write(*,*)'lock error #1 in finction fokkerplanck'
                write(*,*)'j=', 123, ' v=', x(i+1)
                write(*,*)'vj(1)=',this%v(1),' vj(i0)=',this%v(this%i0)
                pause
                stop
            end if
            call linf(this%v,this%f,x(i+1),y(i),klo,khi)
        end do
  
        ybeg = this%f(1)  !boundary conditions
        yend = this%f(this%i0) !zero
        !print *, ' yend =', yend 
        !!!!!!!!!!!!   solve problem   !!!!!!!!!!!!!!!!!!!!!!!!!!
        call savelyev_solver(this%alfa2, nt, this%h, dt, this%n, ybeg, yend, this%d1, this%d2, this%d3, y)
        !call teplova_khavin_solver(this%alfa2, nt, this%h, dt, this%n, ybeg, yend, this%d1,this%d2,this%d3, y)
        allocate(fj(this%n+2))
        fj(1)=ybeg
        fj(this%n+2)=yend
        do i=1,this%n
            fj(i+1)=y(i)
        end do

        do i=2,this%i0-1
          if(this%v(i).lt.this%v_lim) then
              call lock(x,this%n+2,this%v(i),klo,khi,ierr)
              if(ierr.eq.1) then
                  write(*,*)'lock error #2 in finction fokkerplanck'
                  write(*,*)'j=', 123 ,' vij=',this%v(i)
                  write(*,*)'x(1)=', x(1),' x(n+2)=',x(this%n+2)
                  pause
                  stop   
              end if
              call linf(x,fj,this%v(i),this%f(i),klo,khi)
          else
            this%f(i)=zero
          end if
      end do
      deallocate(fj)
  
      !if (present(dfj0)) then
      !    call burying_procedure(vj, fj0, dfj0)
      !else 
      !    call burying_procedure(vj, fj0)
      !end if
  
  
    end subroutine FokkerPlanck1D_solve_time_step


    subroutine burying_procedure(v, f0, df0)
        !! процедура закапывания
        use math_module
        implicit none
        real*8,  intent(in)     :: v(:)        
        real*8,  intent(inout)  :: f0(:)
        real*8,  intent (inout), optional :: df0(:)
        integer i, ii,  i0, ibeg
        real*8, allocatable  :: f(:), df(:)
        real*8 fout1, fout2
    
        i0 = size(f0)
        allocate(f(i0), df(i0))
        
        f(:)=f0(:)
        df(:)=0d0
    
        do i=2, i0-1
            df(i)=0.5d0*(f(i+1)-f(i-1))/v(2)
        end do
        df(1)=0d0
        df(i0)=(f(i0)-f(i0-1))/v(2)
    
    !   сдвиг расределения вправо. зачем-то ???
        ii=0
        ibeg=0
        do i=i0-1,1,-1
            if(df(i).gt.0d0) then
    !          write(*,*) '#1 positive derivs'
    !          write(*,*) '#1 df>0: i,j,k=',i,j,k
    !          write(*,*) '#1 dfj(i),i,j,k=',dfj(i),i,j,k
    !          write(*,*)
                f0(i)=f0(i+1)
                if (present(df0)) then
                    df0(i)=df0(i+1)
                end if
                ii=i
            end if
            if(f0(i).lt.f0(i+1)) then 
                f0(i)=f0(i+1)
                if (present(df0)) then
                    df0(i)=df0(i+1)
                end if            
                ii=i
            end if
        end do
    
        if(ibeg.gt.0) then
            call integral(ibeg,i0,v,f,fout1)
            f(:) = f0(:)
            if (present(df0)) then
                df(:) = df0(:)
            end if
            
            call integral(ibeg,i0,v,f,fout2)
            f0(ibeg:i0) = f(ibeg:i0)*fout1/fout2
            if (present(df0)) then
                df0(ibeg:i0) = df(ibeg:i0)*fout1/fout2
            end if            
    !      write(*,*)'#1 j,k,ibeg=',j,k,ibeg
    !      write(*,*)'#1 v(ibeg)=',vj(ibeg),' f1/f2=',fout1/fout2
        end if
    
        deallocate(f,df)
        ibeg=ii
    end subroutine
 end module FokkerPlanck1D_mod