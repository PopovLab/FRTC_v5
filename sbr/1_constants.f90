module constants 
    !! модуль с математическими и физическими константами
    !! nt_001 тестовый комментарий
    use kind_module
    implicit none

    real(wp), parameter :: zero = 0.0_wp 
    real(wp), parameter :: one  = 1.0_wp
    real(wp), parameter :: two  = 2.0_wp
    real(wp), parameter :: one_third = 1.0_wp/3.0_wp

    real(wp), parameter :: tiny = 1.e-100_wp
    real(wp), parameter :: tiny1  = 1e-7_wp
    

    real(wp), parameter :: pi = acos(-one) !! число Пи = 3.1415....
    real(wp), parameter :: pi2 = 2.0_wp*pi
    real(wp), parameter :: pi4 = 4.0_wp*pi
    real(wp), parameter :: piq = sqrt(pi)
    real(wp), parameter :: pisqrt = sqrt(pi)
    real(wp), parameter :: pi2sqrt = sqrt(pi2)
    !pi2sqrt=2.506628274631d0,pisqrt=1.77245385090552d0)
    real(wp), parameter :: talfa  = 3.5_wp    !! alpha particles' birth energy, MeV
    real(wp), parameter :: zalfa  = 2.0_wp    !! alpha particles' electrical charge
    real(wp), parameter :: xmalfa = 4.0_wp  !! alpha particles' atomic mass

    real(wp), parameter :: clt = 3.0e+10_wp  !! скорость света

    real(wp), parameter :: pme  = 9.11e-28_wp
    real(wp), parameter :: pme_e  = 9.11e-28
    real(wp), parameter :: pqe  = 4.803e-10_wp
    real(wp), parameter :: xlog = 16.0_wp + dlog(16.0_wp)
    real(wp), parameter :: c0 = sqrt(pi4*pqe**2/pme)
    real(wp), parameter :: c1 = pqe/pme/clt
    real(wp), parameter :: xsgs = 1e+13_wp
    !! 1MW = 1e13 erg/s
    real(wp), parameter :: xwtt = 1e-7_wp

    real(wp), parameter :: cnst1 = 0.2965924106e-6_wp
    !! cnst1=(m_e/m_p)**2, CGS
    real(wp), parameter :: cnst2 = 0.359680922e-35_wp
    !! cnst2=(m_e/e)**2,  CGS

contains
    
    subroutine show_constants() 

        print *, '---------------------------------------' 
        print *, "zero       = ", zero          
        print *, "one        = ", one     
        print *, "two        = ", two
        print *, "one_third  = ", one_third
        print *, "tiny       = ", tiny     
        print *, "tiny1      = ", tiny1     
        print *
        print *, "4*ATAN     = ", 4.d0*datan(1.d0)
        print *, "pi         = ", pi
        print *, "pi2        = ", pi2     
        print *, "pi4        = ", pi4     
        print *, "piq        = ", piq             
        print *
        print *, "talfa      = ", talfa   
        print *, "zalfa      = ", zalfa     
        print *, "xmalfa     = ", xmalfa     
        print *, "clt        = ", clt     
        print *     
        print *, "pme_e      = ", pme_e
        print *, "pme        = ", pme     
        print *, "9.11d-28   = ", 9.11d-28     
       
        print *, "pqe        = ", pqe   
        print *, "xlog       = ", xlog     
        print *, "c0         = ", c0   
        print *, "c1         = ", c1     
        print *, "xsgs       = ", xsgs
        print *, "xwtt       = ", xwtt
        print *
        print *, "cnst1  = ", cnst1
        print *, "cnst2  = ", cnst2
        print *, '---------------------------------------'

    end subroutine show_constants

end module constants 