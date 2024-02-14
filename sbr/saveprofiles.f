      subroutine saveprofiles
       use plasma
       use rt_parameters
       use maxwell
      implicit none
      include 'for/parameter.inc'
      include 'for/const.inc'
      include 'for/status.inc'
      logical, save :: first_time = .TRUE.

      if (first_time) then
            call read_parameters('lhcd/ray_tracing.dat')
      end if

      call init_plasma(NA1,ABC,BTOR,RTOR,UPDWN,GP2,
     & AMETR,RHO,SHIF,ELON,TRIA,MU,NE,TE,TI,ZEF,UPL)

      call calc_enorm
     
      if (first_time) then
            call init_maxwell
            first_time = .FALSE.
      end if

      end

