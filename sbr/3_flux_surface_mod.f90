module FluxSurface_mod
    !! все что связанно с магнитными поверхностями
    use kind_module
    type FluxSurface 
    !! класс магнитной поверхности
        integer       :: index
        !! номер магнитной поверхности
        real(wp) :: r
        !! радиус
        real(wp) :: vmax
        !! vmax=cltn/vto
        real(wp) :: vt
        !! наверно тепловая скорость электронов????? vt=fvt(r)
        integer       :: ipt
        !! размер vgrid
        real(wp), allocatable :: vgrid(:)
        !! 
        real(wp), allocatable :: vr_grid(:)
        !! бываший vrj
        real(wp), allocatable :: diffusion(:)
        !! бывший dijk(i,j,k) или dj(i)
        !   complex         :: inst_field1
        contains
        !procedure :: set   => set_e
        !procedure :: print => e_print

    end type FluxSurface  

end module FluxSurface_mod