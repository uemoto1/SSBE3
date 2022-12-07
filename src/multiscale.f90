module multiscale
    implicit none
contains

subroutine multiscale_main(icomm)
    use mpi
    use communication
    use omp_lib
    use sbe_gs
    use sbe_solver
    use input_parameter
    use test
    use em_field
    use fdtd_weyl
    use math_constants, only: pi
    use phys_constants, only: cspeed_au
    use util
    implicit none
    integer, intent(in) :: icomm

    type(s_sbe_gs) :: gs
    type(s_sbe_bloch_solver), allocatable :: sbe(:)
    real(8) :: t
    real(8), allocatable :: Ac_ext_t(:, :)
    integer :: it, i
    integer :: nproc, irank, icomm_macro, nproc_macro, irank_macro
    real(8), allocatable :: Ac_macro(:, :)
    real(8), allocatable :: Jmat_macro(:, :)
    integer, allocatable :: itbl_macro_coord(:, :)
    integer :: nmacro, nmacro_max
    integer :: imacro_min, imacro_max
    integer :: ix, iy, iz, mt, imacro

    type(s_fdtd_system) :: fs
    type(ls_fdtd_weyl) :: fw
    character(256) :: tmp

    call comm_get_groupinfo(icomm, irank, nproc)

    ! FDTD setup
    fs%mg%nd = 1
    fs%mg%is(1) = 1 - abs(nxvac_m(1))
    fs%mg%ie(1) = nx_m + abs(nxvac_m(2))
    fs%mg%is(2) = 1 - abs(nyvac_m(1))
    fs%mg%ie(2) = ny_m + abs(nyvac_m(2))
    fs%mg%is(3) = 1 - abs(nzvac_m(1))
    fs%mg%ie(3) = nz_m + abs(nzvac_m(2))
    fs%mg%is_array(1:3) = fs%mg%is(1:3) - fs%mg%nd
    fs%mg%ie_array(1:3) = fs%mg%ie(1:3) + fs%mg%nd
    fs%hgs(1:3) = (/ hx_m, hy_m, hz_m /)
    fw%dt = dt
    call weyl_init(fs, fw)

    ! Prepare external pulse
    mt = max(nt, int(abs(nxvac_m(1)) * fs%hgs(1) / cspeed_au / dt))
    allocate(Ac_ext_t(1:3, -1:mt+1))
    call calc_Ac_ext_t(0.0d0, dt, 0, mt, Ac_ext_t)    
    call set_incident_field(mt, Ac_ext_t, fs, fw)

    ! Macropoint and media setup
    nmacro_max = nx_m * ny_m * nz_m
    allocate(itbl_macro_coord(3, nmacro_max))
    if (irank == 0) then
        call read_media_info(nmacro_max, itbl_macro_coord, nmacro, fw)
    end if
    call comm_bcast(itbl_macro_coord, icomm, 0)
    call comm_bcast(nmacro, icomm, 0)

    if (nmacro > 0) then
        if (0.0d0 < al(1)) al_vec1(1:3) = (/ al(1), 0.0d0, 0.0d0 /)
        if (0.0d0 < al(2)) al_vec2(1:3) = (/ 0.0d0, al(2), 0.0d0 /)
        if (0.0d0 < al(3)) al_vec3(1:3) = (/ 0.0d0, 0.0d0, al(3) /)
        if (nstate_sbe < 1) nstate_sbe = nstate

        ! Read ground state electronic system:
        call init_sbe_gs(gs, sysname, gs_directory, &
            & nkgrid, nstate, nelec, &
            & al_vec1, al_vec2, al_vec3, &
            & .false., icomm)
        
        ! Distribute
        call distribute_macropoints(irank, nmacro, nproc, imacro_min, imacro_max)
        icomm_macro = comm_create_group(icomm, imacro_min, 0)
        call comm_get_groupinfo(icomm_macro, irank_macro, nproc_macro) 
        allocate(Ac_macro(1:3, nmacro))
        allocate(Jmat_macro(1:3, nmacro))
        allocate(sbe(imacro_min:imacro_max))    

        ! Initialization of SBE solver and density matrix:
        do i = imacro_min, imacro_max
            call init_sbe(sbe(i), gs, nstate_sbe, icomm_macro)
        end do    

        if (irank_macro == 0) then
            do imacro = imacro_min, imacro_max
                write(tmp, "(a,a,a,i6.6,a)") trim(base_directory), trim(sysname), "_sbe_macro_",imacro,"_rt.data"
                open(1000+imacro, file=trim(tmp), action="write")
                write(1000+imacro, '(4a)') "# 1:Time[a.u.] 2:Ac_ext_x[a.u.] 3:Ac_ext_y[a.u.] 4:Ac_ext_z[a.u.] ", &
                    & "5:E_ext_x[a.u.] 6:E_ext_y[a.u.] 7:E_ext_z[a.u.] 8:Ac_tot_x[a.u.] ", &
                    & "9:Ac_tot_y[a.u.] 10:Ac_tot_z[a.u.] 11:E_tot_x[a.u.] 12:E_tot_y[a.u.] ", &
                    & "13:E_tot_z[a.u.]  14:Jm_x[a.u.] 15:Jm_y[a.u.] 16:Jm_z[a.u.]"
            end do
        end if
    end if

    do it = 1, nt
        t = dt * it

        if (nmacro > 0) then
            do imacro = imacro_min, imacro_max
                ix = itbl_macro_coord(1, imacro)
                iy = itbl_macro_coord(2, imacro)
                iz = itbl_macro_coord(3, imacro)
                Ac_macro(1:3, imacro) = fw%vec_Ac_new%v(1:3, ix, iy, iz)
                call dt_evolve_bloch(sbe(imacro), gs, Ac_macro(1:3, imacro), dt)
                call calc_current_bloch(sbe(imacro), gs, Ac_macro(1:3, imacro), Jmat_macro(1:3, imacro), icomm_macro)
                fw%vec_j_em%v(:, ix, iy, iz) = -Jmat_macro(:, imacro)

                if (mod(it, 10) == 0) then
                    if (irank_macro == 0) then
                        write(1000+imacro, '(f12.6,15(es24.15e3))') t, Ac_macro(1:3, imacro), fw%vec_e%v(1:3, ix, iy, iz), &
                            & Ac_macro(1:3, imacro), fw%vec_e%v(1:3, ix, iy, iz), Jmat_macro(:, imacro)
                    end if    
                end if
            end do
        end if

        call weyl_calc(fs, fw)
        
        if ((mod(it, 100) == 0)) then
            if (irank == 0) call write_Ac_field(it, fs, fw)
        end if
    end do

    if (nmacro > 0) then
        if (irank_macro == 0) then
            do imacro = imacro_min, imacro_max
                close(1000+imacro)
            end do
        end if
    end if

    return    


contains

end subroutine multiscale_main


subroutine distribute_macropoints(irank, nmacro, nproc, imacro_min, imacro_max)
    use util, only: split_range
    use communication, only: comm_create_group, comm_get_groupinfo
    implicit none
    integer, intent(in) :: irank
    integer, intent(in) :: nmacro
    integer, intent(in) :: nproc
    integer, intent(out) :: imacro_min
    integer, intent(out) :: imacro_max

    integer :: i
    integer :: itbl_macro_min(0:nproc-1)
    integer :: itbl_macro_max(0:nproc-1)
    integer :: itbl_rank_min(1:nmacro)
    integer :: itbl_rank_max(1:nmacro)

    if (nproc <= nmacro) then
        call split_range(1, nmacro, nproc, itbl_macro_min, itbl_macro_max)
        imacro_min = itbl_macro_min(irank)
        imacro_max = itbl_macro_max(irank)    
    else
        call split_range(0, nproc-1, nmacro, itbl_rank_min, itbl_rank_max)
        do i = 1, nmacro
            if ((itbl_rank_min(i) <= irank) .and. (irank <= itbl_rank_max(i))) then
                imacro_min = i
                imacro_max = i
            end if
        end do
    end if

    return
end subroutine distribute_macropoints


subroutine read_media_info(nmacro_max, itbl_macro_coord, nmacro, fw)
    use input_parameter, only: file_ms_shape, epsilon_em, nx_m, ny_m, nz_m
    use fdtd_weyl, only: ls_fdtd_weyl
    implicit none
    integer, intent(in) :: nmacro_max
    integer, intent(out) :: itbl_macro_coord(1:3, nmacro_max)
    integer, intent(out) :: nmacro
    type(ls_fdtd_weyl), intent(inout) :: fw

    integer :: i, imacro, n, ix, iy, iz, itype

    fw%epsilon%f(:, :, :) = 1.0d0
    imacro = 0
    if (len_trim(file_ms_shape) > 0) then
        open(99, file=trim(file_ms_shape), action="read")
        read(99, *) n
        do i = 1, n
            read(99, *) ix, iy, iz, itype
            if (ix < 1) stop "Error: invalid range!"
            if (iy < 1) stop "Error: invalid range!"
            if (iz < 1) stop "Error: invalid range!"
            if (ix > nx_m) stop "Error: invalid range!"
            if (iy > ny_m) stop "Error: invalid range!"
            if (iz > nz_m) stop "Error: invalid range!"
            if (itype > 0) then
                fw%epsilon%f(ix, iy, iz) = epsilon_em(itype)
            else if (itype == 0) then
                fw%epsilon%f(ix, iy, iz) = 1.0d0
            else
                imacro = imacro + 1
                if (imacro > nmacro_max) stop "Error: number of macropoints is too large!" 
                itbl_macro_coord(1:3, imacro) = (/ ix, iy, iz /)
            end if
        end do
        close(99)
    else
        do iz = 1, nz_m
        do iy = 1, ny_m
        do ix = 1, nx_m
            imacro = imacro + 1
            if (imacro > nmacro_max) stop "Error: number of macropoints is too large!" 
            itbl_macro_coord(1:3, imacro) = (/ ix, iy, iz /)
        end do
        end do
        end do
    end if
    nmacro = imacro
end subroutine 

subroutine set_incident_field(mt, Ac, fs, fw)
    use fdtd_weyl, only: s_fdtd_system, ls_fdtd_weyl
    use phys_constants, only: cspeed_au
    implicit none
    integer, intent(in) :: mt
    real(8), intent(in) :: Ac(1:3, 0:mt)
    type(s_fdtd_system), intent(in) :: fs
    type(ls_fdtd_weyl), intent(inout) :: fw

    integer :: ix, iy, iz, it
    real(8) :: r_it, w, x

    do ix = fs%mg%is_array(1), 0
        do iy = fs%mg%is_array(2), fs%mg%ie_array(2)
            do iz = fs%mg%is_array(3), fs%mg%ie_array(3)
                x = ix * fs%hgs(1)
                r_it = (-x / cspeed_au / fw%dt)
                it = int(r_it)
                if ((0 <= it) .and. (it < mt)) then
                    w = r_it - it
                    fw%vec_Ac_new%v(:, ix, iy, iz) = w * Ac(:, it+1) + (1.0d0-w) * Ac(:, it)
                end if
                r_it = r_it - 1
                it = int(r_it)
                if ((0 <= it) .and. (it < mt)) then
                    w = r_it - it
                    fw%vec_Ac%v(:, ix, iy, iz) = w * Ac(:, it+1) + (1.0d0-w) * Ac(:, it)
                end if
            end do
        end do
    end do
end subroutine set_incident_field

subroutine write_Ac_field(iit, fs, fw)
    use input_parameter, only: sysname, out_ms_ix, out_ms_iy, out_ms_iz
    use fdtd_weyl, only: s_fdtd_system, ls_fdtd_weyl
    use phys_constants, only: cspeed_au
    implicit none
    integer, intent(in) :: iit
    character(256) :: file_Ac_data
    type(s_fdtd_system), intent(in) :: fs
    type(ls_fdtd_weyl), intent(in) :: fw
    integer :: ix, iy, iz

    write(file_Ac_data, "(a, a,i6.6,a)") sysname, "_Ac_", iit, ".data"
    open(999, file=trim(file_Ac_data), action="write")
    write(999, '(a)') "# Multiscale TDDFT calculation"
    write(999, '(a)') "# IX, IY, IZ: FDTD Grid index"
    write(999, '(a)') "# x, y, z: Coordinates"
    write(999, '(a)') "# Ac: Vector potential field"
    write(999, '(a)') "# E: Electric field"
    write(999, '(a)') "# J_em: Electromagnetic current density"
    write(999, '("#",99(1X,I0,":",A,"[",A,"]"))') &
        & 1, "IX", "none", &
        & 2, "IY", "none", &
        & 3, "IZ", "none", &
        & 4, "Ac_x", "[au]", &
        & 5, "Ac_y", "[au]", &
        & 6, "Ac_z", "[au]", &
        & 7, "E_x", "[au]", &
        & 8, "E_y", "[au]", &
        & 9, "E_z", "[au]", &
        & 10, "B_x", "a.u.", &
        & 11, "B_y", "a.u.", &
        & 12, "B_z", "a.u.", &
        & 13, "Jem_x", "[au]", &
        & 14, "Jem_y", "[au]", &
        & 15, "Jem_z", "[au]", &
        & 16, "E_em", "[au]" // "/vol", &
        & 17, "E_abs", "[au]" //  "/vol"
    do iz = max(fs%mg%is(3), out_ms_iz(1)), min(fs%mg%ie(3), out_ms_iz(2))
    do iy = max(fs%mg%is(2), out_ms_iy(1)), min(fs%mg%ie(2), out_ms_iy(2))
    do ix = max(fs%mg%is(1), out_ms_ix(1)), min(fs%mg%ie(1), out_ms_ix(2))
        write(999, "(3i9,99es25.15e4)") &
        ix, iy, iz, &
        fw%vec_Ac%v(:, ix, iy, iz), &
        fw%vec_e%v(:, ix, iy, iz), &
        fw%vec_h%v(:, ix, iy, iz), &
        fw%vec_j_em%v(:, ix, iy, iz)
    end do
    end do
    end do
    close(999)
end subroutine write_Ac_field



end module multiscale
