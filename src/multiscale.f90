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
    real(8), allocatable :: Jmat_macro_tmp(:, :), Jmat_macro(:, :)
    integer, allocatable :: itbl_macro_coord(:, :)
    integer :: nmacro, nmacro_max
    integer :: imacro_min, imacro_max
    integer :: ix, iy, iz, mt, imacro, iobs
    real(8) :: jmat(3)

    type(s_fdtd_system) :: fs
    type(ls_fdtd_weyl) :: fw
    character(256) :: tmp

    logical :: flag_1d_model

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

    flag_1d_model = ((ny_m == 1) .and. (nz_m == 1))

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
            & num_kgrid, nstate, nelec, &
            & al_vec1, al_vec2, al_vec3, &
            & .false., icomm)
        
        ! Distribute
        call distribute_macropoints(irank, nmacro, nproc, imacro_min, imacro_max)
        icomm_macro = comm_create_group(icomm, imacro_min, 0)
        call comm_get_groupinfo(icomm_macro, irank_macro, nproc_macro) 
        allocate(Ac_macro(1:3, nmacro))
        allocate(Jmat_macro_tmp(1:3, nmacro))
        allocate(Jmat_macro(1:3, nmacro))
        allocate(sbe(imacro_min:imacro_max))    

        ! Initialization of SBE solver and density matrix:
        do i = imacro_min, imacro_max
            call init_sbe(sbe(i), gs, nstate_sbe, icomm_macro)
        end do
    end if

    if (nmacro > 0) then
        if (irank == 0) then
            write(tmp, "(a,a,a,a)") "mkdir ", trim(base_directory), trim(sysname), "_sbe_RT_Ac"
            call system(trim(tmp))
            write(tmp, "(a,a,a,a)") "mkdir ", trim(base_directory), trim(sysname), "_sbe_m"
            call system(trim(tmp))
            do imacro = 1, nmacro
                write(tmp, "(a,a,a,a,i6.6)") "mkdir ", trim(base_directory), trim(sysname), "_sbe_m/m", imacro
                call system(trim(tmp))
            end do
        end if
    end if

    call comm_sync_all(icomm)

    if (nmacro > 0) then
        if (irank == 0) then
            if (flag_1d_model) then
                ! _sbe_wave.data
                write(tmp, "(a,a,a)") trim(base_directory), trim(sysname), "_sbe_wave.data"
                open(888, file=trim(tmp), action="write")
                write(888, '("#",99(1X,I0,":",A,"[",A,"]"))') &
                    & 1, "Time", "[a.u.]", &
                    & 2, "E_inc_x", "[a.u.]", &
                    & 3, "E_inc_y", "[a.u.]", &
                    & 4, "E_inc_z", "[a.u.]", &
                    & 5, "E_ref_x", "[a.u.]", &
                    & 6, "E_ref_y", "[a.u.]", &
                    & 7, "E_ref_z", "[a.u.]", &
                    & 8, "E_tra_x", "[a.u.]", &
                    & 9, "E_tra_y", "[a.u.]", &
                    & 10, "E_tra_z", "[a.u.]"
            end if
            ! _obs_sbe_rt.data
            do iobs = 1, obs_num_em
                write(tmp, "(a,a,a,i3.3,a)") trim(base_directory), trim(sysname), "_sbe_obs_", obs_num_em, "_at_point_rt.data"
                open(100+iobs, file=trim(tmp), action="write")
                write(100+iobs,'("#",99(1X,I0,":",A))') &
                1, "Time[a.u.]",                &
                2, "E_x[a.u.]",                 &
                3, "E_y[a.u.]",                 &
                4, "E_z[a.u.]",                 &
                5, "H_x[a.u.]",                 &
                6, "H_y[a.u.]",                 &
                7, "H_z[a.u.]"
            end do
        end if
        if (irank_macro == 0) then
            ! _sbe_rt.data
            do imacro = imacro_min, imacro_max
                write(tmp, "(a,a,a,i6.6,a,a,a)") trim(base_directory), trim(sysname), "_sbe_m/m", imacro, &
                    &  "/", trim(sysname), "_sbe_rt.data"
                ! write(*,*) trim(tmp)
                open(1000+imacro, file=trim(tmp), action="write")
                write(1000+imacro, '(a)') "# Real time calculation:"
                write(1000+imacro, '(a)') "# Ac_ext: External vector potential field"
                write(1000+imacro, '(a)') "# E_ext: External electric field"
                write(1000+imacro, '(a)') "# Ac_tot: Total vector potential field"
                write(1000+imacro, '(a)') "# E_tot: Total electric field"
                write(1000+imacro, '(a)') "# Jm: Matter current density (electrons)"
                write(1000+imacro, '(4a)') "# 1:Time[a.u.] 2:Ac_ext_x[a.u.] 3:Ac_ext_y[a.u.] 4:Ac_ext_z[a.u.] ", &
                    & "5:E_ext_x[a.u.] 6:E_ext_y[a.u.] 7:E_ext_z[a.u.] 8:Ac_tot_x[a.u.] ", &
                    & "9:Ac_tot_y[a.u.] 10:Ac_tot_z[a.u.] 11:E_tot_x[a.u.] 12:E_tot_y[a.u.] ", &
                    & "13:E_tot_z[a.u.]  14:Jm_x[a.u.] 15:Jm_y[a.u.] 16:Jm_z[a.u.]"
            end do
        end if
    end if

    call comm_sync_all(icomm)

    do it = 1, nt
        t = dt * it

        if (nmacro > 0) then
            ! Get vector potential field
            if (irank == 0) then
                do imacro = 1, nmacro
                    ix = itbl_macro_coord(1, imacro)
                    iy = itbl_macro_coord(2, imacro)
                    iz = itbl_macro_coord(3, imacro)
                    Ac_macro(1:3, imacro) = fw%vec_Ac_new%v(1:3, ix, iy, iz)
                end do
            end if
            call comm_bcast(Ac_macro, icomm, 0)
            
            Jmat_macro_tmp = 0.0d0
            do imacro = imacro_min, imacro_max
                call dt_evolve_bloch(sbe(imacro), gs, Ac_macro(1:3, imacro), dt)
                call calc_current_bloch(sbe(imacro), gs, Ac_macro(1:3, imacro), jmat, icomm_macro)

                if (mod(it, 10) == 0) then
                    if (irank_macro == 0) then
                        write(1000+imacro, '(f12.6,15(es24.15e3))') t, Ac_macro(1:3, imacro), fw%vec_e%v(1:3, ix, iy, iz), &
                            & Ac_macro(1:3, imacro), fw%vec_e%v(1:3, ix, iy, iz), jmat(1:3)
                    end if    
                end if

                if (irank_macro == 0) then
                    Jmat_macro_tmp(1:3, imacro) = jmat(1:3)
                end if
            end do
            call comm_summation(Jmat_macro_tmp, Jmat_macro, 3*nmacro, icomm)

            do imacro = 1, nmacro
                ix = itbl_macro_coord(1, imacro)
                iy = itbl_macro_coord(2, imacro)
                iz = itbl_macro_coord(3, imacro)
                fw%vec_j_em%v(:, ix, iy, iz) = -Jmat_macro(:, imacro)
            end do
        end if

        call weyl_calc(fs, fw)

        if (irank == 0) then
            if (mod(it, out_ms_step) == 0) then
                call write_Ac_field(it, fs, fw)
            end if
            if (mod(it, 10) == 0) then
                write(*, "(a,i6)") "Time step = ", it
                if (flag_1d_model) call write_wave_data_file(it, fs, fw)
                call write_obs_data_file(it, fs, fw)
            end if
            if (mod(it, 1000) == 0) then
                if (flag_1d_model) flush(888)
                do iobs = 1, obs_num_em
                    flush(100+iobs)
                end do
            end if
        end if


    end do

    if (nmacro > 0) then
        if (irank_macro == 0) then
            if (flag_1d_model) close(99)
            do imacro = imacro_min, imacro_max
                close(1000+imacro)
            end do
            do iobs = 1, obs_num_em
                close(100+iobs)
            end do
            close(888)
        end if
    end if

    call comm_sync_all(icomm)

    return
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


! subroutine read_shape_info(fs, fw)
!     use input_parameter, only: file_ms_shape, epsilon_em, nx_m, ny_m, nz_m
!     use fdtd_weyl, only: ls_fdtd_weyl
!     implicit none
!     integer, intent(in) :: nmacro_max
!     integer, intent(out) :: itbl_macro_coord(1:3, nmacro_max)
!     integer, intent(out) :: nmacro
!     type(s_fdtd_system), intent(in) :: fs
!     type(ls_fdtd_weyl), intent(inout) :: fw

!     integer :: itbl_media( &
!         & fw%mg%is(1):fw%mg%ie(1), &
!         & fw%mg%is(2):fw%mg%ie(2), &
!         & fw%mg%is(3):fw%mg%ie(3), &
!         & )

!     do i = 1, n_s
!         select case trim(typ_s(i))
!         case "ellipsoid"
!         end select
!     end do
! end subroutine read_shape_info
        

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

    write(file_Ac_data, "(a,a,a,a,i6.6,a)") trim(sysname), "_sbe_RT_Ac/", trim(sysname), "_Ac_", iit, ".data"
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

 ! Experimetal Implementation of Incident/Reflection/Transmit field output
subroutine write_wave_data_file(iit, fs, fw)
    use input_parameter, only: sysname, nx_m
    use fdtd_weyl, only: s_fdtd_system, ls_fdtd_weyl
    use phys_constants, only: cspeed_au
    implicit none
    integer, intent(in) :: iit
    type(s_fdtd_system), intent(in) :: fs
    type(ls_fdtd_weyl), intent(in) :: fw
    real(8) :: dt

    real(8) :: e_inc(3)
    real(8) :: e_ref(3)
    real(8) :: e_tra(3)
    real(8) :: dt_Ac(3)
    real(8) :: dx_Ac(3)
    integer :: iiy, iiz

    iiy = fs%mg%is(2)
    iiz = fs%mg%is(3)
    dt = fw%dt

    ! Left side boundary:
    dx_Ac(:) = (fw%vec_Ac%v(:,0,iiy,iiz) - fw%vec_Ac%v(:,-1,iiy,iiz)) / fs%hgs(1)
    dt_Ac(:) = (0.5d0 * (fw%vec_Ac_new%v(:,0,iiy,iiz) + fw%vec_Ac_new%v(:,-1,iiy,iiz)) & 
        & - 0.5d0 * (fw%vec_Ac_old%v(:,0,iiy,iiz) + fw%vec_Ac_old%v(:,-1,iiy,iiz))) / (2 * dt)
    
    e_inc(:) = -0.5d0 * (dt_Ac - cspeed_au * dx_Ac)
    e_ref(:) = -0.5d0 * (dt_Ac + cspeed_au * dx_Ac)

    ! Right side boundary:
    dx_Ac(:) = (fw%vec_Ac%v(:,nx_m+2,iiy,iiz) - fw%vec_Ac%v(:,nx_m+1,iiy,iiz)) / fs%hgs(1)
    dt_Ac(:) = (0.5d0 * (fw%vec_Ac_new%v(:,nx_m+2,iiy,iiz) + fw%vec_Ac_new%v(:,nx_m+1,iiy,iiz)) & 
        & - 0.5d0 * (fw%vec_Ac_old%v(:,nx_m+2,iiy,iiz) + fw%vec_Ac_old%v(:,nx_m+1,iiy,iiz))) / (2 * dt)
    
    e_tra(:) = -0.5d0 * (dt_Ac - cspeed_au * dx_Ac)

    write(888, '(99(e23.15e3, 1x))')  &
        & iit * dt * 1.0d0, &
        & e_inc(1) * 1.0d0, &
        & e_inc(2) * 1.0d0, &
        & e_inc(3) * 1.0d0, &
        & e_ref(1) * 1.0d0, &
        & e_ref(2) * 1.0d0, &
        & e_ref(3) * 1.0d0, &
        & e_tra(1) * 1.0d0, &
        & e_tra(2) * 1.0d0, &
        & e_tra(3) * 1.0d0
    return
end subroutine write_wave_data_file

subroutine write_obs_data_file(iit, fs, fw)
    use input_parameter, only: sysname, nx_m, obs_num_em, obs_loc_em
    use fdtd_weyl, only: s_fdtd_system, ls_fdtd_weyl
    use phys_constants, only: cspeed_au
    implicit none
    integer, intent(in) :: iit
    type(s_fdtd_system), intent(in) :: fs
    type(ls_fdtd_weyl), intent(in) :: fw
    real(8) :: dt
    real(8) :: e(3)
    integer :: iix, iiy, iiz, iobs

    dt = fw%dt
    do iobs = 1, obs_num_em
        iix = int(obs_loc_em(iobs, 1) / fs%hgs(1))
        iiy = int(obs_loc_em(iobs, 2) / fs%hgs(2))
        iiz = int(obs_loc_em(iobs, 3) / fs%hgs(3))
        iix = min(fs%mg%ie(1), max(fs%mg%is(1), iix))
        iiy = min(fs%mg%ie(2), max(fs%mg%is(2), iiy))
        iiz = min(fs%mg%ie(3), max(fs%mg%is(3), iiz))
        e(:) =  -(fw%vec_Ac_new%v(:, iix, iiy, iiz) - fw%vec_Ac_old%v(:, iix, iiy, iiz)) / (2 * dt)
        write(100+iobs, "(f12.6,99es25.15e4)") iit * dt, e(1:3)
    end do
    return 
end subroutine write_obs_data_file

end module multiscale
