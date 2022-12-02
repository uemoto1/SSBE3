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
    real(8) :: t,  E(3), jmat(3)
    real(8), allocatable :: Ac_ext_t(:, :)
    integer :: it, i, j, n
    real(8) :: energy0, energy
    real(8) :: tr_all, tr_vb
    integer :: nproc, irank, ierr, icomm_macro
    integer :: nthread, nmacro_proc, nproc_macro
    integer, allocatable :: itbl_macro_min(:)
    integer, allocatable :: itbl_macro_max(:)
    integer, allocatable :: itbl_proc_min(:)
    integer, allocatable :: itbl_proc_max(:)
    real(8), allocatable :: vec_Ac_macro(:, :)
    real(8), allocatable :: vec_Jem_macro(:, :)
    integer, parameter :: nmacro_max = 10000
    integer :: nmacro, imacro, itype
    integer :: itbl_macro_coord(3, nmacro_max)
    integer :: imacro_min, imacro_max
    integer :: irank_macro
    real(8) :: x, r_it, w, rtmp
    integer :: ix, iy, iz
    integer :: mt

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

    if (irank == 0) then
        call weyl_init(fs, fw)

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
    end if
    call comm_bcast(nmacro, icomm, 0)
    call comm_bcast(itbl_macro_coord, icomm, 0)
    allocate(vec_Ac_macro(1:3, nmacro), vec_Jem_macro(1:3, nmacro))
 
    if (nmacro > 0) then
        if (0.0d0 < al(1)) al_vec1(1:3) = (/ al(1), 0.0d0, 0.0d0 /)
        if (0.0d0 < al(2)) al_vec2(1:3) = (/ 0.0d0, al(2), 0.0d0 /)
        if (0.0d0 < al(3)) al_vec3(1:3) = (/ 0.0d0, 0.0d0, al(3) /)
        if (nstate_sbe < 1) nstate_sbe = nstate

        ! Distribute
        allocate(itbl_macro_min(0:nproc-1))
        allocate(itbl_macro_max(0:nproc-1))
        allocate(itbl_proc_min(1:nmacro))
        allocate(itbl_proc_max(1:nmacro))

        if (nproc <= nmacro) then
            call split_range(1, nmacro, nproc, itbl_macro_min, itbl_macro_max)
        else
            call split_range(0, nproc-1, nmacro, itbl_proc_min, itbl_proc_max)
            do imacro = 1, nmacro
                itbl_macro_min(itbl_proc_min(imacro):itbl_proc_max(imacro)) = imacro
                itbl_macro_max(itbl_proc_min(imacro):itbl_proc_max(imacro)) = imacro
            end do
        end if

        if (irank == 0) then
            do i = 0, nproc-1
                write(*, "(a,i9,a,2i9)") "# rank=", i, " imacro=", itbl_macro_min(i), itbl_macro_max(i)
            end do
        end if

        imacro_min = itbl_macro_min(irank)
        imacro_max = itbl_macro_max(irank)
        icomm_macro = comm_create_group(icomm, imacro_min, 0)
        ! call MPI_COMM_SPLIT(MPI_COMM_WORLD, imacro_min, 0, icomm_macro, ierr)
        call comm_get_groupinfo(icomm_macro ,irank_macro, nproc_macro)
        
        ! Read ground state electronic system:
        call init_sbe_gs(gs, sysname, gs_directory, &
            & nkgrid, nstate, nelec, &
            & al_vec1, al_vec2, al_vec3, &
            & .false., MPI_COMM_WORLD)
        ! Initialization of SBE solver and density matrix:
        allocate(sbe(imacro_min:imacro_max))
        do i = imacro_min, imacro_max
            call init_sbe(sbe(i), gs, nstate_sbe, icomm_macro)
        end do    
    end if

    ! Prepare external pulse
    mt = max(nt, int(abs(nxvac_m(1)) * fs%hgs(1) / cspeed_au / dt))
    allocate(Ac_ext_t(1:3, -1:mt+1))
    call calc_Ac_ext_t(0.0d0, dt, 0, mt, Ac_ext_t)    

    do ix = fs%mg%is_array(1), 0
        do iy = fs%mg%is_array(2), fs%mg%ie_array(2)
            do iz = fs%mg%is_array(3), fs%mg%ie_array(3)
                x = ix * fs%hgs(1)
                r_it = (-x / cspeed_au / fw%dt)
                it = int(r_it)
                if ((0 <= it) .and. (it < mt)) then
                    w = r_it - it
                    fw%vec_Ac_new%v(:, ix, iy, iz) = w * Ac_ext_t(:, it+1) + (1.0d0-w) * Ac_ext_t(:, it)
                end if
                r_it = r_it - 1
                it = int(r_it)
                if ((0 <= it) .and. (it < mt)) then
                    w = r_it - it
                    fw%vec_Ac%v(:, ix, iy, iz) = w * Ac_ext_t(:, it+1) + (1.0d0-w) * Ac_ext_t(:, it)
                end if
            end do
        end do
    end do


    if (nmacro > 0) then
        do i = imacro_min, imacro_max
            energy0 = calc_energy(sbe(i), gs, Ac_ext_t(:, 0), icomm_macro)
        end do

        if (irank_macro == 0) then
            do i = imacro_min, imacro_max
                write(tmp, "(a,a,a,i6.6,a)") trim(base_directory), trim(sysname), "_sbe_macro_",i,"_rt.data"
                open(1000+i, file=trim(tmp), action="write")
                write(1000+i, '(4a)') "# 1:Time[a.u.] 2:Ac_ext_x[a.u.] 3:Ac_ext_y[a.u.] 4:Ac_ext_z[a.u.] ", &
                    & "5:E_ext_x[a.u.] 6:E_ext_y[a.u.] 7:E_ext_z[a.u.] 8:Ac_tot_x[a.u.] ", &
                    & "9:Ac_tot_y[a.u.] 10:Ac_tot_z[a.u.] 11:E_tot_x[a.u.] 12:E_tot_y[a.u.] ", &
                    & "13:E_tot_z[a.u.]  14:Jm_x[a.u.] 15:Jm_y[a.u.] 16:Jm_z[a.u.]"
            end do
        end if
    end if

    do it = 1, nt
        t = dt * it

        if (nmacro > 0) then
            do i = imacro_min, imacro_max
                call dt_evolve_bloch(sbe(i), gs, Ac_ext_t(:, it), dt)
            end do

            if (mod(it, 10) == 0) then
                E(:) = (Ac_ext_t(:, it + 1) - Ac_ext_t(:, it - 1)) / (2 * dt)
                do i = imacro_min, imacro_max
                    call calc_current_bloch(sbe(i), gs, Ac_ext_t(:, it), Jmat, icomm_macro)
                    tr_all = calc_trace(sbe(i), gs, nstate_sbe, icomm_macro)
                    if (irank_macro == 0) then
                        write(1000+i, '(f12.6,15(es24.15e3))') t, Ac_ext_t(:, it), E(:), Ac_ext_t(:, it), E(:), Jmat(:)
                    end if    
                end do

                if (irank == 0) then
                    write(*, '(f12.6,es24.15e3)') t, tr_all
                end if
            end if
        end if


        if (irank == 0) then
            call weyl_calc(fs, fw)

            if ((mod(it, 100) == 0)) then
                write(tmp, "(a,i6.6,a)") "test", it, ".txt"
                open(999, file=trim(tmp), action="write")
                write(999, "(a)") "# 1:ix 2:iy 3:iz 4:x 5:y 6:z 7:Acx 8:Acy " & 
                // "9:Acz 10:Ex 11:Ey 12:Ez 13:Hx 14:Hy 15:Hz"
                do iz = max(fs%mg%is(3), out_ms_iz(1)), min(fs%mg%ie(3), out_ms_iz(2))
                do iy = max(fs%mg%is(2), out_ms_iy(1)), min(fs%mg%ie(2), out_ms_iy(2))
                do ix = max(fs%mg%is(1), out_ms_ix(1)), min(fs%mg%ie(1), out_ms_ix(2))
                    write(999, "(3i9,3f15.3,99es25.15e4)") &
                    ix, iy, iz, &
                    ix*fs%hgs(1), iy*fs%hgs(2), iz*fs%hgs(3), &
                    fw%vec_Ac%v(:, ix, iy, iz), &
                    fw%vec_e%v(:, ix, iy, iz), &
                    fw%vec_h%v(:, ix, iy, iz)
                end do
                end do
                end do
                close(999)
            end if
        end if
    end do

    if (nmacro > 0) then
        if (irank_macro == 0) then
            do i = imacro_min, imacro_max
                close(1000+i)
            end do
        end if
    end if

    return    
end subroutine multiscale_main
end module multiscale
