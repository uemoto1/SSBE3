program main
    use mpi
    use omp_lib
    use sbe_gs
    use sbe_solver
    use input_parameter
    use test
    use em_field
    implicit none

    type(s_sbe_bloch_solver) :: sbe
    type(s_sbe_gs) :: gs
    real(8) :: t,  E(3), jmat(3)
    real(8), allocatable :: Ac_ext_t(:, :)
    integer :: it, i, ikey
    real(8) :: energy0, energy
    real(8) :: tr_all, tr_vb
    integer :: nproc, irank, ierr, icomm_macro
    integer :: nthread, nmacro_proc, nproc_macro
    integer, allocatable :: itbl_macro_min(:)
    integer, allocatable :: itbl_macro_max(:)
    integer :: imacro_min, imacro_max
    integer :: irank_macro, isize_macro

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)

    !$omp parallel
    !$omp master
    nthread = omp_get_num_threads()
    !$omp end master
    !$omp end parallel

    if (irank == 0) then
        write(*, "(a)") "# Parallelization:", nproc
        write(*, "(a,i6)") "# number of MPI processes=", nproc
        write(*, "(a,i6)") "# number of OMP threads=", nthread
    end if

    call read_input()

    ! Distribute
    allocate(itbl_macro_min(0:nproc-1))
    allocate(itbl_macro_max(0:nproc-1))

    if (nproc <= nmacro) then
        if (mod(nmacro, nproc) == 0) then
            nmacro_proc = nmacro / nproc
            do i = 0, nproc-1
                itbl_macro_min(i) = i * nmacro_proc + 1
                itbl_macro_max(i) = itbl_macro_min(i) + (nmacro_proc - 1)
            end do
        else
            call MPI_FINALIZE(ierr)
            stop "ERROR: mod(nmacro, nproc) != 0"
        end if
    else
        if (mod(nproc, nmacro) == 0) then
            nproc_macro = nproc / nmacro
            do i = 0, nproc-1
                itbl_macro_min(i) = (i / nproc_macro) + 1
                itbl_macro_max(i) = itbl_macro_min(i)
            end do
        else
            call MPI_FINALIZE(ierr)
            stop "ERROR: mod(nproc, nmacro) != 0"
        end if
    end if
    if (irank == 0) then
        do i = 0, nproc-1
            write(*, "(a,i9,a,2i9)") "# rank=", i, " imacro=", itbl_macro_min(i), itbl_macro_max(i)
        end do
    end if
    imacro_min = itbl_macro_min(irank)
    imacro_max = itbl_macro_max(irank)
    
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, imacro_min, 0, icomm_macro, ierr)
    ! call MPI_COMM_RANK(icomm_macro, irank_macro, ierr)
    ! call MPI_COMM_SIZE(icomm_macro, isize_macro, ierr)
    ! write(*, *) irank, icomm_macro, irank_macro, isize_macro, imacro_min
    ! stop "---"

    if (0.0d0 < al(1)) al_vec1(1:3) = (/ al(1), 0.0d0, 0.0d0 /)
    if (0.0d0 < al(2)) al_vec2(1:3) = (/ 0.0d0, al(2), 0.0d0 /)
    if (0.0d0 < al(3)) al_vec3(1:3) = (/ 0.0d0, 0.0d0, al(3) /)
    if (nstate_sbe < 1) nstate_sbe = nstate

    ! Read ground state electronic system:
    call init_sbe_gs(gs, sysname, gs_directory, &
        & nkgrid, nstate, nelec, &
        & al_vec1, al_vec2, al_vec3, &
        & .false., MPI_COMM_WORLD)        

    ! stop "--"


    ! Initialization of SBE solver and density matrix:
    call init_sbe(sbe, gs, nstate_sbe, icomm_macro)

    ! Prepare external pulse
    allocate(Ac_ext_t(1:3, -1:nt+1))
    call calc_Ac_ext_t(0.0d0, dt, 0, nt, Ac_ext_t)

    energy0 = calc_energy(sbe, gs, Ac_ext_t(:, 0), icomm_macro)

    ! Realtime calculation
    if (irank == 0) then
        open(unit=100, file=trim(base_directory)//trim(sysname)//"_sbe_rt.data")
        write(100, '(4a)') "# 1:Time[a.u.] 2:Ac_ext_x[a.u.] 3:Ac_ext_y[a.u.] 4:Ac_ext_z[a.u.] ", &
            & "5:E_ext_x[a.u.] 6:E_ext_y[a.u.] 7:E_ext_z[a.u.] 8:Ac_tot_x[a.u.] ", &
            & "9:Ac_tot_y[a.u.] 10:Ac_tot_z[a.u.] 11:E_tot_x[a.u.] 12:E_tot_y[a.u.] ", &
            & "13:E_tot_z[a.u.]  14:Jm_x[a.u.] 15:Jm_y[a.u.] 16:Jm_z[a.u.]"

        open(unit=101, file=trim(base_directory)//trim(sysname)//"_sbe_rt_energy.data")
        write(101, '(a)') "# 1:Time[a.u.] 2:Eall[a.u.] 3:Eall-Eall0[a.u.]"

        open(unit=102, file=trim(base_directory)//trim(sysname)//"_sbe_nex.data")
        write(102, '(a)') "# 1:Time[a.u.] 2:nelec[a.u.] 3:nhole[a.u.]"

        write(101, '(f12.6,2(es24.15e3))') 0.0d0, energy0, 0.0d0
    end if



    do it = 1, nt
        t = dt * it
        call dt_evolve_bloch(sbe, gs, Ac_ext_t(:, it), dt)

        if (mod(it, 10) == 0) then
            E(:) = (Ac_ext_t(:, it + 1) - Ac_ext_t(:, it - 1)) / (2 * dt)
            call calc_current_bloch(sbe, gs, Ac_ext_t(:, it), Jmat, icomm_macro)
            energy = calc_energy(sbe, gs, Ac_ext_t(:, it), icomm_macro)
            tr_all = calc_trace(sbe, gs, nstate_sbe, icomm_macro)
            tr_vb = calc_trace(sbe, gs, nelec / 2, icomm_macro)
            
            if (irank == 0) then
                write(100, '(f12.6,15(es24.15e3))') t, Ac_ext_t(:, it), E(:), Ac_ext_t(:, it), E(:), Jmat(:)
                write(101, '(f12.6,2(es24.15e3))') t, energy, energy-energy0
                write(102, '(f12.6,2(es24.15e3))') t, tr_all - tr_vb, nelec - tr_vb 
                write(*, '(f12.6,es24.15e3)') t, tr_all
            end if

            write(1000+irank, '(f12.6,15(es24.15e3))') t, Ac_ext_t(:, it), E(:), Ac_ext_t(:, it), E(:), Jmat(:)
        end if

    end do

    if (irank == 0) then
        close(100)
        close(101)
        close(102)
    end if


    call MPI_FINALIZE(ierr)
    stop "Bye!"

    
end program 
