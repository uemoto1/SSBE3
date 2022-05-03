! This file is automatically created by input_parameter.py
module input_parameter
    implicit none

    character(256) :: theory
    character(256) :: sysname
    character(256) :: base_directory
    character(256) :: gs_directory
    character(256) :: read_sbe_gs_bin
    character(256) :: write_sbe_gs_bin
    real(8) :: al(3)
    real(8) :: al_vec1(3)
    real(8) :: al_vec2(3)
    real(8) :: al_vec3(3)
    integer :: nstate
    integer :: nelec
    integer :: nstate_sbe
    integer :: nkgrid(3)
    integer :: nt
    real(8) :: dt
    real(8) :: e_impulse
    character(256) :: ae_shape1
    character(256) :: ae_shape2
    real(8) :: epdir_re1(3)
    real(8) :: epdir_re2(3)
    real(8) :: epdir_im1(3)
    real(8) :: epdir_im2(3)
    real(8) :: phi_cep1
    real(8) :: phi_cep2
    real(8) :: E_amplitude1
    real(8) :: E_amplitude2
    real(8) :: I_wcm2_1
    real(8) :: I_wcm2_2
    real(8) :: tw1
    real(8) :: tw2
    real(8) :: omega1
    real(8) :: omega2
    real(8) :: t1_t2
    real(8) :: t1_start
    integer :: nenergy
    real(8) :: de
    real(8) :: gamma
    character(256) :: fdtddim
    character(256) :: twod_shape
    integer :: nx_m
    integer :: ny_m
    integer :: nz_m
    integer :: nmacro
    real(8) :: hx_m
    real(8) :: hy_m
    real(8) :: hz_m
    integer :: nxvac_m
    integer :: nyvac_m
    integer :: nzvac_m
    character(256) :: file_macropoint


contains

    subroutine read_input()
        use mpi
        implicit none
        integer :: ret, irank, ierr
        character(256) :: tmp

        namelist/calculation/ &
        & theory
        namelist/control/ &
        & sysname, &
        & base_directory, &
        & gs_directory, &
        & read_sbe_gs_bin, &
        & write_sbe_gs_bin
        namelist/system/ &
        & al, &
        & al_vec1, &
        & al_vec2, &
        & al_vec3, &
        & nstate, &
        & nelec, &
        & nstate_sbe
        namelist/kgrid/ &
        & nkgrid
        namelist/tgrid/ &
        & nt, &
        & dt
        namelist/emfield/ &
        & e_impulse, &
        & ae_shape1, &
        & ae_shape2, &
        & epdir_re1, &
        & epdir_re2, &
        & epdir_im1, &
        & epdir_im2, &
        & phi_cep1, &
        & phi_cep2, &
        & E_amplitude1, &
        & E_amplitude2, &
        & I_wcm2_1, &
        & I_wcm2_2, &
        & tw1, &
        & tw2, &
        & omega1, &
        & omega2, &
        & t1_t2, &
        & t1_start
        namelist/analysis/ &
        & nenergy, &
        & de, &
        & gamma
        namelist/multiscale/ &
        & fdtddim, &
        & twod_shape, &
        & nx_m, &
        & ny_m, &
        & nz_m, &
        & nmacro, &
        & hx_m, &
        & hy_m, &
        & hz_m, &
        & nxvac_m, &
        & nyvac_m, &
        & nzvac_m, &
        & file_macropoint


        file_macropoint = ''


        call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)

        if (irank == 0) then
            open(99, file='.namelist.tmp', action='write')
            do while (.true.)
                read(*, '(a)', iostat=ret) tmp
                if (ret < 0) exit ! End of file
                write(99, '(a)') trim(tmp)
            end do
            close(99)
            open(99, file='.namelist.tmp', action='read')
            rewind(99); read(99, nml=calculation, iostat=ret)
            rewind(99); read(99, nml=control, iostat=ret)
            rewind(99); read(99, nml=system, iostat=ret)
            rewind(99); read(99, nml=kgrid, iostat=ret)
            rewind(99); read(99, nml=tgrid, iostat=ret)
            rewind(99); read(99, nml=emfield, iostat=ret)
            rewind(99); read(99, nml=analysis, iostat=ret)
            rewind(99); read(99, nml=multiscale, iostat=ret)

            close(99)
        end if
        write(*,'(a,a)') 'calculation: theory = ', trim(theory)
        write(*,'(a,a)') 'control: sysname = ', trim(sysname)
        write(*,'(a,a)') 'control: base_directory = ', trim(base_directory)
        write(*,'(a,a)') 'control: gs_directory = ', trim(gs_directory)
        write(*,'(a,a)') 'control: read_sbe_gs_bin = ', trim(read_sbe_gs_bin)
        write(*,'(a,a)') 'control: write_sbe_gs_bin = ', trim(write_sbe_gs_bin)
        write(*,'(a,99es25.15e3)') 'system: al = ', al
        write(*,'(a,99es25.15e3)') 'system: al_vec1 = ', al_vec1
        write(*,'(a,99es25.15e3)') 'system: al_vec2 = ', al_vec2
        write(*,'(a,99es25.15e3)') 'system: al_vec3 = ', al_vec3
        write(*,'(a,99i9)') 'system: nstate = ', nstate
        write(*,'(a,99i9)') 'system: nelec = ', nelec
        write(*,'(a,99i9)') 'system: nstate_sbe = ', nstate_sbe
        write(*,'(a,99i9)') 'kgrid: nkgrid = ', nkgrid
        write(*,'(a,99i9)') 'tgrid: nt = ', nt
        write(*,'(a,99es25.15e3)') 'tgrid: dt = ', dt
        write(*,'(a,99es25.15e3)') 'emfield: e_impulse = ', e_impulse
        write(*,'(a,a)') 'emfield: ae_shape1 = ', trim(ae_shape1)
        write(*,'(a,a)') 'emfield: ae_shape2 = ', trim(ae_shape2)
        write(*,'(a,99es25.15e3)') 'emfield: epdir_re1 = ', epdir_re1
        write(*,'(a,99es25.15e3)') 'emfield: epdir_re2 = ', epdir_re2
        write(*,'(a,99es25.15e3)') 'emfield: epdir_im1 = ', epdir_im1
        write(*,'(a,99es25.15e3)') 'emfield: epdir_im2 = ', epdir_im2
        write(*,'(a,99es25.15e3)') 'emfield: phi_cep1 = ', phi_cep1
        write(*,'(a,99es25.15e3)') 'emfield: phi_cep2 = ', phi_cep2
        write(*,'(a,99es25.15e3)') 'emfield: E_amplitude1 = ', E_amplitude1
        write(*,'(a,99es25.15e3)') 'emfield: E_amplitude2 = ', E_amplitude2
        write(*,'(a,99es25.15e3)') 'emfield: I_wcm2_1 = ', I_wcm2_1
        write(*,'(a,99es25.15e3)') 'emfield: I_wcm2_2 = ', I_wcm2_2
        write(*,'(a,99es25.15e3)') 'emfield: tw1 = ', tw1
        write(*,'(a,99es25.15e3)') 'emfield: tw2 = ', tw2
        write(*,'(a,99es25.15e3)') 'emfield: omega1 = ', omega1
        write(*,'(a,99es25.15e3)') 'emfield: omega2 = ', omega2
        write(*,'(a,99es25.15e3)') 'emfield: t1_t2 = ', t1_t2
        write(*,'(a,99es25.15e3)') 'emfield: t1_start = ', t1_start
        write(*,'(a,99i9)') 'analysis: nenergy = ', nenergy
        write(*,'(a,99es25.15e3)') 'analysis: de = ', de
        write(*,'(a,99es25.15e3)') 'analysis: gamma = ', gamma
        write(*,'(a,a)') 'multiscale: fdtddim = ', trim(fdtddim)
        write(*,'(a,a)') 'multiscale: twod_shape = ', trim(twod_shape)
        write(*,'(a,99i9)') 'multiscale: nx_m = ', nx_m
        write(*,'(a,99i9)') 'multiscale: ny_m = ', ny_m
        write(*,'(a,99i9)') 'multiscale: nz_m = ', nz_m
        write(*,'(a,99i9)') 'multiscale: nmacro = ', nmacro
        write(*,'(a,99es25.15e3)') 'multiscale: hx_m = ', hx_m
        write(*,'(a,99es25.15e3)') 'multiscale: hy_m = ', hy_m
        write(*,'(a,99es25.15e3)') 'multiscale: hz_m = ', hz_m
        write(*,'(a,99i9)') 'multiscale: nxvac_m = ', nxvac_m
        write(*,'(a,99i9)') 'multiscale: nyvac_m = ', nyvac_m
        write(*,'(a,99i9)') 'multiscale: nzvac_m = ', nzvac_m
        write(*,'(a,a)') 'multiscale: file_macropoint = ', trim(file_macropoint)

        call MPI_BCAST(theory, 256, MPI_CHARACTER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(sysname, 256, MPI_CHARACTER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(base_directory, 256, MPI_CHARACTER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(gs_directory, 256, MPI_CHARACTER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(read_sbe_gs_bin, 256, MPI_CHARACTER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(write_sbe_gs_bin, 256, MPI_CHARACTER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(al, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(al_vec1, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(al_vec2, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(al_vec3, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(nstate, 1, MPI_INTEGER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(nelec, 1, MPI_INTEGER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(nstate_sbe, 1, MPI_INTEGER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(nkgrid, 3, MPI_INTEGER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(nt, 1, MPI_INTEGER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(e_impulse, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(ae_shape1, 256, MPI_CHARACTER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(ae_shape2, 256, MPI_CHARACTER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(epdir_re1, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(epdir_re2, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(epdir_im1, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(epdir_im2, 3, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(phi_cep1, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(phi_cep2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(E_amplitude1, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(E_amplitude2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(I_wcm2_1, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(I_wcm2_2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(tw1, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(tw2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(omega1, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(omega2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(t1_t2, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(t1_start, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(nenergy, 1, MPI_INTEGER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(de, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(gamma, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(fdtddim, 256, MPI_CHARACTER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(twod_shape, 256, MPI_CHARACTER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(nx_m, 1, MPI_INTEGER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(ny_m, 1, MPI_INTEGER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(nz_m, 1, MPI_INTEGER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(nmacro, 1, MPI_INTEGER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(hx_m, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(hy_m, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(hz_m, 1, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(nxvac_m, 1, MPI_INTEGER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(nyvac_m, 1, MPI_INTEGER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(nzvac_m, 1, MPI_INTEGER, MPI_COMM_WORLD, 0, ierr)
        call MPI_BCAST(file_macropoint, 256, MPI_CHARACTER, MPI_COMM_WORLD, 0, ierr)


    end subroutine read_input
end module input_parameter

