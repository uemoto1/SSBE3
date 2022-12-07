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
    integer :: num_kgrid(3)
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
    real(8) :: hx_m
    real(8) :: hy_m
    real(8) :: hz_m
    integer :: nxvac_m(2)
    integer :: nyvac_m(2)
    integer :: nzvac_m(2)
    character(256) :: file_ms_shape
    integer :: out_ms_ix(2)
    integer :: out_ms_iy(2)
    integer :: out_ms_iz(2)
    integer :: out_ms_it(2)
    integer :: media_num
    real(8) :: epsilon_em(9)


contains

    subroutine read_input(icomm)
        use communication
        implicit none
        integer, intent(in) :: icomm
        integer :: ret, irank, nproc, ierr
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
        & num_kgrid
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
        & hx_m, &
        & hy_m, &
        & hz_m, &
        & nxvac_m, &
        & nyvac_m, &
        & nzvac_m, &
        & file_ms_shape, &
        & out_ms_ix, &
        & out_ms_iy, &
        & out_ms_iz, &
        & out_ms_it
        namelist/maxwell/ &
        & media_num, &
        & epsilon_em


        theory = 'perturb_dielec'
        sysname = 'test'
        base_directory = './'
        gs_directory = './'
        read_sbe_gs_bin = 'n'
        write_sbe_gs_bin = 'y'
        al = (/0.0, 0.0, 0.0/)
        al_vec1 = (/0.0, 0.0, 0.0/)
        al_vec2 = (/0.0, 0.0, 0.0/)
        al_vec3 = (/0.0, 0.0, 0.0/)
        nstate = 0
        nelec = 0
        nstate_sbe = 0
        num_kgrid = (/0, 0, 0/)
        nt = 1000
        dt = 1.0d-2
        e_impulse = 0.0d0
        ae_shape1 = 'none'
        ae_shape2 = 'none'
        epdir_re1 = (/0.0, 0.0, 0.0/)
        epdir_re2 = (/0.0, 0.0, 0.0/)
        epdir_im1 = (/0.0, 0.0, 0.0/)
        epdir_im2 = (/0.0, 0.0, 0.0/)
        phi_cep1 = 0.0d0
        phi_cep2 = 0.0d0
        E_amplitude1 = 0.0d0
        E_amplitude2 = 0.0d0
        I_wcm2_1 = 0.0d0
        I_wcm2_2 = 0.0d0
        tw1 = 0.0d0
        tw2 = 0.0d0
        omega1 = 0.0d0
        omega2 = 0.0d0
        t1_t2 = 0.0d0
        t1_start = 0.0d0
        nenergy = 1000
        de = 1.0d-3
        gamma = 5.0d-3
        fdtddim = ''
        twod_shape = ''
        nx_m = 0
        ny_m = 0
        nz_m = 0
        hx_m = 0.0d0
        hy_m = 0.0d0
        hz_m = 0.0d0
        nxvac_m = 0
        nyvac_m = 0
        nzvac_m = 0
        file_ms_shape = ''
        out_ms_ix = (/-1000000, 1000000/)
        out_ms_iy = (/-1000000, 1000000/)
        out_ms_iz = (/-1000000, 1000000/)
        out_ms_it = (/-1000000, 1000000/)
        media_num = 0
        epsilon_em = 1.0d0


        call comm_get_groupinfo(icomm, irank, nproc)

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
            rewind(99); read(99, nml=maxwell, iostat=ret)

            close(99)
            write(*,'(a,a)') '# calculation: theory = ', trim(theory)
            write(*,'(a,a)') '# control: sysname = ', trim(sysname)
            write(*,'(a,a)') '# control: base_directory = ', trim(base_directory)
            write(*,'(a,a)') '# control: gs_directory = ', trim(gs_directory)
            write(*,'(a,a)') '# control: read_sbe_gs_bin = ', trim(read_sbe_gs_bin)
            write(*,'(a,a)') '# control: write_sbe_gs_bin = ', trim(write_sbe_gs_bin)
            write(*,'(a,99es25.15e3)') '# system: al = ', al
            write(*,'(a,99es25.15e3)') '# system: al_vec1 = ', al_vec1
            write(*,'(a,99es25.15e3)') '# system: al_vec2 = ', al_vec2
            write(*,'(a,99es25.15e3)') '# system: al_vec3 = ', al_vec3
            write(*,'(a,99i9)') '# system: nstate = ', nstate
            write(*,'(a,99i9)') '# system: nelec = ', nelec
            write(*,'(a,99i9)') '# system: nstate_sbe = ', nstate_sbe
            write(*,'(a,99i9)') '# kgrid: num_kgrid = ', num_kgrid
            write(*,'(a,99i9)') '# tgrid: nt = ', nt
            write(*,'(a,99es25.15e3)') '# tgrid: dt = ', dt
            write(*,'(a,99es25.15e3)') '# emfield: e_impulse = ', e_impulse
            write(*,'(a,a)') '# emfield: ae_shape1 = ', trim(ae_shape1)
            write(*,'(a,a)') '# emfield: ae_shape2 = ', trim(ae_shape2)
            write(*,'(a,99es25.15e3)') '# emfield: epdir_re1 = ', epdir_re1
            write(*,'(a,99es25.15e3)') '# emfield: epdir_re2 = ', epdir_re2
            write(*,'(a,99es25.15e3)') '# emfield: epdir_im1 = ', epdir_im1
            write(*,'(a,99es25.15e3)') '# emfield: epdir_im2 = ', epdir_im2
            write(*,'(a,99es25.15e3)') '# emfield: phi_cep1 = ', phi_cep1
            write(*,'(a,99es25.15e3)') '# emfield: phi_cep2 = ', phi_cep2
            write(*,'(a,99es25.15e3)') '# emfield: E_amplitude1 = ', E_amplitude1
            write(*,'(a,99es25.15e3)') '# emfield: E_amplitude2 = ', E_amplitude2
            write(*,'(a,99es25.15e3)') '# emfield: I_wcm2_1 = ', I_wcm2_1
            write(*,'(a,99es25.15e3)') '# emfield: I_wcm2_2 = ', I_wcm2_2
            write(*,'(a,99es25.15e3)') '# emfield: tw1 = ', tw1
            write(*,'(a,99es25.15e3)') '# emfield: tw2 = ', tw2
            write(*,'(a,99es25.15e3)') '# emfield: omega1 = ', omega1
            write(*,'(a,99es25.15e3)') '# emfield: omega2 = ', omega2
            write(*,'(a,99es25.15e3)') '# emfield: t1_t2 = ', t1_t2
            write(*,'(a,99es25.15e3)') '# emfield: t1_start = ', t1_start
            write(*,'(a,99i9)') '# analysis: nenergy = ', nenergy
            write(*,'(a,99es25.15e3)') '# analysis: de = ', de
            write(*,'(a,99es25.15e3)') '# analysis: gamma = ', gamma
            write(*,'(a,a)') '# multiscale: fdtddim = ', trim(fdtddim)
            write(*,'(a,a)') '# multiscale: twod_shape = ', trim(twod_shape)
            write(*,'(a,99i9)') '# multiscale: nx_m = ', nx_m
            write(*,'(a,99i9)') '# multiscale: ny_m = ', ny_m
            write(*,'(a,99i9)') '# multiscale: nz_m = ', nz_m
            write(*,'(a,99es25.15e3)') '# multiscale: hx_m = ', hx_m
            write(*,'(a,99es25.15e3)') '# multiscale: hy_m = ', hy_m
            write(*,'(a,99es25.15e3)') '# multiscale: hz_m = ', hz_m
            write(*,'(a,99i9)') '# multiscale: nxvac_m = ', nxvac_m
            write(*,'(a,99i9)') '# multiscale: nyvac_m = ', nyvac_m
            write(*,'(a,99i9)') '# multiscale: nzvac_m = ', nzvac_m
            write(*,'(a,a)') '# multiscale: file_ms_shape = ', trim(file_ms_shape)
            write(*,'(a,99i9)') '# multiscale: out_ms_ix = ', out_ms_ix
            write(*,'(a,99i9)') '# multiscale: out_ms_iy = ', out_ms_iy
            write(*,'(a,99i9)') '# multiscale: out_ms_iz = ', out_ms_iz
            write(*,'(a,99i9)') '# multiscale: out_ms_it = ', out_ms_it
            write(*,'(a,99i9)') '# maxwell: media_num = ', media_num
            write(*,'(a,99es25.15e3)') '# maxwell: epsilon_em = ', epsilon_em

        end if
        call comm_bcast(theory, icomm, 0)
        call comm_bcast(sysname, icomm, 0)
        call comm_bcast(base_directory, icomm, 0)
        call comm_bcast(gs_directory, icomm, 0)
        call comm_bcast(read_sbe_gs_bin, icomm, 0)
        call comm_bcast(write_sbe_gs_bin, icomm, 0)
        call comm_bcast(al, icomm, 0)
        call comm_bcast(al_vec1, icomm, 0)
        call comm_bcast(al_vec2, icomm, 0)
        call comm_bcast(al_vec3, icomm, 0)
        call comm_bcast(nstate, icomm, 0)
        call comm_bcast(nelec, icomm, 0)
        call comm_bcast(nstate_sbe, icomm, 0)
        call comm_bcast(num_kgrid, icomm, 0)
        call comm_bcast(nt, icomm, 0)
        call comm_bcast(dt, icomm, 0)
        call comm_bcast(e_impulse, icomm, 0)
        call comm_bcast(ae_shape1, icomm, 0)
        call comm_bcast(ae_shape2, icomm, 0)
        call comm_bcast(epdir_re1, icomm, 0)
        call comm_bcast(epdir_re2, icomm, 0)
        call comm_bcast(epdir_im1, icomm, 0)
        call comm_bcast(epdir_im2, icomm, 0)
        call comm_bcast(phi_cep1, icomm, 0)
        call comm_bcast(phi_cep2, icomm, 0)
        call comm_bcast(E_amplitude1, icomm, 0)
        call comm_bcast(E_amplitude2, icomm, 0)
        call comm_bcast(I_wcm2_1, icomm, 0)
        call comm_bcast(I_wcm2_2, icomm, 0)
        call comm_bcast(tw1, icomm, 0)
        call comm_bcast(tw2, icomm, 0)
        call comm_bcast(omega1, icomm, 0)
        call comm_bcast(omega2, icomm, 0)
        call comm_bcast(t1_t2, icomm, 0)
        call comm_bcast(t1_start, icomm, 0)
        call comm_bcast(nenergy, icomm, 0)
        call comm_bcast(de, icomm, 0)
        call comm_bcast(gamma, icomm, 0)
        call comm_bcast(fdtddim, icomm, 0)
        call comm_bcast(twod_shape, icomm, 0)
        call comm_bcast(nx_m, icomm, 0)
        call comm_bcast(ny_m, icomm, 0)
        call comm_bcast(nz_m, icomm, 0)
        call comm_bcast(hx_m, icomm, 0)
        call comm_bcast(hy_m, icomm, 0)
        call comm_bcast(hz_m, icomm, 0)
        call comm_bcast(nxvac_m, icomm, 0)
        call comm_bcast(nyvac_m, icomm, 0)
        call comm_bcast(nzvac_m, icomm, 0)
        call comm_bcast(file_ms_shape, icomm, 0)
        call comm_bcast(out_ms_ix, icomm, 0)
        call comm_bcast(out_ms_iy, icomm, 0)
        call comm_bcast(out_ms_iz, icomm, 0)
        call comm_bcast(out_ms_it, icomm, 0)
        call comm_bcast(media_num, icomm, 0)
        call comm_bcast(epsilon_em, icomm, 0)


    end subroutine read_input
end module input_parameter
