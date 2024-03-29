#!/usr/bin/env python

setting = [
    ("calculation", [
        ("theory", "character(256)", None, "'perturb_dielec'"),
    ]),
    ("control", [
        ("sysname", "character(256)", None, "'test'"),
        ("base_directory", "character(256)", None, "'./'"),
        ("gs_directory", "character(256)", None, "'./'"),
        ("read_sbe_gs_bin", "character(256)", None, "'n'"),
        ("write_sbe_gs_bin", "character(256)", None, "'y'"),
    ]),
    ("system", [
        ("al", "real(8)", 3, "(/0.0, 0.0, 0.0/)"),
        ("al_vec1", "real(8)", 3, "(/0.0, 0.0, 0.0/)"),
        ("al_vec2", "real(8)", 3, "(/0.0, 0.0, 0.0/)"),
        ("al_vec3", "real(8)", 3, "(/0.0, 0.0, 0.0/)"),
        ("nstate", "integer", None, "0"),
        ("nelec", "integer", None, "0"),
        ("nstate_sbe", "integer", None, "0"),
    ]),
    ("kgrid", [
        ("num_kgrid", "integer", 3, "(/0, 0, 0/)"),
    ]),
    ("tgrid", [
        ("nt", "integer", None, "1000"),
        ("dt", "real(8)", None, "1.0d-2"),
    ]),
    ("emfield", [
        ("e_impulse", "real(8)", None, "0.0d0"),
        ("ae_shape1", "character(256)", None, "'none'"),
        ("ae_shape2", "character(256)", None, "'none'"),
        ("epdir_re1", "real(8)", 3, "(/0.0, 0.0, 0.0/)"),
        ("epdir_re2", "real(8)", 3, "(/0.0, 0.0, 0.0/)"),
        ("epdir_im1", "real(8)", 3, "(/0.0, 0.0, 0.0/)"),
        ("epdir_im2", "real(8)", 3, "(/0.0, 0.0, 0.0/)"),
        ("phi_cep1", "real(8)", None, "0.0d0"),
        ("phi_cep2", "real(8)", None, "0.0d0"),
        ("E_amplitude1", "real(8)", None, "0.0d0"),
        ("E_amplitude2", "real(8)", None, "0.0d0"),
        ("I_wcm2_1", "real(8)", None, "0.0d0"),
        ("I_wcm2_2", "real(8)", None, "0.0d0"),
        ("tw1", "real(8)", None, "0.0d0"),
        ("tw2", "real(8)", None, "0.0d0"),
        ("omega1", "real(8)", None, "0.0d0"),
        ("omega2", "real(8)", None, "0.0d0"),
        ("t1_t2", "real(8)", None, "0.0d0"),
        ("t1_start", "real(8)", None, "0.0d0"),
    ]),
    ("analysis", [
        ("nenergy", "integer", None, "1000"),
        ("de", "real(8)", None, "1.0d-3"),
        ("gamma", "real(8)", None, "5.0d-3"),
        ("out_ms_step", "integer", None, "100"),
        ("out_ms_ix", "integer", 2, "(/-1000000, 1000000/)"),
        ("out_ms_iy", "integer", 2, "(/-1000000, 1000000/)"),
        ("out_ms_iz", "integer", 2, "(/-1000000, 1000000/)"),
        ("out_ms_it", "integer", 2, "(/-1000000, 1000000/)"),
    ]),
    ("multiscale", [
        ("fdtddim", "character(256)", None, "''"),
        ("twod_shape", "character(256)", None, "''"),
        ("nx_m", "integer", None, "0"),
        ("ny_m", "integer", None, "0"),
        ("nz_m", "integer", None, "0"),
        ("hx_m", "real(8)", None, "0.0d0"),
        ("hy_m", "real(8)", None, "0.0d0"),
        ("hz_m", "real(8)", None, "0.0d0"),
        ("nxvac_m", "integer", 2, "0"),
        ("nyvac_m", "integer", 2, "0"),
        ("nzvac_m", "integer", 2, "0"),
        ("file_ms_shape", "character(256)", None, "''"),
    ]),
    ("maxwell", [
        # ("al_em", "real(8)", 3, "0"),
        # ("num_rgrid_em", "integer", 3, "0"),
        # ("at_em", "real(8)", None, "0.0d0"),
        #("media_num", "integer", 9, "0"),
        #("media_type", "character(256)", 9, "'vacuum'"),
        ("epsilon_em", "real(8)", 9, "1.0d0"),
        # ("sigma_em", "real(8)", 9, "0.0d0"),
        # ("pole_num_ld", "integer", 9, "1"),
        # ("omega_p_ld", "real(8)", 9, "1.0d0"),
        # ("gamma_ld", "real(8)", 9, "1.0d0"),
        # ("omega_ld", "real(8)", 9, "1.0d0"),
        # ("ek_dir1", "real(8)", 3, "(/1.0,0.0,0.0/)"),
        # ("obs_num_em", "integer", None, "0"),
        # ("obs_loc_em", "real(8)", (3, 9), "(/1.0,0.0,0.0/)"),
        #("yn_make_shape", "character(1)", 9, "'n'"),
        #("yn_output_shape", "character(1)", 9, "'n'"),
        #("n_s", "integer", None, "0"),
        #("id_s", "integer", 9, "0"),
        #("typ_s", "character(256)", 9, "''"),
        #("inf_s", "real(8)", [9, 3], "0.0d0"),
        #("ori_s", "real(8)", [9, 3], "0.0d0"),
        ("obs_num_em", "integer", None, "0"),
        ("obs_loc_em", "real(8)", [9, 3], "0.0d0"),
    ]),
]






import os


template = r"""! This file is automatically created by {PROGRAM}
module input_parameter
    implicit none

{CODE_DEFINE}

contains

    subroutine read_input(icomm)
        use communication
        implicit none
        integer, intent(in) :: icomm
        integer :: ret, irank, nproc
        character(256) :: tmp

{CODE_NAMELIST}

{CODE_DEFAULT}

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
{CODE_REWIND}
            close(99)
{CODE_DUMP}
        end if
{CODE_BCAST}

    end subroutine read_input
end module input_parameter
"""

tbl_fmt = {
    "real(8)": "es25.15e3",
    "integer": "i9",
}

tbl_mpi = {
    "real(8)": "MPI_DOUBLE_PRECISION",
    "integer": "MPI_INTEGER",
}

ind1 = " "*4
ind2 = " "*8
ind3 = " "*12

code_define = ""
for group_name, group_data in setting:
    for var_name, var_type, var_dim, var_defval in group_data:
        if var_dim:
            if type(var_dim) is int:
                var_dim = str(var_dim)
            else:
                var_dim = ",".join([str(x) for x in var_dim])
            code_define += ind1 + "%s :: %s(%s)\n" % (var_type, var_name,var_dim)
        else:
            code_define += ind1 + "%s :: %s\n" % (var_type, var_name)

code_namelist = ""
for group_name, group_data in setting:
    code_namelist += ind2 + "namelist/%s/ &\n" % (group_name)
    for i, (var_name, var_type, var_dim, var_defval) in enumerate(group_data):
        if i > 0:
            code_namelist += ", &\n" 
        code_namelist += ind2 + "& " + var_name
    code_namelist += "\n"

code_default = ""
for group_name, group_data in setting:
    for var_name, var_type, var_dim, var_defval in group_data:
        code_default += ind2 + "%s = %s\n" % (var_name, var_defval)

code_rewind = ""
for group_name, group_data in setting:
    code_rewind += ind3 + "rewind(99); read(99, nml=%s, iostat=ret)\n" % group_name

code_dump = ""
for group_name, group_data in setting:
    for var_name, var_type, var_dim, var_defval in group_data:
        if var_type.startswith("character"):
            if var_dim:
                if type(var_dim) is int:
                    for i in range(1, var_dim+1):
                        code_dump += ind3 + "write(*,'(a,a)') '# %s: %s = ', trim(%s(%d))\n" % (group_name, var_name, var_name, i)
            else:
                code_dump += ind3 + "write(*,'(a,a)') '# %s: %s = ', trim(%s)\n" % (group_name, var_name, var_name)
        else:
            code_dump += ind3 + "write(*,'(a,99%s)') '# %s: %s = ', %s\n" % (tbl_fmt[var_type], group_name, var_name, var_name)

code_bcast = ""
for group_name, group_data in setting:
    for var_name, var_type, var_dim, var_defval in group_data:
        code_bcast += ind2 + "call comm_bcast(%s, icomm, 0)\n" % (var_name)

program = os.path.split(__file__)[-1]
f90file = os.path.splitext(__file__)[0] + ".f90"

with open(f90file, "w") as fh:
    fh.write(template.format(
        PROGRAM=program,
        CODE_DEFINE=code_define,
        CODE_NAMELIST=code_namelist,
        CODE_DEFAULT=code_default,
        CODE_REWIND=code_rewind,
        CODE_DUMP=code_dump,
        CODE_BCAST=code_bcast,
    ))
    print("Exported: %s" % fh.name)
