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
        ("nkgrid", "integer", 3, "(/0, 0, 0/)"),
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
    ]),
    ("multiscale", [
        ("fdtddim", "character(256)", None, "''"),
        ("twod_shape", "character(256)", None, "''"),
        ("nx_m", "integer", None, "0"),
        ("ny_m", "integer", None, "0"),
        ("nz_m", "integer", None, "0"),
        ("nmacro", "integer", None, "0"),
        ("hx_m", "real(8)", None, "0.0d0"),
        ("hy_m", "real(8)", None, "0.0d0"),
        ("hz_m", "real(8)", None, "0.0d0"),
        ("nxvac_m", "integer", None, "0"),
        ("nyvac_m", "integer", None, "0"),
        ("nzvac_m", "integer", None, "0"),
        ("file_macropoint", "character(256)", None, "''"),
    ]),
]









template = r"""! This file is automatically created by {PROGRAM}
module input_parameter
    implicit none

{CODE_DEFINE}

contains

    subroutine read_input()
        use mpi
        implicit none
        integer :: ret, irank, ierr
        character(256) :: tmp

{CODE_NAMELIST}

{CODE_DEFAULT}

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
            code_dump += ind3 + "write(*,'(a,a)') '%s: %s = ', trim(%s)\n" % (group_name, var_name, var_name)
        else:
            code_dump += ind3 + "write(*,'(a,99%s)') '%s: %s = ', %s\n" % (tbl_fmt[var_type], group_name, var_name, var_name)

code_bcast = ""
for group_name, group_data in setting:
    for var_name, var_type, var_dim, var_defval in group_data:
        if var_type.startswith("character"):
            mpi_len = var_type.split("(")[1].split(")")[0]
            mpi_type = "MPI_CHARACTER"
        else:
            if var_dim:
                mpi_len = var_dim
            else:
                mpi_len = 1
            mpi_type = tbl_mpi[var_type]
        code_bcast += ind2 + "call MPI_BCAST(%s, %s, %s, MPI_COMM_WORLD, 0, ierr)\n" % (var_name, mpi_len, mpi_type)

print(template.format(
    PROGRAM=__file__.split("/")[-1],
    CODE_DEFINE=code_define,
    CODE_NAMELIST=code_namelist,
    CODE_DEFAULT=code_default,
    CODE_REWIND=code_rewind,
    CODE_DUMP=code_dump,
    CODE_BCAST=code_bcast,
))