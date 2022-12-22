program shape
    implicit none
    integer, parameter :: nx_m = 100
    integer, parameter :: ny_m = 100
    integer, parameter :: nz_m = 1
    real(8), parameter :: hx_m = 100.0d0
    real(8), parameter :: hy_m = 100.0d0
    real(8), parameter :: hz_m = 100.0d0

    real(8), parameter :: x0 = 5000.0
    real(8), parameter :: y0 = 5000.0
    real(8), parameter :: r0 = 3000.0

    integer :: ix, iy, iz, imacro, nmacro
    real(8) :: x, y, z, r
    integer :: itbl_macro_coord(1:3, nx_m*ny_m*nz_m)

    imacro = 0

    do iz = 1, nz_m
        do iy = 1, ny_m
            do ix = 1, nx_m
                x = hx_m * ix
                y = hy_m * iy
                z = hz_m * iz

                r = sqrt((x-x0)**2 + (y-y0)**2)

                if (r < r0) then
                    imacro = imacro + 1
                    itbl_macro_coord(1, imacro) = ix
                    itbl_macro_coord(2, imacro) = iy
                    itbl_macro_coord(3, imacro) = iz
                end if
            end do
        end do
    end do

    nmacro = imacro

    open(100, file="shape.txt", action="write")
    write(100, "(i6)") nmacro
    do imacro = 1, nmacro
        write(100, "(4i6)") itbl_macro_coord(1:3, imacro), -1
    end do

    stop "bye!"
end program shape

