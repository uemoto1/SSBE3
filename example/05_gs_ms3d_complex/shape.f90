program main
    implicit none
    integer, parameter :: nx_m = 100
    integer, parameter :: ny_m = 100
    integer, parameter :: nz_m = 1
    real(8), parameter :: hx_m = 100.0d0
    real(8), parameter :: hy_m = 100.0d0
    real(8), parameter :: hz_m = 100.0d0

    real(8), parameter :: x0 = 5000.0d0, y0 = 5000.0d0, r0 = 3000.0d0

    integer :: ix_m, iy_m, iz_m, icount, ncount
    integer :: ibuf(4, nx_m*ny_m*nz_m)
    real(8) :: x, y, z, r

    icount = 0
    do ix_m = 1, nx_m
        do iy_m = 1, ny_m
            do iz_m = 1, nz_m
                x = ix_m * hx_m
                y = iy_m * hy_m
                z = iz_m * hz_m

                r = sqrt((x-x0)**2 + (y-y0)**2)
                if (r <= r0) then
                    icount = icount + 1
                    ibuf(1, icount) = ix_m
                    ibuf(2, icount) = iy_m
                    ibuf(3, icount) = iz_m
                    ibuf(4, icount) = -1
                else if (x > x0) then
                    icount = icount + 1
                    ibuf(1, icount) = ix_m
                    ibuf(2, icount) = iy_m
                    ibuf(3, icount) = iz_m
                    ibuf(4, icount) = 1
                end if
            end do
        end do
    end do
    ncount = icount

    open(99, file="shape.txt", action="write")
    write(99, "(99(i6))") ncount
    do icount = 1, ncount
        write(99, "(99(i6))") ibuf(:, icount)
    end do
    close(99)
    
    stop
end program main
