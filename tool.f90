module tool
    implicit none
contains
subroutine calc_dist_num(n, m, dist_n)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(out) :: dist_n(1:m)
    integer :: i, mod_n_m

    mod_n_m = mod(n, m)
    do i = 1, m
        if (i < mod_n_m) then
            dist_n(i) = n / m + 1
        else
            dist_n(i) = n / m
        end if
    end do
    return
end subroutine calc_dist_num

subroutine calc_dist_range(nmin, nmax, m, itbl_min, itbl_max)
    implicit none
    integer, intent(in) :: nmin, nmax
    integer, intent(in) :: m
    integer, intent(out) :: itbl_min(m), itbl_max(m)
    integer :: ntmp(m), itmp, i

    call calc_dist_num(nmax-nmin+1, m, ntmp)
    itmp = nmin
    do i = 1, m
        itbl_min(i) = itmp
        itbl_max(i) = itmp + (ntmp(i) - 1)
        itmp = itmp + 1
    end do
    return
end subroutine calc_dist_range

end module tool
