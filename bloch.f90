module time_evolution
    implicit none

    type :: rt_info
        integer :: nk
        integer :: nstate
        real(8) :: dt
        complex(8), allocatable :: rho(:, :, :)
    end type

contains

subroutine init_bloch(rt, gs)
    implicit none
    type(rt_info), intent(out) :: rt
    type(gs_info), intent(in) :: gs
    rt%nk = gs%nk
    rt%nstate = gs%nstate
    allocate(rt%rho(rt%nstate, rt%nstate, rt%nk))
    return
end subroutine init_bloch

subroutine dt_evolve_bloch(rt, gs, dt, Ac0)
    use ground_state
    implicit none
    type(rt_info), intent(inout) :: rt
    type(gs_info), intent(in) :: gs
    real(8), intent(in) :: dt, Ac0(3)

contains

subroutine calc_drho(drho, rho, ik, Ac)
    implicit none
    complex(8), intent(out) :: drho(rt%nstate, rt%nstate, rt%nk)
    complex(8), intent(in) :: rho(rt%nstate, rt%nstate, rt%nk)
    integer, intent(in) :: ik
    real(8), intent(in) :: Ac(3)
    integer :: i, j, l
    
    drho(:, :, :) = 0.0d0
    do ik = 1, rt%nk
        do i = 1, rt%nstate
            do j = 1, rt%nstate
                do l = 1, rt%nstate
                    drho(i, j, ik) = drho_k(i, j, ik) + dcmplx(0.0, -1.0) * Ac(n) &
                        & * ((gs%pmatrix(n, i, l, ik) + gs%rvnl(n, i, l, ik)) * rho_k(l, j) &
                        & -  rho_k(i, l) * (gs%pmatrix(n, l, j, ik) + gs%rvnl(n, l, j, ik)))
                end do
            end do
        end do
    end do
    return
end subroutine
end subroutine


end module time_evolution