module sbe_solver
    use sbe_gs
    use math_constants, only: pi, zi
    use communication, only: comm_get_groupinfo, comm_summation
    use util, only: split_range
    implicit none



    type s_sbe_bloch_solver
        !k-points for real-time SBE calculation
        integer :: nk, nb
        integer :: ik_max, ik_min
        complex(8), allocatable :: rho(:, :, :)
    end type



contains



subroutine init_sbe(sbe, gs, nb_sbe, icomm)
    use util
    use communication
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    integer, intent(in) :: nb_sbe
    integer, intent(in) :: icomm
    integer :: ik, ib, nk_proc, irank, nproc, ierr
    integer, allocatable :: itbl_min(:), itbl_max(:)

    call comm_get_groupinfo(icomm, irank, nproc)

    sbe%nk = gs%nk
    sbe%nb = nb_sbe

    allocate(itbl_min(0:nproc-1), itbl_max(0:nproc-1))

    call split_range(1, sbe%nk, nproc, itbl_min, itbl_max)
    sbe%ik_min = itbl_min(irank)
    sbe%ik_max = itbl_max(irank)

    allocate(sbe%rho(1:sbe%nb, 1:sbe%nb, sbe%ik_min:sbe%ik_max))
    
    sbe%rho(:, :, :) = 0d0
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            sbe%rho(ib, ib, ik) = gs%occup(ib, ik)
        end do
    end do
end subroutine


subroutine calc_current_bloch(sbe, gs, Ac, jmat, icomm)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(out) :: jmat(1:3)
    integer, intent(in) :: icomm
    integer :: ik, idir, ib, jb, ierr
    complex(8) :: tmp1(1:3), tmp(1:3)

    tmp1(1:3) = 0d0

    !$omp parallel do default(shared) private(ik,ib,jb,idir) reduction(+:tmp1)
    do ik = sbe%ik_min, sbe%ik_max
        do idir = 1, 3
            do ib = 1, sbe%nb
                do jb = 1, sbe%nb
                    tmp1(idir) = tmp1(idir) + gs%kweight(ik) * sbe%rho(jb, ib, ik) * ( &
                        & gs%p_matrix(ib, jb, idir, ik) & !+ zI * gs%rvnl_matrix(ib, jb, idir, ik) &
                        & )
                end do
            end do
        end do
    end do
    !$omp end parallel do
    call comm_summation(tmp1, tmp, 3, icomm)

    jmat(:) = (real(tmp(1:3)) / sum(gs%kweight(:)) &
        & + Ac * calc_trace(sbe, gs, sbe%nb, icomm)) / gs%volume    

    return
end subroutine calc_current_bloch


subroutine dt_evolve_bloch(sbe, gs, Ac, dt)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs), intent(inout) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(in) :: dt
    complex(8), parameter :: zi = dcmplx(0d0, 1d0)
    integer :: nb, nk, ik

    complex(8) :: hrho1_k(1:sbe%nb, 1:sbe%nb)
    complex(8) :: hrho2_k(1:sbe%nb, 1:sbe%nb)
    complex(8) :: hrho3_k(1:sbe%nb, 1:sbe%nb)
    complex(8) :: hrho4_k(1:sbe%nb, 1:sbe%nb)

    nb = sbe%nb
    nk = sbe%nk

    !$omp parallel do default(shared) private(ik, hrho1_k, hrho2_k, hrho3_k, hrho4_k)
    do ik = sbe%ik_min, sbe%ik_max
        call calc_hrho_bloch_k(ik, sbe%rho(:, :, ik), hrho1_k)
        call calc_hrho_bloch_k(ik, hrho1_k, hrho2_k)
        call calc_hrho_bloch_k(ik, hrho2_k, hrho3_k)
        call calc_hrho_bloch_k(ik, hrho3_k, hrho4_k)

        sbe%rho(:, :, ik) = sbe%rho(:, :, ik) + hrho1_k * (- zi * dt)
        sbe%rho(:, :, ik) = sbe%rho(:, :, ik) + hrho2_k * (- zi * dt) ** 2 * (1d0 / 2d0)
        sbe%rho(:, :, ik) = sbe%rho(:, :, ik) + hrho3_k * (- zi * dt) ** 3 * (1d0 / 6d0)
        sbe%rho(:, :, ik) = sbe%rho(:, :, ik) + hrho4_k * (- zi * dt) ** 4 * (1d0 / 24d0)
    end do
    return

contains


    !Calculate [H, rho] commutation:
    subroutine calc_hrho_bloch_k(ik, rho_k, hrho_k)
        implicit none
        integer, intent(in) :: ik
        complex(8), intent(in) :: rho_k(nb, nb)
        complex(8), intent(out) :: hrho_k(nb, nb)
        integer :: idir
        !hrho = hrho + Ac(t) * (p * rho - rho * p)
        hrho_k(1:nb, 1:nb) = gs%delta_omega(1:nb, 1:nb, ik) * rho_k(1:nb, 1:nb)
        do idir=1, 3 !1:x, 2:y, 3:z
            ! hrho(1:nb, 1:nb, ik) = hrho(1:nb, 1:nb, ik) + Ac(idir) * (&
            ! & + matmul(gs%p_matrix(1:nb, 1:nb, idir, ik), rho(1:nb, 1:nb, ik)) &
            ! & - matmul(rho(1:nb, 1:nb, ik), gs%p_matrix(1:nb, 1:nb, idir, ik)) &
            ! & )

            call ZGEMM("N","N", sbe%nb, sbe%nb, sbe%nb, &
                dcmplx(+Ac(idir), 0d0), &
                gs%p_matrix(:, :, idir, ik),sbe%nb, &
                rho_k(:, :), sbe%nb, &
                dcmplx(1d0, 0d0), hrho_k(:, :),sbe%nb)

            call ZGEMM("N","N", sbe%nb, sbe%nb, sbe%nb, &
                dcmplx(-Ac(idir), 0d0), &
                rho_k(:, :), sbe%nb, &
                gs%p_matrix(:, :, idir, ik),sbe%nb, &
                dcmplx(1d0, 0d0), hrho_k(:, :), sbe%nb)
                    
        end do !idir
        return
    end subroutine calc_hrho_bloch_k
end subroutine

function calc_trace(sbe, gs, nb_max, icomm) result(tr)
    use communication
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    integer, intent(in) :: icomm
    integer, intent(in) :: nb_max
    real(8) :: tr

    integer :: ik, ib
    real(8) :: tmp, tmp1
    
    tmp1 = 0d0
    !$omp parallel do default(shared) private(ik, ib) reduction(+: tmp1) collapse(2) 
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, nb_max
            tmp1 = tmp1 + real(sbe%rho(ib, ib, ik)) * gs%kweight(ik)
        end do
    end do
    !$omp end parallel do
    call comm_summation(tmp1, tmp, icomm)
    tr = tmp / sum(gs%kweight)
    
    return
end function calc_trace


function calc_energy(sbe, gs, Ac, icomm) result(energy)
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    integer, intent(in) :: icomm
    real(8), intent(in) :: Ac(1:3)
    integer :: ik, ib, jb, idir
    real(8) :: tmp1, tmp, energy
    ! real(8) :: kvec(1:3)
    tmp1 = 0d0
    !$omp parallel do default(shared) private(ik, ib, jb, idir) reduction(+: tmp1)
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            do idir = 1, 3
                do jb = 1, sbe%nb
                    tmp1 = tmp1 &
                        & + Ac(idir) * real(sbe%rho(ib, jb, ik) * gs%p_matrix(jb, ib, idir, ik)) * gs%kweight(ik)
                end do
            end do
            tmp1 = tmp1 &
                & + real(sbe%rho(ib, ib, ik)) * ( &
                & + gs%eigen(ib, ik) &
                !& + dot_product(kvec(:), Ac(:))
                & + 0.5 * dot_product(Ac, Ac) &
                & ) * gs%kweight(ik)
        end do
    end do
    !$omp end parallel do
    call comm_summation(tmp1, tmp, icomm)
    energy = tmp / sum(gs%kweight)

    return 
end function calc_energy


end module



