program main
    use mpi
    use omp_lib
    use communication
    use input_parameter
    use multiscale
    use realtime
    implicit none
    integer :: ierr, icomm, irank, nthread, nproc

    call MPI_INIT(ierr)
    call comm_get_globalinfo(icomm, irank, nproc)

    !$omp parallel
    !$omp master
    nthread = omp_get_num_threads()
    !$omp end master
    !$omp end parallel

    if (irank == 0) then
        write(*, "(a)") "# SSBE3"
        write(*, "(a,x,a,x,a)") "# built:", __DATE__, __TIME__
        write(*, "(a)") "# -----"
        write(*, "(a)") "# Parallelization:"
        write(*, "(a,i6)") "# number of MPI processes=", nproc
        write(*, "(a,i6)") "# number of OMP threads=", nthread
    end if

    call read_input(icomm)

    select case (trim(theory))
    case ("sbe_rt")
        call realtime_main(icomm)
    case ("sbe_ms")
        call multiscale_main(icomm)
    end select

    if (irank == 0) then
        write(*, "(a)") "# Successfully finished"
    end if

    call MPI_FINALIZE(ierr)

    stop "Bye!"
end program 
