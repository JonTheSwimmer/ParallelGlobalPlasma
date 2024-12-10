!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main code to run and store photon trajectory data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    program main
    use parameters_base
    use parameters
    use routines
    use spectrum
    use grid
    use trajectory
    use emission
    implicit none

    integer :: icount, fnum
    character(len=8) :: filenum
    integer :: runnum
    real*8 :: time
    
    icount = iargc()
    if (icount.eq.1) then
        call getarg(1, filenum)
        read (filenum, '(I8)') fnum
    else
        write (*,*) "file number not specified"
        stop
    end if
    time = 3.6d3 * 1.6d1 !expected runtime in seconds
    runnum = floor(time * 450)
    call initialize(fnum)
    call propagate(runnum)
    call output_data(fnum)
    write (*, '(A, I0, A)') "Run ", fnum, " complete"
    end program main
