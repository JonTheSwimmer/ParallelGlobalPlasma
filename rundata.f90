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
    use parallel_prop
    implicit none

    integer :: icount, fnum
    character(len=8) :: filenum, thread_count
    integer :: runnum
    real*8 :: time
    !time tracking 
    real*8 :: start, finish , time_para
    integer :: count1, count2, count_rate, count_max
    
    icount = iargc()
    if (icount.eq.2) then
        call getarg(1, filenum)
        read (filenum, '(I8)') fnum
        call getarg(2, thread_count)
        read (thread_count, "(I8)") N_threads
    else
        write (*,*) "wrong number of arguments"
        stop
    end if
    time = 3.6d3 * 7.d0 !expected runtime in seconds
    runnum = floor(N_threads * time * 350)
    !set B_dip
    B_dip = 2.d0
    call initialize(fnum)

    call cpu_time(start)
    call system_clock(count1, count_rate, count_max)
    print *, "writing"
    !call propagate_parallel(runnum)
    call propagate_write(runnum)
    call cpu_time(finish)
    call system_clock(count2, count_rate, count_max)
    print *, sum(cell_time)
    print *, sum(cell_bins)
    print *, stat_emit 
    print *, stat_surf
    print *, stat_conv
    print *, sum(surf_flux)
    print *, N_prop
    print *, "series"
    call print_time(finish-start)
    print *, "parallel"
    time_para = 1.d0 * (count2 - count1)/count_rate
    call print_time(time_para)
    call output_data(fnum)
    write (*, '(A, I0, A)') "Run ", fnum, " complete"
    end program main
