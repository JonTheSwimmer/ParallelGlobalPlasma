    program main
    use parameters_base
    use parameters
    use routines
    use spectrum
    use grid
    use trajectory
    use emission
    implicit none

    integer :: icount, fnum1, fnum2
    character(len=8) :: file1, file2, outstr
    integer :: iD, outnum
    real*8 :: time

    icount = iargc()
    if (icount.eq.3)then
        call getarg(1, file1)
        call getarg(2, file2)
        call getarg(3, outstr)
        read (file1, '(I8)') fnum1
        read (file2, '(I8)') fnum2
        read (outstr, '(I8)') outnum
    else
        write (*,*) "file numbers not specified"
        stop
    end if
    print *, fnum1, fnum2
    do iD = fnum1, fnum2
        call load_data(iD)
    end do
    print *, N_prop

    call output_data(outnum)
     
    end program main
