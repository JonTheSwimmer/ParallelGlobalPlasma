
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Shorthand for opening and checking if a file already exists
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine file_open(unitnum, filename)
    implicit none
    integer, intent(in) :: unitnum
    character(len=*), intent(in) :: filename
    logical :: exist
    inquire(file=filename, exist = exist)
    if (exist) then
        open(unit = unitnum, file = filename,Status='old', position="append", action="write")
    else
        open(unit = unitnum, file = filename, status = 'new')
    endif
    return
    end subroutine file_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   output the data array (in my case cell_time) to a file
!   can adjust name by runnum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine output_data(runnum)
    implicit none
    integer, intent(in) :: runnum
    integer :: iR, iT, iP, iK, iW
    character(len=8) :: num ! format descriptor
    write (num, "(I0)") runnum
    call file_open(2, 'data/cell_time'//trim(num)//'.dat')
101 format(*(g0, :, ", "))
    do iR = 1, Nr+1 !iterate over array dimensions
        do iT = 1, 2*Ntheta
            write (2, 101) cell_time(iR, iT, :)
        end do
    end do
    close (2)
    return
    end subroutine output_data

 
