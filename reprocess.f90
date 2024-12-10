!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    program main
    use parameters_base
    use parameters
    use routines
    use spectrum
    use grid
    use trajectory
    use emission
    use photon
    use OMP_LIB
    implicit none
    !for both n_pair and N_create, the southern polar cap is stored at Nphi = 2 
    real*8, dimension(Nr+1, Nphi) :: n_psurf !density of pair plasma at the surface
    real*8, dimension(Nr+1, Nphi) :: N_create, C_surf, C_out, C_ann !pair creation/loss rate for each bundle
    real*8, dimension(Nr+1, 2*Ntheta, Nphi) :: threshold_arr, phot_highE
    real*8, dimension(Nr+1, 2, Nphi) :: surf_temp
    real*8, dimension(Nomega) :: lum_arr
    real*8 :: E_ploss, lfactor, beta !total erg/s of photons being converted to pairs
    real*8 :: run_time !"real" time for simuation
    integer :: iR, iTS, iT, iW, iP, i_Omin, ifilt, i
    integer :: iD, batchnum, readnum, rawnum
    real*8 :: a, b, c, l_tot
    character(len=8) :: num, readdata, rawdata
    character (len=1) :: text_ignore
    logical :: combined 
    n_psurf = 0
    N_create = 0
    C_out = 0
    C_surf = 0
    C_ann = 0
    threshold_arr = 0
    phot_highE = 0
    E_ploss = 0
    call initialize(2)
    call init_photon()
    print *, "initialized"
    batchnum = 710
    readnum = 710
    !read in summary data
    write (num, "(I0)") readnum
    open(2, file='data/'//trim(num)//"_summary_stats.dat")
    read (2, *) (text_ignore, i=1, 6), E_ploss
    print *, E_ploss
    read (2, *) (text_ignore, i=1, 5), lfactor
    print *, lfactor
    read (2, *) (text_ignore, i=1, 3), beta
    print *, beta
    close(2)
    if (lfactor .lt. 1.d0) then
        print *, "energy is wrong!"
    end if
    beta = sqrt(1.d0 - lfactor**(-2.d0))

    open(5, file = "data/"//trim(num)//"_N_create.dat")
    do iR = 1, Nr+1
        read(5, *) N_create(iR, :)
        C_surf(iR, :) = ann_surface(iR, 1, 1.d0)
        C_out(iR, :) = ann_outflow(iR, 1, 1.d0, beta)
        C_ann(iR, :) = ann_average(iR, 1, beta)
    end do
    close(5)
    do iR = 1, Nr+1
        do iP = 1, Nphi
            a = C_ann(iR, iP)
            b = C_surf(iR, iP) + C_out(iR, iP)
            c = N_create(iR, iP)
            !density at surface in n/R_NS**3
            n_psurf(iR, iP) = (-b + sqrt(b**2.d0 + 4.d0*a*c))/(2.d0*a) 
        end do
    end do
101 format(*(g0, :, ", "))
    write (num, "(I0)") batchnum
    call file_open(15, trim(num)//"_N_create.dat")
    call file_open(16, trim(num)//"_C_surf.dat")
    call file_open(17, trim(num)//"_C_ann.dat")
    call file_open(18, trim(num)//"_n_surf.dat")
    call file_open(19, trim(num)//"_C_out.dat")
    !call file_open(24, trim(num)//"_lum_ratio.dat")
    call file_open(25, trim(num)//"_process_params.dat")
    !write(24, 101) lum_arr
    write(25, 101) luminosity, B_dip
    close(25)
    !close(24)
    do iR = 1, Nr+1
        write(15, 101) N_create(iR, :)
        write(16, 101) C_surf(iR, :)
        write(19, 101) C_out(iR, :)
        write(17, 101) C_ann(iR, :)
        write(18, 101) n_psurf(iR, :)
    end do
    close(15); close(16); close(17); close(18)
    close(19)
    !summary statistics
    call file_open(22, trim(num)//"_summary_stats.dat")
    write (22, 101) "Total E_dot of gamma ann = ", E_ploss
    write (22, 101) "Average gamma of pair = ", lfactor
    write (22, 101) "Plasma beta =", beta 
    write (22, 101) "Total photon heating =", sum(surf_flux)/run_time
    write (22, 101) "Total positron heating =", sum(n_psurf*C_surf*2*lfactor*m_e)
    write (22, 101) "Volumetric annihilation =", sum(n_psurf**2 *C_ann*2*lfactor*m_e)
    close(22)
    end program main
