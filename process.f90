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
    use photon_split
    use OMP_LIB
    implicit none
    !for both n_pair and N_create, the southern polar cap is stored at Nphi = 2 
    real*8, dimension(Nr+1, Nphi) :: n_psurf !density of pair plasma at the surface
    real*8, dimension(Nr+1, Nphi) :: N_create, C_surf, C_out, C_ann !pair creation/loss rate for each bundle
    real*8, dimension(Nr+1, 2*Ntheta, Nphi) :: threshold_arr, phot_highE
    real*8, dimension(Nr+1, 2, Nphi) :: surf_temp
    real*8, dimension(Nr+1, Nphi, Nomega) :: create_splits
    real*8, dimension(Nomega) :: loss_arr, emit_arr, phot_ratio
    real*8 :: E_ploss, lfactor, beta !total erg/s of photons being converted to pairs
    real*8 :: run_time !"real" time for simuation
    real*8, dimension(Nomega) :: Edot
    integer :: iR, iTS, iT, iW, iP, i_Omin, ifilt
    integer :: iD, batchnum
    real*8 :: a, b, c, l_tot
    character(len=8) :: num
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
    print *, "initialized"
    batchnum = 710
    combined = .false.
    call load_data(1700)
    !do iD = 101, 120
    !    call load_data(iD)
    !end do
    !do iD = 201, 220
    !    call load_data(iD)
    !end do
    if (combined) then
        call output_data(2000+batchnum)
    end if
    call init_photon()
    print *, "pre-filter"
    print *, sum(cell_time)
    print *, sum(cell_bins)
    !determine how much pair production is coming from each energy bin, normalized to luminosity of each bin
    call rate_split_para(Edot)
    ifilt = Nomega
    do iW = 1, Nomega
        loss_arr(iW) = sum(cell_loss_split(:, :, :, iW))
        emit_arr(iW) = rate_total(omega01(iW), omega02(iW))
        phot_ratio(iW) = loss_arr(iW) / emit_arr(iW)
        if (phot_ratio(iW) .gt. 1.5d0) then
            ifilt = min(ifilt, iW-1)
        end if
    end do
    print *, "filtered at", ifilt
    E_ploss = sum(Edot(1:ifilt))
    cell_create = sum(cell_create_split(:, :, :, 1:ifilt), DIM=4) 
    !save pair creation rate for each flux bundle, split by energy in case of re-cutting
    create_splits = sum(cell_create_split, DIM = 2)
    !determine average energy of particles
    lfactor = (E_ploss / (sum(cell_create) * 2.d0* m_e))
    run_time = N_prop / N_emission
    !determine bin when energy > m_e c^2
    do iW = 1, Nomega
        if (omega_mid(iW) .lt. m_e) then
            i_Omin = min(Nomega, iW+1)
        end if
    end do
    if (lfactor .lt. 1.d0) then
        print *, "energy is wrong!"
    end if
    beta = sqrt(1.d0 - lfactor**(-2.d0))

    do iR = 1, Nr+1
        if (iR .lt. Nr+1) then
            do iP = 1, Nphi
                N_create(iR, iP) = sum(cell_create(iR, :, iP))
            end do
        else
            N_create(iR, 1) = sum(cell_create(iR, 1:Ntheta, 1))
            N_create(iR, 2) = sum(cell_create(iR, (Ntheta+1):2*Ntheta, 1))
        end if
        C_surf(iR, :) = ann_surface(iR, 1, 1.d0)
        C_out(iR, :) = ann_outflow(iR, 1, 1.d0, beta)
        C_ann(iR, :) = ann_average(iR, 1, beta)
    end do
    do iR = 1, Nr+1
        do iP = 1, Nphi
            a = C_ann(iR, iP)
            b = C_surf(iR, iP) + C_out(iR, iP)
            c = N_create(iR, iP)
            !density at surface in n/R_NS**3
            n_psurf(iR, iP) = (-b + sqrt(b**2.d0 + 4.d0*a*c))/(2.d0*a) 
        end do
    end do
    do iR = 1, Nr+1
        do iTS = 1, 2*Ntheta
            do iP = 1, Nphi
                if (cell_volume(iR, iTS, iP) .gt. 0) then
                    threshold_arr(iR, iTS, iP) = threshold(iR, iTS, iP)
                    l_tot = sum(cell_bins(iR, iTS, iP, :, :, i_Omin:))
                    phot_highE(iR, iTS, iP) = (l_tot / (c_light * cell_volume(iR, iTS, iP)))
                end if 
            end do
        end do
    end do
    phot_highE = phot_highE / (run_time * R_NS**3.d0) !convert to cm^-3
    do iR = 1, Nr+1
        do iTS = 1, 2
            do iP = 1, Nphi
                surf_temp(iR, iTS, iP) = flux2temp(iR, iTS, iP)     
            end do
        end do
    end do
101 format(*(g0, :, ", "))
    write (num, "(I0)") batchnum
    call file_open(15, trim(num)//"_N_create.dat")
    call file_open(16, trim(num)//"_C_surf.dat")
    call file_open(17, trim(num)//"_C_ann.dat")
    call file_open(18, trim(num)//"_n_surf.dat")
    call file_open(19, trim(num)//"_C_out.dat")
    call file_open(20, trim(num)//"_threshold.dat")
    call file_open(21, trim(num)//"_surf_temp.dat")
    call file_open(23, trim(num)//"_phot_highE.dat")
    call file_open(24, trim(num)//"_lum_ratio.dat")
    call file_open(25, trim(num)//"_process_params.dat")
    call file_open(26, trim(num)//"_create_splits.dat")
    write(24, 101) phot_ratio
    write(25, 101) luminosity, B_dip
    close(25)
    close(24)
    do iR = 1, Nr+1
        write(15, 101) N_create(iR, :)
        write(16, 101) C_surf(iR, :)
        write(19, 101) C_out(iR, :)
        write(17, 101) C_ann(iR, :)
        write(18, 101) n_psurf(iR, :)
        write(21, 101) surf_temp(iR, 1, :)
        write(21, 101) surf_temp(iR, 2, :)
        do iTS = 1, 2*Ntheta
            write(20, 101) threshold_arr(iR, iTS, :)
            write(23, 101) phot_highE(iR, iTS, :)
        end do
        do iP = 1, Nphi
            write(26, 101) create_splits(iR, iP, :)
        end do
    end do
    close(15); close(16); close(17); close(18); 
    close(19); close(20); close(21); close(23); close(26)
    !summary statistics
    call file_open(22, trim(num)//"_summary_stats.dat")
    write (22, 101) "Total E_dot of gamma ann = ", E_ploss
    write (22, 101) "Average gamma of pair = ", lfactor
    write (22, 101) "Plasma beta =", beta 
    write (22, 101) "Total photon heating =", sum(surf_flux)/run_time
    write (22, 101) "Total positron heating =", sum(n_psurf*C_surf*2*lfactor*m_e)
    write (22, 101) "Volumetric annihilation =", sum(n_psurf**2 *C_ann*2*lfactor*m_e)
    close(22)
    call output_rate(batchnum)
    call output_emit_param(batchnum) 
    call output_grid(batchnum)
    call single_photon_filter()
    print *, "ratio after filtering = ", sum(cell_bins)/sum(cell_time)
    end program main
