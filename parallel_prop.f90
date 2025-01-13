!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Module for parallel photon emission and assignment on the grid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module parallel_prop
    use parameters_base
    use parameters
    use routines
    use grid
    use spectrum
    use trajectory
    use emission
    use OMP_LIB
    implicit none 

    contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine path_create(l_phot, coord_phot, k_phot, l0, ln, Npath, omega, coord,side_ind, &
         ind_phot, time_phot, entry_index, ind_surf, flux_surf, surf_index)
    !assign a photon path onto the grid using bisection
    implicit none
    !input variables
    real*8, intent(in) :: l0, ln, omega
    real*8, dimension(2*N_spoints) :: l_phot
    real*8, dimension(2*N_spoints, 3), intent(in) :: k_phot, coord_phot
    integer, intent(in) :: Npath, side_ind
    real*8, dimension(5) :: coord
    !in/out variables
    integer :: ind_phot(N_para_list, 6), entry_index, ind_surf(N_para_list, 4), surf_index
    real*8 :: time_phot(N_para_list), flux_surf(N_para_list)
    !local variables
    real*8, dimension(2*N_spoints, 3) :: coord2_phot, k2_phot
    integer :: N
    real*8 :: start_l, end_l, lower, upper, mid, tol, scale_l
    !tolerance parameters
    real*8 :: rtol, atol
    !energy parameters
    real*8 :: omega_start, omega_end, r_e, r_i
    real*8, dimension(3) :: k_start, k_end, cpt, pt
    integer :: iR, iP, indices(3), cell_ind(3), calls, k_ind(3)
    logical :: cont, cont2
    real*8 :: drdL
    101 format(*(g0, :, ", "))
    rtol = 1.d-4
    atol = 1.d-6
    N = Npath
    do iR = 1, 3
        call spline(l_phot, coord_phot(:, iR), N, k_phot(1, iR), k_phot(N, iR), coord2_phot(:, iR))
        call spline(l_phot, k_phot(:, iR), N, coord2_phot(1, iR), coord2_phot(N, iR), k2_phot(:, iR))
    end do
    !set initial cell and photon direction
    indices = 0
    start_l = l0
    call traj_place(l_phot, coord_phot, coord2_phot, N, start_l, indices, r_e)
    call triple_splint(l_phot, k_phot, k2_phot, N, start_l, k_start)
    do while ((sum(indices) .lt. 0) .and. (start_l .lt. ln))
        print *, "moving forward"
        start_l = start_l + rtol * (ln-l0)
        do iR = 1, 3
            call splint(l_phot, coord_phot(:, iR), coord2_phot(:, iR), N, start_l, cpt(iR), .true.)
        end do
        call spherical(cpt(1), cpt(2), cpt(3), pt(1), pt(2), pt(3))
        drdL = sum(cpt(:) * k_start(:)) / pt(1)
        call traj_place(l_phot, coord_phot, coord2_phot, N, start_l, indices, r_e)
        print *, start_l, r_e, pt, coord(1), drdL, coord(4)
    end do
    call triple_splint(l_phot, k_phot, k2_phot, N, start_l, k_start)
    omega_start = omega
    cont = .true.
    calls = 0
    do while (cont)
        cell_ind = indices
        cont2 = .true.
        lower = start_l 
        upper = ln
        !bisection method to determine end_l of cell
        do while (cont2)
            calls = calls + 1
            if (calls .ge. 1d6) then
                cont = .false.
                cont2 = .false.    
                call file_open(17, "error.txt")
                write (17, *) "bisection failed"
                write (17, 101) coord
                write (17, *) cell_ind
                write (17, *) upper, lower, start_l
                close (17)
            end if
            end_l = (upper + lower) / 2.d0
            call traj_place(l_phot, coord_phot, coord2_phot, N, end_l, indices, r_i)
            if (all(cell_ind(:) .eq. indices(:))) then
                lower = end_l
            else
                upper = end_l
            end if
            if (((upper - lower)/(upper-start_l) .lt. rtol) .or. ((upper - lower) .lt. atol)) then 
                cont2 = .false.
            end if
        end do
        omega_end = omega_shift(r_e, r_i, omega)
        call triple_splint(l_phot, k_phot, k2_phot, N, end_l, k_end)
        call place_k_arr((k_start + k_end), (omega_start + omega_end)/2.d0, k_ind)
        !assign path length to grid array
        if (inBounds_arr(cell_ind, k_ind)) then 
            ind_phot(entry_index, :) = (/cell_ind(1), cell_ind(2), cell_ind(3), &
                                         k_ind(1), k_ind(2), k_ind(3)/)
            time_phot(entry_index) = end_l - start_l
            entry_index = entry_index + 1
        end if
        !check if photon is above single photon conversion threshhold
        if (single_convert(omega_start, omega_end, k_start, &
                            k_end, cell_ind(1), cell_ind(2), cell_ind(3))) then
            cont = .false. 
            ind_surf(surf_index, :) = (/0, 0, 0, side_ind /)
            flux_surf(surf_index) = 0.d0
            surf_index = surf_index + 1
            return
        end if
        call traj_place(l_phot, coord_phot, coord2_phot, N, upper, indices, r_i)
        !check if path has left the grid or end of path has been reached
        if ((sum(indices) < 0) .or. (abs((ln - upper) / (ln-l0)) < tol)) then
            cont = .false.
        end if
        !if all ending criteria are not met, set up for next cell
        start_l = end_l
        k_start = k_end
        omega_start = omega_end
        cell_ind = indices
    end do
    !check if end point of path is at the surface of the star
    call triple_splint(l_phot, coord_phot, coord2_phot, N, ln, cpt)
    call spherical(cpt(1), cpt(2), cpt(3), pt(1), pt(2), pt(3))
    if (pt(1) .lt. 1.1d0) then
        !check that photon is travelling inward
        call place_surf(pt(2), pt(3), indices(1), indices(2), indices(3))
        call triple_splint(l_phot, k_phot, k2_phot, N, ln, k_end)
        drdL = sum(cpt(:) * k_end(:)) / pt(1)
        if (drdL .lt. 0.d0) then
            ind_surf(surf_index, :) = (/indices(1), indices(2), indices(3), side_ind/)
            flux_surf(surf_index) = omega_shift(r_e, 1.d0, omega)
            surf_index = surf_index + 1
        end if
    end if
    return 
    end subroutine path_create
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   propagate photons directly to grid, with OMP CRITICAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine grid_prop(l_phot, coord_phot, k_phot, l0, ln, Npath, omega, coord, side_ind)
    !assign a photon path onto the grid using bisection
    implicit none
    !input variables
    real*8, intent(in) :: l0, ln, omega
    real*8, dimension(2*N_spoints) :: l_phot
    real*8, dimension(2*N_spoints, 3), intent(in) :: k_phot, coord_phot
    integer, intent(in) :: Npath, side_ind
    real*8, dimension(5) :: coord
    !local variables
    real*8, dimension(2*N_spoints, 3) :: coord2_phot, k2_phot
    integer :: N
    real*8 :: start_l, end_l, lower, upper, mid, tol, scale_l
    !energy parameters
    real*8 :: omega_start, omega_end, r_e, r_i
    real*8, dimension(3) :: k_start, k_end, cpt, pt
    integer :: iR, iP, indices(3), cell_ind(3), calls, k_ind(3)
    logical :: cont, cont2
    !surface check parameters
    real*8 ::  T01, T02, P1, P2, mu, n_hat(3)
    real*8 :: drdL
    101 format(*(g0, :, ", "))
    tol = 1.d-4
    N = Npath
    do iR = 1, 3
        call spline(l_phot, coord_phot(:, iR), N, k_phot(1, iR), k_phot(N, iR), coord2_phot(:, iR))
        !call spline(l_phot, k_phot(:, iR), N, dk0(iR), dkn(iR), k2_phot(:, iR))
        call spline(l_phot, k_phot(:, iR), N, coord2_phot(1, iR), coord2_phot(N, iR), k2_phot(:, iR))
    end do
    !print *, sum(k2_phot)
    !set initial cell and photon direction
    indices = 0
    start_l = l0
    call traj_place(l_phot, coord_phot, coord2_phot, N, start_l, indices, r_e)
    do while ((sum(indices) .eq. 0) .and. (start_l .lt. ln))
        print *, "moving forward"
        start_l = start_l + tol * (ln-l0)
        call traj_place(l_phot, coord_phot, coord2_phot, N, start_l, indices, r_e)
    end do
    call triple_splint(l_phot, k_phot, k2_phot, N, start_l, k_start)
    omega_start = omega
    cont = .true.
    calls = 0
    do while (cont)
        cell_ind = indices
        cont2 = .true.
        lower = start_l 
        upper = ln
        !bisection method to determine end_l of cell
        do while (cont2)
            calls = calls + 1
            if (calls .ge. 1d6) then
                cont = .false.
                cont2 = .false.    
                call file_open(17, "error.txt")
                write (17, *) "bisection failed"
                write (17, 101) coord
                write (17, *) cell_ind
                close (17)
            end if
            end_l = (upper + lower) / 2.d0
            call traj_place(l_phot, coord_phot, coord2_phot, N, end_l, indices, r_i)
            if (all(cell_ind(:) .eq. indices(:))) then 
                lower = end_l
            else
                upper = end_l
            end if
            if (abs((upper - lower)/(upper-start_l)) .lt. tol) then 
                cont2 = .false.
            end if
        end do
        omega_end = omega_shift(r_e, r_i, omega)
        call triple_splint(l_phot, k_phot, k2_phot, N, end_l, k_end)
        call place_k_arr((k_start + k_end), (omega_start + omega_end)/2.d0, k_ind)
        !assign path length to grid array
        if (inBounds_arr(cell_ind, k_ind)) then 
            !$OMP CRITICAL(GRID_CRITICAL)
            associate (iR => cell_ind(1), iT => cell_ind(2), &
                iP => cell_ind(3), iKT => k_ind(1), iKP => k_ind(2), &
                iW => k_ind(3))
            cell_time(iR, iT, iP) = cell_time(iR, iT, iP) + (end_l - start_l)
            cell_bins(iR, iT, iP, iKT, iKP, iW) = cell_bins(iR, iT, iP, iKT, iKP, iW) &
                                        + (end_l - start_l)
            end associate 
            !$OMP END CRITICAL(GRID_CRITICAL)
        end if
        !check if photon is above single photon conversion threshhold
        if (single_convert(omega_start, omega_end, k_start, &
                            k_end, cell_ind(1), cell_ind(2), cell_ind(3))) then
            cont = .false. 
            !$OMP CRITICAL(CONVERT_CRITICAL)
            stat_conv(side_ind) = stat_conv(side_ind) + 1
            !$OMP END CRITICAL(CONVERT_CRITICAL)
            return
        end if
        call traj_place(l_phot, coord_phot, coord2_phot, N, upper, indices, r_i)
        !check if path has left the grid or end of path has been reached
        if ((sum(indices) < 0) .or. (abs((ln - upper) / (ln-l0)) < tol)) then
            cont = .false.
        end if
        !if all ending criteria are not met, set up for next cell
        start_l = end_l
        k_start = k_end
        omega_start = omega_end
        cell_ind = indices
    end do
    !check if end point of path is at the surface of the star
    call triple_splint(l_phot, coord_phot, coord2_phot, N, ln, cpt)
    call spherical(cpt(1), cpt(2), cpt(3), pt(1), pt(2), pt(3))
    if (pt(1) .lt. 1.1d0) then
        !check that photon is travelling inward
        call place_surf(pt(2), pt(3), indices(1), indices(2), indices(3))
        call triple_splint(l_phot, k_phot, k2_phot, N, ln, k_end)
        drdL = sum(cpt(:) * k_end(:)) / pt(1)
        !print *, pt, k_end, drdL
        if (drdL .lt. 0.d0) then
            !$OMP CRITICAL(SURFACE_CRITICAL)
            omega_end = omega_shift(r_e, 1.d0, omega)
            surf_flux(indices(1), indices(2), indices(3)) = &
                        surf_flux(indices(1), indices(2), indices(3)) + omega_end
            stat_surf(side_ind) = stat_surf(side_ind) + 1
            !$OMP END CRITICAL(SURFACE_CRITICAL)
        end if
    end if
    return 
    end subroutine grid_prop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! helper function for splining across 3D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine triple_splint(X, Y, Y2, N, X0, Yout)
    !input variables
    integer :: N
    real*8, dimension(2*N_spoints) :: X 
    real*8, dimension(2*N_spoints, 3) :: Y, Y2 
    real*8 :: X0 
    !output variables
    real*8, dimension(3) :: Yout
    !local variables
    integer :: iX 
    do iX = 1, 3
        call splint(X, Y(:, iX), Y2(:, iX), N, X0, Yout(iX), .true.)
    end do
    return
    end subroutine triple_splint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign list of entries to the grid 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_path_to_grid(ind_phot, time_phot, entry_index)
    !input variables
    integer :: ind_phot(N_para_list, 6), entry_index
    real*8 :: time_phot(N_para_list)
    !local variables
    integer :: iX, iR, iT, iP, iKT, iKP, iW
    do iX = 1, (entry_index - 1)
        iR = ind_phot(iX, 1); iT = ind_phot(iX, 2); iP = ind_phot(iX, 3)
        iKT = ind_phot(iX, 4); iKP = ind_phot(iX, 5); iW = ind_phot(iX, 6)
        cell_time(iR, iT, iP) = cell_time(iR, iT, iP) + time_phot(iX)
        cell_bins(iR, iT, iP, iKT, iKP, iW) = cell_bins(iR, iT, iP, iKT, iKP, iW) + time_phot(iX)
    end do
    end subroutine write_path_to_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Assign list of entries to the surface grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine write_path_to_surface(ind_surf, flux_surf, surf_index)
    !input variables
    integer :: ind_surf(N_para_list, 4), surf_index
    real*8 :: flux_surf(N_para_list)
    !local variables
    integer :: iX, iS1, iS2, iS3, iS4 
    do iX = 1, (surf_index - 1)
        iS1 = ind_surf(iX, 1); iS2 = ind_surf(iX, 2); iS3 = ind_surf(iX, 3); iS4 = ind_surf(iX, 4)
        if (iS1 .ne. 0) then 
            surf_flux(iS1, iS2, iS3) = surf_flux(iS1, iS2, iS3) + flux_surf(iX)
            stat_surf(iS4) = stat_surf(iS4) + 1
        else 
            stat_conv(iS4) = stat_conv(iS4) + 1
        end if
    end do
    end subroutine write_path_to_surface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Propagate photons in the grid in parallel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine propagate_parallel(Np)
    implicit none
    !input variable
    integer, intent(in) :: Np !number of photons to propagate on rope
    integer :: ind_phot(Ncells, 3), N_cells, gcalls, pcalls
    real*8 :: time_phot(Ncells)
    !parallelized variables
    integer :: iP, side_ind, Npath
    real*8 :: R, theta, phi, eta, zeta, omega, l0, ln
    real*8, dimension(2*N_spoints) :: l_phot
    real*8, dimension(3) :: dk0, dkn
    real*8, dimension(5) :: coord
    real*8, dimension(2*N_spoints, 3) :: coord_phot, k_phot
    call omp_set_num_threads(N_threads)
    print *, N_threads, " threads used for propagation"
    !$OMP PARALLEL PRIVATE(iP, side_ind, R, theta, phi, eta, zeta, omega) &
    !$OMP& PRIVATE(l_phot, coord_phot, k_phot, dk0, dkn, l0, ln, Npath, coord)
    !$OMP DO
    do iP = 1, Np
        !$OMP CRITICAL
        call emit_point(R, theta, phi, eta, zeta, omega, side_ind)
        !$OMP END CRITICAL
        call traj(R, theta, phi, eta, zeta, l_phot, coord_phot, k_phot, dk0, dkn, l0, ln, Npath)  
        coord = [R, theta, phi, eta, zeta]
        if (ln > l0) then
            call grid_prop(l_phot, coord_phot, k_phot, l0, ln, &
                                            Npath, omega, coord, side_ind)
            !$OMP CRITICAL
            stat_emit(side_ind) = stat_emit(side_ind) + 1
            N_prop = N_prop + 1
            !$OMP END CRITICAL
        end if
        !check 
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    return
    end subroutine propagate_parallel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Propagate photons in the grid in parallel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine propagate_write(Np)
    implicit none
    !input variable
    integer, intent(in) :: Np !number of photons to propagate on rope
    !parallelized variables
    integer :: iP, side_ind, Npath, iThread
    real*8 :: R, theta, phi, eta, zeta, omega, l0, ln
    real*8, dimension(2*N_spoints) :: l_phot
    real*8, dimension(3) :: dk0, dkn
    real*8, dimension(5) :: coord
    real*8, dimension(2*N_spoints, 3) :: coord_phot, k_phot
    !entry variables
    integer :: ind_phot(N_para_list, 6), entry_index, ind_surf(N_para_list, 4), surf_index
    real*8 :: time_phot(N_para_list), flux_surf(N_para_list)
    call omp_set_num_threads(N_threads)
    print *, N_threads, " threads used for propagation"
    !$OMP PARALLEL PRIVATE(iP, side_ind, R, theta, phi, eta, zeta, omega, iThread) &
    !$OMP& PRIVATE(l_phot, coord_phot, k_phot, dk0, dkn, l0, ln, Npath, coord) &
    !$OMP& PRIVATE(ind_phot, entry_index, ind_surf, surf_index, time_phot, flux_surf)
    ind_phot = 0; ind_surf = 0; time_phot = 0.d0; flux_surf = 0.d0
    entry_index = 1; surf_index = 1
    !$OMP DO
    do iP = 1, Np
        iThread = omp_get_thread_num() + 1
        !$OMP CRITICAL(EMIT_CRIT)
        call emit_point(R, theta, phi, eta, zeta, omega, side_ind)
        !$OMP END CRITICAL(EMIT_CRIT)
        call traj(R, theta, phi, eta, zeta, l_phot, coord_phot, k_phot, dk0, dkn, l0, ln, Npath)  
        coord = [R, theta, phi, eta, zeta]
        if ((ln > l0) .and. ((R .gt. 1.001d0) .or. (cos(eta) .gt. 0.1d0))) then
            call path_create(l_phot, coord_phot, k_phot, l0, ln, Npath, omega, coord, side_ind, &
                    ind_phot, time_phot, entry_index, ind_surf, flux_surf, surf_index)
            !$OMP CRITICAL
            stat_emit(side_ind) = stat_emit(side_ind) + 1
            N_prop = N_prop + 1
            !$OMP END CRITICAL
        end if
        if (entry_index .gt. 7*N_para_list/8) then 
            !$OMP CRITICAL
            !print *, "thread", iThread, "ran", entry_index, "grid sections"
            call write_path_to_grid(ind_phot, time_phot, entry_index)
            ind_phot = 0; time_phot = 0.d0; entry_index = 1
            !$OMP END CRITICAL
        end if
        if (surf_index .gt. 7*N_para_list/8) then 
            !$OMP CRITICAL
            !print *, "thread", iThread, "ran", surf_index, "surface sections"
            call write_path_to_surface(ind_surf, flux_surf, surf_index)
            ind_surf = 0; flux_surf = 0.d0; surf_index = 1
            !$OMP END CRITICAL
        end if
    end do
    !$OMP END DO
    !$OMP CRITICAL
    !print *, "thread", iThread, "ran", entry_index, "grid sections"
    call write_path_to_grid(ind_phot, time_phot, entry_index)
    !print *, "thread", iThread, "ran", surf_index, "surface sections"
    call write_path_to_surface(ind_surf, flux_surf, surf_index)
    !$OMP END CRITICAL
    !$OMP END PARALLEL
    return
    end subroutine propagate_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end module parallel_prop
