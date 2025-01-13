!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Module for photon emission and assignment on the grid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module emission
    use parameters_base
    use parameters
    use routines
    use grid
    use spectrum
    use trajectory
    use OMP_LIB
    implicit none 
    real*8, dimension(N_espline) :: emit_Ain_spl, emit_thetain_spl, emit_thetain2_spl
    real*8, dimension(N_espline) :: emit_Aout_spl, emit_thetaout_spl, emit_thetaout2_spl
    real*8, dimension(N_espline) :: emit_Aside_spl, emit_Rside_spl, emit_Rside2_spl
    real*8, dimension(N_espline) :: emit_intside_spl, emit_thetaside_spl, emit_thetaside2_spl
    real*8, save :: area_in, area_left, area_right, area_out, area_sum
    real*8 :: Rm_in, Rm_out
    integer:: idum
    contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine initialize(input)
    real*8 :: theta_step, yp, dphi, A0, Rm_step, theta, Rm, A(1)
    integer :: iT, input
    integer :: draws
    draws = 1d5
    call create_grid
    call create_spline_array
    
    idum = -input
    do iT = 1, draws
        yp = ran2(idum)
    end do
    Rm_out = 1.d0/sin(theta_out)**2.d0
    print *, Rm_out
    dphi = E_phi1 - E_phi0
    Rm_in = 1.d0/sin(theta_in)**2.d0
    print *, Rm_in
    !determine areas of each section
    !inner slice
    A = 0.d0
    call integrate_dlsode(1, A, theta_in, pi-theta_in, deriv_radial)
    area_in = A(1)*Rm_in**2.d0 * dphi
    !area_in = (indef_tan(theta_in, pi-theta_in) - indef_tan(theta_in, theta_in)) * dphi
    !outer slice
    A = 0.d0
    call integrate_dlsode(1, A, theta_out, pi-theta_out, deriv_radial)
    area_out = A(1)*Rm_out**2.d0 * dphi
    !area_out = (indef_tan(theta_out, pi-theta_out) - indef_tan(theta_out, theta_out)) * dphi
     
    area_left = (Rm_out**2.d0) * (intsin4(theta_in) - intsin4(theta_out)) - (theta_in-theta_out) &
                + ((Rm_out**2 - Rm_in**2)/2.d0) * (intsin4(pi - theta_in) - intsin4(theta_in))
    area_right = area_left
    area_sum = area_in + area_out + area_left + area_right
    print *, area_in, area_out, area_left, area_right
    !spline for outer tangent surface
    theta_step = (pi - 2*theta_out) / (N_espline - 1)
    yp = ((Rm_out**2.d0) * sin(theta_out)**(4.d0) * sqrt(3*cos(theta_out)**2 + 1)*dphi)**(-1.d0)
    !set up spline for uniform area emission
    do iT = 1, N_espline
        A = 0.d0
        theta = theta_out + (iT - 1) * theta_step
        emit_thetaout_spl(iT) = theta
        call integrate_dlsode(1, A, theta_out, theta, deriv_radial)
        emit_Aout_spl(iT) = A(1)*Rm_out**2.d0 * dphi
    end do
    call spline(emit_Aout_spl, emit_thetaout_spl, N_espline, yp, yp, emit_thetaout2_spl)
    
    !spline for inner tangent surface
    theta_step = (pi - 2*theta_in) / (N_espline - 1)
    yp =(dphi * (Rm_in**2.d0) * sin(theta_in)**(4.d0) * sqrt(3*cos(theta_in)**2 + 1))**(-1.d0)
    !set up spline for uniform area emission
    do iT = 1, N_espline
        A = 0.d0
        theta = theta_in + (iT - 1) * theta_step
        emit_thetain_spl(iT) = theta
        call integrate_dlsode(1, A, theta_in, theta, deriv_radial)
        emit_Ain_spl(iT) = A(1) * Rm_in**2.d0 * dphi
    end do
    call spline(emit_Ain_spl, emit_thetain_spl, N_espline, yp, yp, emit_thetain2_spl)
    !spline for area to rmax for side surface
    Rm_step = (Rm_out - Rm_in) / (N_espline-1)
    emit_Rside_spl(1) = Rm_in
    emit_Aside_spl(1) = 0.d0
    do iT = 2, N_espline
        Rm = Rm_in + Rm_step * (iT-1)
        emit_Rside_spl(iT) = Rm
        emit_Aside_spl(iT) = emit_Aside_spl(iT-1) + dAsdR(Rm) * Rm_step 
    end do
    call spline(emit_Aside_spl, emit_Rside_spl, N_espline, 1.d0/dAsdR(Rm_in), &
            1.d0/dAsdR(Rm_out), emit_Rside2_spl)
    !spline for rmax to theta for side surface
    theta_step = (pi - 2*theta_out) / (N_espline-1)
    do iT = 1, N_espline
        theta = theta_out + theta_step * (iT-1)
        emit_thetaside_spl(iT) = theta
        emit_intside_spl(iT) = intsin4(theta)
    end do
    call spline(emit_intside_spl, emit_thetaside_spl, N_espline, sin(theta_out)**-4.d0, &
            sin(theta_out)**(-4.d0), emit_thetaside2_spl)
    
    return
    end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function indef_tan(theta0, theta)
    !indefinite integral for the area of tangent surfaces, as a function of theta
    real*8, intent(in) :: theta0, theta    
    real*8 :: indef_tan, Rm
    Rm = 1/sin(theta0)**(2.d0)
    indef_tan = ((Rm**2.d0) / 288.d0) * ((3*sqrt(6*cos(2*theta) + 10))*(3*cos(3*theta)-13*cos(theta)) &
                -52*sqrt(3.d0)*log(sqrt(6.d0)*cos(theta)+sqrt(3*cos(2*theta)+5)))
    return
    end function indef_tan

! integrand for area wrt theta for inner/outer surface
    subroutine deriv_radial(neq, T_in, X, dXdT)
    integer, intent(in) :: neq
    real*8, intent(in) :: T_in, X(neq)
    real*8, intent(out) :: dXdT(neq)
    dXdT = sin(T_in)**4.d0 * (1.d0 + 3.d0*cos(T_in)**2)**(0.5d0)
    return
    end subroutine deriv_radial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine sphere_path(coord_path, k_path, points, sph_path, ksph_path) 
    ! convert the cartesian path from trajectories to a (R, theta, phi) path for grid assignment
    implicit none
    integer, intent(in) :: points !number of points in path
    real*8, intent(in) :: coord_path(N_smax, 3), k_path(N_smax, 3)
    real*8, intent(out) :: sph_path(N_smax, 3), ksph_path(N_smax, 3)
    integer :: iP
    do iP = 1, points
        call spherical(coord_path(iP, 1),coord_path(iP, 2), coord_path(iP, 3), sph_path(iP, 1), sph_path(iP, 2), sph_path(iP, 3))
    end do

    do iP = 1, points
        !dr/dl
        ksph_path(iP, 1) = sum(coord_path(iP, :) * k_path(iP, :)) / sph_path(iP, 1)
        !dtheta/dl
        ksph_path(iP, 2) = ((coord_path(iP, 3) / sph_path(iP, 1)) * ksph_path(iP, 1) - k_path(iP, 3)) &
                            / (sph_path(iP, 1) * sin(sph_path(iP, 2)))
        !dphi/dl
        ksph_path(iP, 3) = (coord_path(iP,1) * k_path(iP, 2) - coord_path(iP, 2) * k_path(iP, 1)) &
                            / (coord_path(iP, 1)**2 + coord_path(iP, 2)**2)
    end do
    end subroutine sphere_path

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tan2sphere(R, theta, phi, alpha, beta, eta, zeta)
    !converts emission from tangent plane to spherical coordinates
    implicit none
    real*8, intent(in) :: R, theta, phi, alpha, beta
    real*8, intent(out) :: eta, zeta
    real*8 :: sin_zeta
    eta = acos((1.d0 + 3.d0 * cos(theta)**2)**(-0.5) * &
                        (cos(alpha)*sin(theta) + 2.d0*sin(alpha)*cos(beta)*cos(theta)))
    if (abs(cos(eta)) == 1) then
        zeta = 0
        return
    end if
    zeta = acos((1.d0 + 3.d0*cos(theta)**2)**(-0.5) * &
                (sin(alpha)*cos(beta)*sin(theta) - 2.d0*cos(alpha)*cos(theta))/sin(eta))
    sin_zeta = sin(alpha) * sin(beta) / sin(eta)
    if (sin_zeta < 0) then
        zeta = twopi - zeta
    end if
    return
    end subroutine tan2sphere


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine tangentnorm(theta, phi, nx, ny, nz)
    !returns the norm vector for the tangent surface to the dipole field
    implicit none
    real*8, intent(in) ::  theta, phi
    real*8, intent(out) :: nx, ny, nz
    nx = 3.d0 * cos(theta) * sin(theta) * cos(phi)
    ny = 3.d0 * cos(theta) * sin(theta) * sin(phi)
    nz = 2.d0 * cos(theta)**2 - sin(theta)**2
    return
    end subroutine tangentnorm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine emit_point(R, theta, phi, eta, zeta, omega, side_ind)
    implicit none
    !output variables
    real*8, intent(out) :: R, theta, phi, eta, zeta, omega
    integer, intent(out) :: side_ind
    !local variables
    real*8 :: area
    logical :: cont 
    cont = .true. 
    do while (cont)
        area = area_sum*ran2(idum)
        if (area .lt. area_in+area_left) then
            if (area .lt. area_in) then
                call emit_point_in(R, theta, phi, eta, zeta)
                side_ind = 1
            else
                call emit_point_left(R, theta, phi, eta, zeta)
                side_ind = 2
            end if
        else
            area = area - (area_in+area_left)
            if (area .lt. area_out) then
                call emit_point_out(R, theta, phi, eta, zeta)
                side_ind = 3
            else 
                call emit_point_right(R, theta, phi, eta, zeta)
                side_ind = 4
            end if
        end if
        if (valid_emission(R, theta, phi, eta, zeta)) then
            cont = .false.
        end if
    end do
    call splint(spec_num_spl, spec_omega_spl, spec_omega2_spl, &
                    N_ospline, 1.d0*ran2(idum), omega, .true.)
    end subroutine emit_point

    function valid_emission(R, theta, phi, eta, zeta) result(check)
    !input functions
    real*8 :: R, theta, phi, eta, zeta 
    !output
    logical :: check 
    !local variables
    real*8 :: coord(5)
    coord = (/R, theta, phi, eta, zeta/)
    check = .true.
    if (any(isnan(coord(:)))) then
        check = .false.
        return
    end if
    if ((R .lt. 1.001d0) .and. (cos(eta) .lt. 0.d0)) then 
        check = .false.
    end if
    return 
    end function valid_emission
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine emit_point_out(R, theta, phi, eta, zeta)
    ! emits a photon from the flux surface isotropically
    ! assumes dipole magnetic field
    ! flux surface covers the first slice in phi, with a square cross section.
    ! emission for outer side
    implicit none
    real*8, intent(out) :: R, theta, phi, eta, zeta
    real*8 :: alpha, beta, path_A
    path_A = area_out * ran2(idum)
    call splint(emit_Aout_spl, emit_thetaout_spl, emit_thetaout2_spl,&
                    N_espline, path_A, theta, .true.)
    R = (sin(theta) / sin(theta_out))**2
    alpha = acos(ran2(idum)) !outward emission
    beta = twopi * ran2(idum)
    phi = ran2(idum) * (E_phi1 - E_phi0) + E_phi0 
    phi = mod(phi + twopi, twopi)
    call tan2sphere(R, theta, phi, alpha, beta, eta, zeta)
    return
    end subroutine emit_point_out
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine emit_point_in(R, theta, phi, eta, zeta)
    ! emits a photon from the flux surface isotropically
    ! assumes dipole magnetic field
    ! flux surface covers the first slice in phi, with a square cross section.
    ! emission for inner side
    implicit none
    real*8, intent(out) :: R, theta, phi, eta, zeta
    real*8 :: alpha, beta, path_A
    path_A = area_in * ran2(idum)
    call splint(emit_Ain_spl, emit_thetain_spl, emit_thetain2_spl, N_espline, path_A, theta, .true.)
    R = (sin(theta) / sin(theta_in))**2
    alpha = acos(-ran2(idum)) !inward emission
    beta = twopi * ran2(idum)
    phi = ran2(idum) * (E_phi1 - E_phi0) + E_phi0
    phi = mod(phi + twopi, twopi)
    call tan2sphere(R, theta, phi, alpha, beta, eta, zeta)

    return
    end subroutine emit_point_in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function intsin4(x)
    !indefinite integral for area along a given R_max
    real*8 :: intsin4
    real*8, intent(in) :: x
    intsin4 = (12.d0 * x - 8*sin(2*x) + sin(4*x)) / 32
    return
    end function intsin4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function dAsdR(Rm)
    !derivative of side area w.r.t Rmax
    real*8  :: dAsdR, tmin, tmax
    real*8, intent(in) :: Rm
    tmin = asin(sqrt(1.d0/Rm))
    tmax = pi - tmin
    dAsdR =  Rm * (intsin4(tmax) - intsin4(tmin))
    return
    end function dAsdR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine emit_point_left(R, theta, phi, eta, zeta)
    ! emits a photon from the flux surface isotropically
    ! assumes dipole magnetic field
    ! flux surface covers the first slice in phi, with a square cross section.
    ! left side emission (from top down view)
    implicit none
    real*8, intent(out) :: R, theta, phi, eta, zeta
    real*8 :: alpha, beta, Rm, area, sinmin, sinmax, intsin, sinzeta
    area = ran2(idum) * area_left
    call splint(emit_Aside_spl, emit_Rside_spl, emit_Rside2_spl, N_espline, area, Rm, .true.)
    sinmin = intsin4(asin(sqrt(1/Rm)))
    sinmax = intsin4(pi - asin(sqrt(1/Rm)))
    intsin = sinmin + (sinmax - sinmin) * ran2(idum)
    call splint(emit_intside_spl, emit_thetaside_spl, emit_thetaside2_spl, &
                N_espline, intsin, theta, .true.)
    R = Rm * (sin(theta))**2.d0
    alpha = acos(ran2(idum)) !outward emission
    beta = twopi * ran2(idum)
    phi = mod(E_phi1+twopi, twopi)
    eta = acos(-sin(alpha) * sin(beta))
    if (eta .eq. 0) then
        zeta = 0
    else
        zeta = acos(sin(alpha) * cos(beta) / sin(eta))
    end if
    
    return
    end subroutine emit_point_left

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine emit_point_right(R, theta, phi, eta, zeta)
    ! emits a photon from the flux surface isotropically
    ! assumes dipole magnetic field
    ! flux surface covers the first slice in phi, with a square cross section.
    ! right side emission (from top down view)
    implicit none
    real*8, intent(out) :: R, theta, phi, eta, zeta
    real*8 :: alpha, beta, Rm, area, sinmin, sinmax, intsin
    area = ran2(idum) * area_right
    call splint(emit_Aside_spl, emit_Rside_spl, emit_Rside2_spl, N_espline, area, Rm, .true.)
    sinmin = intsin4(asin(sqrt(1/Rm)))
    sinmax = intsin4(pi - asin(sqrt(1/Rm)))
    intsin = sinmin + (sinmax - sinmin) * ran2(idum)
    call splint(emit_intside_spl, emit_thetaside_spl, emit_thetaside2_spl, &
           N_espline, intsin, theta, .true.)
    R = Rm * (sin(theta))**2.d0
    alpha = acos(ran2(idum)) !outward emission
    beta = twopi * ran2(idum)
    phi = mod(E_phi0+twopi, twopi)
    eta = acos(sin(alpha) * sin(beta))
    if (eta .eq. 0) then
        zeta = 0
    else
        zeta = acos(sin(alpha) * cos(beta) / sin(eta))
        zeta = twopi - zeta
    end if
    return 
    end subroutine emit_point_right
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Find redshifted frequency as a function of radius
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function omega_shift(r_e, r_i, omega_e)
    real*8, intent(in) :: r_e, r_i, omega_e
    real*8 :: omega_shift

    omega_shift = omega_e * sqrt(1 - (rg / r_e)) / sqrt(1 - (rg / r_i))


    end function omega_shift



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Given a path length and path, determine cell indices, radius
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine traj_place(l_phot, coord_phot, coord2_phot, Npath, l, indices, r_l)
    implicit none
    real*8, intent(in) :: l !path length
    ! spline arrays
    integer, intent(in) :: Npath !number of valid points in the array
    real*8, dimension(2*N_spoints, 3), intent(in) :: coord_phot, coord2_phot
    real*8, dimension(2*N_spoints), intent(in) :: l_phot
    ! output variables
    integer, dimension(3), intent(out) :: indices !spatial cell indices
    real*8, intent(out) :: r_l !radius at path length l
    ! intermediate variables
    real*8, dimension(3) :: cpt !point in path in cartesian
    real*8, dimension(3) :: pt !point in path in spherical
    integer :: iX 
    do iX = 1, 3
        call splint(l_phot, coord_phot(:, iX), coord2_phot(:, iX), Npath, l, cpt(iX), .true.)
    end do
    call spherical(cpt(1), cpt(2), cpt(3), pt(1), pt(2), pt(3))
    r_l = pt(1)
    call place(pt(1), pt(2), pt(3), indices(1), indices(2), indices(3))
    end subroutine traj_place


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! helper function to determine if photon is past single photon threshold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function single_convert(w_start, w_end, k_start, k_end, indR, indT, indP) result(yn)
    !output variable
    logical :: yn
    !input variables
    real*8 :: w_start, w_end, k_start(3), k_end(3)
    integer :: indR, indT, indP
    !local variables
    real*8:: omega, sin_kb
    real*8, dimension(3) :: k, B, Bhat

    omega = (w_start + w_end)/2.d0
    k = (k_start + k_end)/2.d0
    call cellB(indR, indT, indP, B)
    Bhat = B/norm2(B)
    sin_kb = (1.d0 - sum(k*Bhat)**2)**(0.5d0)
    yn = (omega * sin_kb) .ge. (2 * m_e)
    return
    end function single_convert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine grid_assign(l_phot, coord_phot, k_phot, dk0, dkn, &
                    l0, ln, Np, coord, Npath, omega, side_ind)
    !assign a photon path onto the grid using bisection
    implicit none
    !input variables
    ! start/end path length, start/end photon curvature dk/dl, initial emission energy
    real*8, intent(in) :: l0, ln, dk0(3), dkn(3), omega
    ! spline arrays
    real*8, dimension(2*N_spoints) :: l_phot
    real*8, dimension(2*N_spoints, 3), intent(in) :: k_phot, coord_phot
    integer :: side_ind
    integer, intent(in) :: Np, Npath
    real*8, dimension(5) :: coord
    !local variables
    real*8, dimension(2*N_spoints, 3) :: coord2_phot, k2_phot
    integer :: N
    real*8 :: start_l, end_l, lower, upper, mid, tol, cpt(3), pt(3), scale_l
    real*8, dimension(3) :: k_start, k_end, n_hat
    real*8 :: w_start, w_end, r_e, r_i, T01, T02, P1, P2, mu
    integer :: iR, iP, indices(3), calls, cell_ind(3), indR, indT, indP, k_ind(3)
    logical :: cont, cont2
    character(len=8) :: num ! format descriptor
    write (num, "(I0)") Np
    tol = 0.001
    N = Npath
    do iR = 1, 3
        call spline(l_phot, coord_phot(:, iR), N, k_phot(1, iR), k_phot(N, iR), coord2_phot(:, iR))
        call spline(l_phot, k_phot(:, iR), N, dk0(iR), dkn(iR), k2_phot(:, iR))
    end do
    start_l = l0
    w_start = omega
    call traj_place(l_phot, coord_phot, coord2_phot, Npath, l0, indices, r_e)
    cont2 = .true.
    do iR = 1, 3
        call splint(l_phot(1:N), k_phot(1:N, iR), k2_phot(1:N, iR), N, start_l, k_start(iR), .true.)
    end do
    
101 format(*(g0, :, ", "))
    do while (cont2)
        calls = 0
        cont = .true.
        cell_ind(:) = indices(:) !hold index of current cell
        !print *, cell_ind
        if (sum(cell_ind) > 0) then
            scale_l = cell_length(cell_ind(1), cell_ind(2), cell_ind(3))
            lower = start_l
            if (lower + scale_l > ln) then
                cont = .false.
                upper = ln
            else
                upper = lower + scale_l
            end if
        else
            cont = .false.
            lower = start_l
            upper = ln
        end if
        !shift bounds on endpoint by characteristic cell length until new cell is found 
        do while (cont)
            call traj_place(l_phot, coord_phot, coord2_phot, Npath, upper, indices, r_i)
            calls = calls + 1
            if (calls .ge. 1d6) then
                cont = .false.
                cont2 = .false.    
                call file_open(17, "error.txt")
                write (17, *) "part 1 too long"
                write (17, 101) coord
                write (17, *)  Np
                close(17)
                return
            end if
            if ((all(cell_ind(:) == indices(:))) .and. (abs((ln - upper) / ln) .ge. tol)) then
                lower = upper
                upper = min(lower + scale_l, ln)
            else
                cont = .false.
            end if
        end do
        cont = .true.
        calls = 0
        do while (cont)
            end_l = (upper + lower) / 2.d0
            call traj_place(l_phot, coord_phot, coord2_phot, Npath, end_l, indices, r_i)
            calls = calls + 1
            if (all(cell_ind(:) == indices(:))) then
                lower = end_l    
            else
                upper = end_l
            end if
            if (abs((upper - lower) / (upper - start_l)) < tol) then
                cont = .false.
            end if
            if (calls .ge. 1d6) then
                cont = .false.
                cont2 = .false.    
                call file_open(17, "error.txt")
                write (17, *) "part 2 too long"
                write (17, 101) coord
                write (17, *)  Np
                write (17, *) cell_ind
                close (17)
            end if
        end do
        indR = cell_ind(1); indT = cell_ind(2); indP = cell_ind(3)
        call traj_place(l_phot, coord_phot, coord2_phot, Npath, upper, indices, r_i)
        calls = calls + 1
        w_end = omega_shift(r_e, r_i, omega)
        do iR = 1, 3
            call splint(l_phot(1:N), k_phot(1:N, iR), k2_phot(1:N, iR), N, upper, k_end(iR), .true.)
        end do
        call place_k((k_start + k_end), (w_start + w_end)/2.d0, k_ind(1), k_ind(2), k_ind(3))
        if (inBounds(indR, indT, indP, k_ind(1), k_ind(2))) then
            cell_time(indR, indT, indP) = cell_time(indR, indT, indP) + end_l - start_l
            associate (bin => cell_bins(indR, indT, indP, k_ind(1), k_ind(2), k_ind(3)))
            bin = bin + end_l - start_l 
            end associate
            !check if photon is above single photon conversion threshhold
            if (single_convert(w_start, w_end, k_start, k_end, indR, indT, indP)) then
                cont2 = .false. 
                stat_conv(side_ind) = stat_conv(side_ind) + 1
            return
            end if
        else
            call file_open(17, "error.txt")
            write (17, *) "placement wrong"
            write (17, *) indR, indT, indP, k_ind(1), k_ind(2)
            write (17, 101) coord
            write (17, *)  Np
            close (17)
        end if
        !print *, cell_ind, end_l - start_l
        if ((sum(indices) < 0) .or. (abs((ln - upper) / (ln-l0)) < tol)) then
            cont2 = .false.
        end if
        start_l = upper
        k_start = k_end
        w_start = w_end
    end do
    ! determine if photon hit surface of the star     
    do iR = 1, 3
        call splint(l_phot(1:N), coord_phot(1:N, iR), coord2_phot(1:N, iR), N, ln, cpt(iR), .true.)
        call splint(l_phot(1:N), k_phot(1:N, iR), k2_phot(1:N, iR), N, ln, k_end(iR), .true.)
    end do
    call spherical(cpt(1), cpt(2), cpt(3), pt(1), pt(2), pt(3))
    if (pt(1) .le. 1.1d0) then
        call place_surf(pt(2), pt(3), indices(1), indices(2), indices(3))
        if (indices(1) .lt. Nr+1) then
            T01 = theta01(indices(1))
            T02 = theta02(indices(1))
            P1 = phi1(indices(3))
            P2 = phi2(indices(3))
            call cartesian(1.d0, (T01 + T02)/2.d0, (P1 + P2)/2.d0, n_hat(1), n_hat(2), n_hat(3))
        else
            T01 = theta02(Nr); T02 = 0.d0; P1 = 0.d0; P2 = 0.d0;
            n_hat = 0.d0
            n_hat(3) = 1.d0
        end if
        if (indices(2) .eq. 2) then
            n_hat(3) = -n_hat(3)
        end if
        mu = sum(n_hat*k_end)
        if (mu .lt. 0) then
            associate (bin => surf_flux(indices(1), indices(2), indices(3)))
            w_end = omega_shift(r_e, 1.d0, omega)
            bin = bin + w_end
            stat_surf(side_ind) = stat_surf(side_ind) + 1
            end associate
        end if
    end if
    return
    end subroutine grid_assign
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine path_assign(l_phot, coord_phot, k_phot, dk0, dkn, l0, ln, &
            Npath, ind_phot, time_phot, kind_phot, N_cells, omega)
    !assign a photon path onto the grid using bisection
    implicit none
    real*8, intent(in) :: l0, ln, dk0(3), dkn(3)
    real*8, dimension(2*N_spoints) :: l_phot
    real*8, dimension(2*N_spoints, 3), intent(in) :: k_phot, coord_phot
    real*8, dimension(2*N_spoints, 3) :: coord2_phot, k2_phot
    integer, intent(in) :: Npath
    integer :: N
    real*8 :: start_l, end_l, lower, upper, mid, tol, cpt(3), pt(3), scale_l, omega
    real*8, dimension(3) :: k_start, k_end
    integer :: iR, iP, indices(3), calls, k_ind(3)
    logical :: cont, cont2
    integer, intent(out) :: ind_phot(Ncells, 3), N_cells, kind_phot(Ncells, 3)
    real*8, intent(out) :: time_phot(Ncells)
    
    tol = 0.001
    N = Npath
    do iR = 1, 3
        call spline(l_phot, coord_phot(:, iR), N, k_phot(1, iR), k_phot(N, iR), coord2_phot(:, iR))
        call spline(l_phot, k_phot(:, iR), N, dk0(iR), dkn(iR), k2_phot(:, iR))
        call splint(l_phot, coord_phot(:, iR), coord2_phot(:, iR), N, l0, cpt(iR), .true.)
    end do
    start_l = l0
    call spherical(cpt(1), cpt(2), cpt(3), pt(1), pt(2), pt(3))
    call place(pt(1), pt(2), pt(3), indices(1), indices(2), indices(3))
    calls = 0
    cont2 = .true.
    do iR = 1, 3
        call splint(l_phot(1:N), k_phot(1:N, iR), k2_phot(1:N, iR), N, start_l, k_start(iR), .true.)
    end do
    iP = 1
101 format(*(g0, :, ", "))
    do while (cont2)
        cont = .true.
        ind_phot(iP, :) = indices(:)
        if (sum(indices) > 0) then
            scale_l = 0.5d0*cell_length(ind_phot(iP, 1), ind_phot(iP, 2), ind_phot(iP, 3))
            lower = start_l
            if (lower + scale_l > ln) then
                cont = .false.
                upper = ln
            else
                upper = lower + scale_l
            end if
        else
            cont = .false.
        end if
        do while (cont)
            do iR = 1, 3
                call splint(l_phot(1:N), coord_phot(1:N, iR), coord2_phot(1:N, iR), N, upper, cpt(iR), .true.)
            end do
            calls = calls + 1
            call spherical(cpt(1), cpt(2), cpt(3), pt(1), pt(2), pt(3))
            call place(pt(1), pt(2), pt(3), indices(1), indices(2), indices(3))
            if (calls .ge. 1000) then
                cont = .false.
                cont2 = .false.    
                print *, "part 1 too long"
            end if
            if ((all(ind_phot(iP, :) == indices(:))) .and. ((ln - upper) / ln .ge. tol)) then
                lower = upper
                upper = min(lower + scale_l, ln)
            else
                cont = .false.
            end if
        end do
        cont = .true.
        do while (cont)
            end_l = (upper + lower) / 2.d0
            do iR = 1, 3
                call splint(l_phot(1:N), coord_phot(1:N, iR), coord2_phot(1:N, iR), N, end_l, cpt(iR), .true.)
            end do
            call spherical(cpt(1), cpt(2), cpt(3), pt(1), pt(2), pt(3))
            call place(pt(1), pt(2), pt(3), indices(1), indices(2), indices(3))
            calls = calls + 1
            if (all(ind_phot(iP, :) == indices(:))) then
                lower = end_l    
            else
                upper = end_l
            end if
            if ((upper - lower) / upper < tol) then
                cont = .false.
            end if
            if (calls .ge. 1000) then
                cont = .false.
                cont2 = .false.    
                print *,  "part 2 too long"
            end if
        end do
        do iR = 1, 3
            call splint(l_phot(1:N), coord_phot(1:N, iR), coord2_phot(1:N, iR), N, upper, cpt(iR), .true.)
            call splint(l_phot(1:N), k_phot(1:N, iR), k2_phot(1:N, iR), N, upper, k_end(iR), .true.)
        end do
        time_phot(iP) = end_l - start_l
        call place_k(k_end + k_start, omega, k_ind(1), k_ind(2), k_ind(3))
        kind_phot(iP, :) = k_ind(:)
        call spherical(cpt(1), cpt(2), cpt(3), pt(1), pt(2), pt(3))
        call place(pt(1), pt(2), pt(3), indices(1), indices(2), indices(3))
        calls = calls + 1
        if ((sum(indices) < 0) .or. ((ln - upper) / ln < tol)) then
            cont2 = .false.
        end if
        start_l = upper
        k_start = k_end
        iP = iP + 1
    end do
    N_cells = iP - 1
    return
    end subroutine path_assign
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine propagate(Np)
    implicit none
    integer, intent(in) :: Np !number of photons to propagate on rope
    integer :: iP, point_max, Npath, side_ind
    real*8 :: R, theta, phi, eta, zeta, omega, ln, l0
    real*8, dimension(2*N_spoints) :: l_phot
    real*8, dimension(3) :: dk0, dkn
    real*8, dimension(5) :: coord
    real*8, dimension(2*N_spoints, 3) :: coord_phot, k_phot
    integer :: ind_phot(Ncells, 3), N_cells, gcalls, pcalls
    real*8 :: time_phot(Ncells)
    do iP = 1, Np
        call emit_point(R, theta, phi, eta, zeta, omega, side_ind)
        call traj(R, theta, phi, eta, zeta, l_phot, coord_phot, k_phot, dk0, dkn, l0, ln, Npath)    
        coord = [R, theta, phi, eta, zeta]
        if (ln > l0) then
            call grid_assign(l_phot, coord_phot, k_phot, dk0, dkn, &
                    l0, ln, iP, coord, Npath, omega, side_ind)
            stat_emit(side_ind) = stat_emit(side_ind) + 1
            N_prop = N_prop + 1
        end if
    end do
    return
    end subroutine propagate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine output_emit_param(runnum)
    implicit none
    integer, intent(in) :: runnum
    character(len=8) :: num ! format descriptor
    write (num, "(I0)") runnum
    call file_open(2, trim(num)//'_emit_param.dat')
101 format(*(g0, :, ", "))
    write (2, 101) theta_out, theta_in, E_phi0, E_phi1
    close (2)
    end subroutine output_emit_param

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine emission_test(runnum, pnum)
    implicit none
    real*8 :: area
    real*8 :: R, theta, phi, eta, zeta
    integer, intent(in) :: runnum, pnum
    integer :: iP
    character(len=8) :: num ! format descriptor
    write (num, "(I0)") runnum
    call file_open(11, trim(num)//'_emit_left.dat')
    call file_open(12, trim(num)//'_emit_right.dat')
    call file_open(13, trim(num)//'_emit_in.dat')
    call file_open(14, trim(num)//'_emit_out.dat')
101 format(*(g0, :, ", "))
    do iP = 1, pnum
        area = area_sum*ran2(idum)
        if (area .lt. area_right+area_left) then
            if (area .lt. area_right) then
                call emit_point_right(R, theta, phi, eta, zeta)
                write(12, 101) R, theta, phi, eta, zeta
                stat_emit(4) = stat_emit(4) + 1
            else
                call emit_point_left(R, theta, phi, eta, zeta)
                write(11, 101) R, theta, phi, eta, zeta
                stat_emit(2) = stat_emit(2) + 1
            end if
        else
            area = area - (2.d0*area_left)
            if (area .lt. area_out) then
                call emit_point_out(R, theta, phi, eta, zeta)
                write(14, 101) R, theta, phi, eta, zeta
                stat_emit(3) = stat_emit(3) + 1
            else
                call emit_point_in(R, theta, phi, eta, zeta)
                write(13, 101) R, theta, phi, eta, zeta
                stat_emit(1) = stat_emit(1) + 1
            end if
        end if
    end do
    close(11); close(12); close(13); close(14)
    end subroutine emission_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end module emission
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

