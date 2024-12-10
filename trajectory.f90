!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Module for generating photon trajectories around the neutron star
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module trajectory
    use parameters_base
    use routines
    implicit none

    real*8, dimension(N_spline, N_spoints) :: r_data, psi_data, l_data, rdpsi_data, dr_data
    real*8, dimension(N_spline, N_spoints) :: r2_data, psi2_data, l2_data, rdpsi2_data, dr2_data
    real*8, dimension(N_spline) :: b_data, Rmin_data
    real*8 :: b_var
    contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function b_imp(R, eta)
    implicit none
    real*8, intent(in) :: R, eta
    real*8 :: b_imp
    if (R < rg) then
        b_imp = 0
        return
    end if
    b_imp = R * sin(eta) * (1 - (rg/R))**(-0.5)

    return
    end function b_imp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function R_min(b)
    implicit none
    real*8, intent(in) :: b
    real*8 :: R_min
    if (b < b_inner) then
        R_min = R_inner
    else
        R_min = 2 * (b**2 / 3)**(0.5d0) * cos(acos(-3*rg/2 * (3/b**2)**0.5d0) / 3.d0)
    end if
    return 
    end function R_min

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function dpsi_dr(r)
    implicit none
    real*8, intent(in) :: r
    real*8  :: dpsi_dr
    dpsi_dr = (r**(-2.d0)) * ((b_var**(-2.d0)) - (r**(-2.d0))*(1.d0 - (rg/r)))**(-0.5d0)
    return 
    end function dpsi_dr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function dl_dr(r)
    implicit none
    real*8, intent(in) :: r
    real*8 :: dl_dr
    
    dl_dr = (1.d0 + (r * dpsi_dr(r))**2)**(0.5)

    return
    end function dl_dr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine derivs_r(neq, r_in, X, dXdr)
    implicit none
    integer, intent(in) :: neq
    real*8, intent(in) :: r_in, X(neq)
    real*8, intent(out) :: dXdr(neq)
    !print *, r_in
    dXdr(1) = dpsi_dr(r_in)
    dXdr(2) = dl_dr(r_in)
    !print *, dXdr
    return
    end subroutine derivs_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    function f(r)
    implicit none
    real*8, intent(in) :: r
    real*8 :: f
    f = ( (r / b_var)**(2.d0) - 1 + (rg/r)) ** (-0.5d0)
    end function f


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function dfdr(r)
    implicit none
    real*8, intent(in) :: r
    real*8 :: dfdr

    dfdr = f(r)**(3.d0) * ((rg / (2*r**2.d0)) - (r / b_var**2.d0))

    end function dfdr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine derivs_l(neq, l_in, X, dXdl)
    implicit none
    integer, intent(in) :: neq
    real*8, intent(in) :: l_in, X(neq)
    real*8, intent(out) :: dXdl(neq)
    real*8 :: rdpsi
    ! X(1) = psi, X(2) = r
    rdpsi = f(X(2)) !r dpsi/dr
    dXdl(1) = (rdpsi / X(2)) * (1.d0 + rdpsi**(2.d0))**(-0.5d0)! dpsi/dl
    dXdl(2) = (1.d0 + rdpsi**(2.d0))**(-0.5d0)! dr/dl
    
    return
    end subroutine derivs_l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine spline_data(b, psi_arr, r_arr, l_arr, rdpsi_arr, dr_arr)
    implicit none
    real*8, intent(in) :: b
    real*8, dimension(N_spoints), intent(out) :: psi_arr, r_arr, l_arr, rdpsi_arr, dr_arr
    real*8 :: logr_stop, logr_start, log_step, r1, r2
    real*8 :: X(2), l_max, r_stop, r_start, l_step, dXdl(2)
    integer :: iP
    logr_stop = log10(R_limit)
    r_stop = R_limit
    b_var = b
    if (b == 0) then
        l_step = (r_stop - R_inner) / (N_spoints - 1)
        do iP = 1, N_spoints
            r_arr(iP) = R_inner + l_step * (iP - 1)
            l_arr(iP) = l_step*(iP - 1)
            psi_arr(iP) = pi/2.d0
            rdpsi_arr(iP) = 0.d0
            dr_arr(iP) = 1.d0
        end do

    else
        if (R_min(b) < R_inner) then
            r_start = R_inner  !start spline from half of stellar radius
            X(:) = 0
            call integrate_dlsode(2, X, r_start, r_stop, derivs_r)
            l_max = X(2)
            l_step = l_max / (N_spoints - 1)
            r_arr(1) = r_start
            
            psi_arr(1) = (pi/2.d0) - asin(b*((1 - rg/R_inner)**0.5d0)/R_inner) ; l_arr(1) = 0.d0
            dr_arr(1) = 0; rdpsi_arr(1) = 1.d0
            X(1) = psi_arr(1); X(2) = r_arr(1)
            do iP = 2, N_spoints
                l_arr(iP) = l_arr(1) + l_step * (iP - 1)
                call integrate_dlsode(2, X, l_arr(iP - 1), l_arr(iP), derivs_l)
                psi_arr(iP) = X(1)
                r_arr(iP) = X(2)
                call derivs_l(2, l_arr(iP), X, dXdl)
                dr_arr(iP) = dXdl(2)
                rdpsi_arr(iP) = r_arr(iP) * dXdl(1)
            end do
        else
            r_start = R_min(b)
            r_arr(1) = r_start * (1.d0 + dpsi**2)**0.5
            l_arr(1) = dpsi*r_start
            psi_arr(1) = dpsi
            dr_arr(1) = dl_dr(r_arr(1))**(-1.d0)
            rdpsi_arr(1) = r_arr(1) * dpsi_dr(r_arr(1)) * dr_arr(1)
            X(:) = 0
            call integrate_dlsode(2, X, r_arr(1), R_limit, derivs_r)
            l_max = X(2)

            l_step = (l_max) / (N_spoints - 1)
            X(1) = psi_arr(1); X(2) = r_arr(1)
        
            do iP = 2, N_spoints
                l_arr(iP) = l_arr(1) + l_step * (iP - 1)
                call integrate_dlsode(2, X, l_arr(iP - 1), l_arr(iP), derivs_l)
                psi_arr(iP) = X(1)
                r_arr(iP) = X(2)
                call derivs_l(2, l_arr(iP), X, dXdl)
                dr_arr(iP) = dXdl(2)
                rdpsi_arr(iP) = r_arr(iP) * dXdl(1)
            end do
        end if

    end if
    return
    end subroutine spline_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine create_spline_array
    implicit none
    integer :: iB, n_crit
    real*8 :: logb_step
    real*8 :: drdl2(2), drpsidl2(2), rdpsi(2), df(2)
    b_data(1) = 0.d0
    call spline_data(b_data(1), psi_data(1, :), r_data(1, :), l_data(1, :), rdpsi_data(1, :), dr_data(1, :))
    psi2_data(1, :) = 0; r2_data(1, :) = 0; l2_data(1, :) = 0; rdpsi2_data(1, :) = 0; dr2_data(1, :) = 0
    !calculate steps in b such that it goes from b_min, includes b_inner, and stops slightly past b_max
    n_crit = floor((N_spline-1) * (log10(b_inner) - logb_min) / (log10(b_max) - logb_min))
    logb_step = (log10(b_inner) - logb_min) / n_crit
    !print *, logb_step
    Rmin_data(1) = R_inner
    do iB = 2, N_spline
        b_data(iB) = 1.d1 ** (logb_min + logb_step * (iB - 1))
        Rmin_data(iB) = max(R_inner, R_min(b_data(iB)))
        call spline_data(b_data(iB), psi_data(iB, :), r_data(iB, :), l_data(iB, :), rdpsi_data(iB, :), dr_data(iB, :))
    end do
    !print *, b_data(N_spline) / b_max !should be slightly more than 1
    do iB = 2, N_spline
        b_var = b_data(iB)
        !spline for r(l)
        call spline(l_data(iB, :), r_data(iB, :), N_spoints, dr_data(iB, 1), dr_data(iB, N_spoints), r2_data(iB, :))
        !spline for psi(l)
        call spline(l_data(iB, :), psi_data(iB, :), N_spoints, rdpsi_data(iB, 1)/r_data(iB, 1), &
                        rdpsi_data(iB, N_spoints)/r_data(iB, N_spoints), psi2_data(iB, :))
        !spline for dr(l)/dl
        rdpsi(1) = f(r_data(iB, 1)); rdpsi(2) = f(r_data(iB, N_spoints))
        df(1) = dfdr(r_data(iB, 1)); df(2) = dfdr(r_data(iB, N_spoints))
        drdl2(:) = -rdpsi(:) * (1.d0 + rdpsi(:)**2.d0)**(-2.d0) * df(:)
        call spline(l_data(iB, :), dr_data(iB, :), N_spoints, drdl2(1), drdl2(2), dr2_data(iB, :))
        !spline for rdpsi(l)/dl
        drpsidl2(:) = df(:) * (1 + rdpsi(:)**2.d0)**(-2.d0)
        call spline(l_data(iB, :), rdpsi_data(iB, :), N_spoints, drpsidl2(1), drpsidl2(2), rdpsi2_data(iB, :))
        !spline for l(r)
        call spline(r_data(iB, :), l_data(iB, :), N_spoints, dl_dr(r_data(iB, 1)), dl_dr(r_data(iB, N_spoints)), l2_data(iB, :)) 
    end do
    end subroutine create_spline_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine b_ratio(b, c1, c2, i1, i2)
    implicit none
    real*8, intent(in) :: b
    logical :: arr_increasing
    real*8, intent(out) :: c1, c2
    real*8 :: b1, b2, r, r1, r2
    integer, intent(out) :: i1, i2
    
    arr_increasing = .true.    
    i1 = bisect(b_data, b, n_spline, arr_increasing)
    if (i1 == n_spline) then
        c1 = 1; c2 = 0; i2 = 0; !impact parameter is above b_max
        return
    end if

    i2 = i1 + 1
    b1 = b_data(i1); b2 = b_data(i2)
    r = R_min(b)
    r1 = R_min(b1) ; r2 = R_min(b2)
    if (b .le. b_inner) then
        c1 = 1.d0 * (b2 - b) / (b2 - b1)
        c2 = 1.d0 * (b - b1) / (b2 - b1)
    else
        c1 = 1.d0 * (r2 - r) / (r2 - r1)
        c2 = 1.d0 * (r - r1) / (r2 - r1)
    end if
    return
    end subroutine b_ratio

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine traj(R, theta, phi, eta, zeta, l_phot, coord_phot, k_phot, dk0, dkn, l0, ln, Npath)
    !determine trajectory of a photon emitted at R, theta, phi
    !   with emisison angles eta w.r.t radial vector, zeta w.r.t theta vector
    !trajectory in cartesian
    implicit none
    real*8, intent(in) :: R, theta, phi, eta, zeta
    real*8, dimension(N_spoints*2), intent(out) :: l_phot !points of path length of photon
    real*8, dimension(N_spoints*2, 3), intent(out) :: coord_phot, k_phot !coordinate and direction of photon, in cartesian  
    real*8, dimension(3), intent(out) :: dk0, dkn !dk/dl at start and end point, for splining
    integer, intent(out) :: Npath
    real*8, intent(out) :: l0, ln !starting l and ending l for trajectory
    real*8 :: rm, bs(2), psi0!initial values for linearly interpolated spline
    integer :: iL, iX, iBs(2), iK, i1, i2
    real*8 :: step, b, r_ratio, c1, c2, l, radius, l_stop, l_mid, rad, psi1, psi2, rdpsi
    real*8, dimension(2) :: fdata, df, drdl2, dpsidl2
    real*8, dimension(3) :: r1_hat, psi1_hat, r2_hat, psi2_hat, r0_hat, n_hat, psi0_hat, k_hat
    real*8, dimension(N_spoints) ::  r_spl, psi_spl, l_spl, dr_spl, rdpsi_spl, &
                                    r_spl2, psi_spl2, l_spl2, dr_spl2, rdpsi_spl2
    logical :: cont, abscissa_increasing


    !initial emission direction
    k_hat(1) = cos(eta)*sin(theta)*cos(phi) + sin(eta)*(cos(zeta)*cos(theta)*cos(phi) - sin(zeta)*sin(phi))
    k_hat(2) = cos(eta)*sin(theta)*sin(phi) + sin(eta)*(cos(zeta)*cos(theta)*sin(phi) + sin(zeta)*cos(phi))
    k_hat(3) = cos(eta)*cos(theta) - sin(eta)*cos(zeta)*sin(theta)

    !initial radial direction
    r0_hat(1) = sin(theta)*cos(phi) ; r0_hat(2) = sin(theta)*sin(phi) ; r0_hat(3) = cos(theta)

    l_phot(1) = 0
    coord_phot(1, 1) = R * r0_hat(1) ; coord_phot(1, 2) = R * r0_hat(2) ; coord_phot(1, 3) = R * r0_hat(3)
    k_phot(1, :) = k_hat(:)

    b = b_imp(R, eta)
    b_var = b
    !normal vector to trajectory plane
    if (b == 0) then
        !set n_hat to theta_hat if k and r_hat are parallel
        n_hat(1) = cos(theta)*cos(phi); n_hat(2) = cos(theta)*sin(phi); n_hat(3) = -sin(theta)
    else
        call cross(k_hat, r0_hat, n_hat)
        n_hat = n_hat / norm2(n_hat)
    end if
    call cross(r0_hat, n_hat, psi0_hat)
    call b_ratio(b, c1, c2, iBs(1), iBs(2))
    bs(:)= b_data(iBs(:))
    !print *, c1, c2
    abscissa_increasing = .true.

! initialize the two splines to interpolate with
    r_spl(:) = c1*r_data(iBs(1), :) + c2*r_data(iBs(2),:)
    l_spl(:) = c1*l_data(iBs(1),:) + c2*l_data(iBs(2),:)
    psi_spl(:) = c1*psi_data(iBs(1), :) + c2*psi_data(iBs(2),:)
    dr_spl(:) = c1*dr_data(iBs(1),:) + c2*dr_data(iBs(2),:)
    rdpsi_spl(:) = c1*rdpsi_data(iBs(1),:) + c2*rdpsi_data(iBs(2),:)
    
    call spline(l_spl, r_spl, N_spoints, dr_spl(1), dr_spl(N_spoints), r_spl2)
    call spline(l_spl, psi_spl, N_spoints, rdpsi_spl(1)/r_spl(1), &
                rdpsi_spl(N_spoints)/r_spl(N_spoints), psi_spl2)
    call spline(r_spl, l_spl, N_spoints, 1/dr_spl(1), 1/dr_spl(N_spoints), l_spl2)
    fdata(1) = f(r_spl(1)); fdata(2) = f(r_spl(N_spoints))
    df(1) = dfdr(r_spl(1)); df(2) = dfdr(r_spl(N_spoints))
    drdl2 = -fdata * (1.d0 + fdata**2.d0)**(-2.d0) * df
    call spline(l_spl, dr_spl, N_spoints, drdl2(1), drdl2(2), dr_spl2)
    dpsidl2 = df / (1.d0 + fdata**2.d0)**(2.d0)
    call spline(l_spl, rdpsi_spl, N_spoints, dpsidl2(1), dpsidl2(2), rdpsi_spl2)

    if (R < r_spl(1)) then
        rm = R_min(b)
        l0 = (R**2.d0 - rm**2.d0)**0.5d0
        psi0 = l0 / rm
    else
        call splint(r_spl, l_spl, l_spl2, N_spoints, R, l0, abscissa_increasing)
        call splint(l_spl, psi_spl, psi_spl2, N_spoints, l0, psi0, abscissa_increasing)
    end if
    if (cos(eta) < 0) then
        l0 = -l0
        psi0 = -psi0
    end if
    if (b > b_inner) then
        Npath = 2*N_spoints
        do iL = 1, N_spoints
            i1 = N_spoints + 1 - iL !negative l spline
            i2 = N_spoints + iL !positive l spline
            l_phot(i1) = -l_spl(iL)
            l_phot(i2) = l_spl(iL)
            rad = r_spl(iL)
            psi1 = -psi_spl(iL) - psi0
            psi2 = psi_spl(iL) - psi0
            r1_hat = cos(psi1)*r0_hat + sin(psi1)*psi0_hat
            call cross(r1_hat, n_hat, psi1_hat)
            coord_phot(i1, :) = rad * r1_hat(:)
            k_phot(i1, :) = -dr_spl(iL)*r1_hat(:) + rdpsi_spl(iL)*psi1_hat(:)
            r2_hat = cos(psi2)*r0_hat + sin(psi2)*psi0_hat
            coord_phot(i2, :) = rad * r2_hat
            call cross(r2_hat, n_hat, psi2_hat)
            k_phot(i2, :) = dr_spl(iL)*r2_hat(:) + rdpsi_spl(iL)*psi2_hat(:)
        end do
        !determine dk/dl at start
        rdpsi = f(rad)
        df(1) = dfdr(rad)
        dk0 = ((df(1)/(1 + rdpsi**2.d0) + rdpsi/rad)/(1+rdpsi**2.d0)) * (psi1_hat - rdpsi * r1_hat)
        !determine dk/dl at end
        dkn = ((df(1)/(1 + rdpsi**2.d0) + rdpsi/rad)/(1+rdpsi**2.d0)) * (psi2_hat - rdpsi * r2_hat)
        !shift the middle of the array to be at R_min
        l_phot(N_spoints) = 0
        rad = R_min(b)
        psi1 = -psi0
        r1_hat = cos(psi1)*r0_hat + sin(psi1)*psi0_hat
        call cross(r1_hat, n_hat, psi1_hat)
        coord_phot(N_spoints, :) = rad*r1_hat(:)
        k_phot(N_spoints, :) = psi1_hat(:)
    else
        Npath = N_spoints
        if (cos(eta) < 0) then
            do iL = 1, N_spoints
                i1 = N_spoints + 1 - iL
                l_phot(i1) = -l_spl(iL)
                rad = r_spl(iL)
                psi1 = -psi_spl(iL) - psi0
                r1_hat = cos(psi1)*r0_hat + sin(psi1)*psi0_hat
                call cross(r1_hat, n_hat, psi1_hat)
                coord_phot(i1, :) = rad * r1_hat(:)
                k_phot(i1, :) = -dr_spl(iL)*r1_hat(:) + rdpsi_spl(iL)*psi1_hat(:)
                if (iL .eq. 1) then
                    rdpsi = f(rad)
                    df(1) = dfdr(rad)
                    dk0 = ((df(1)/(1 + rdpsi**2.d0) + rdpsi/rad)/(1+rdpsi**2.d0)) * (psi1_hat - rdpsi * r1_hat)
                else if (iL .eq. N_spoints) then
                    rdpsi = f(rad)
                    df(1) = dfdr(rad)
                    dkn = ((df(1)/(1 + rdpsi**2.d0) + rdpsi/rad)/(1+rdpsi**2.d0)) * (psi1_hat - rdpsi * r1_hat)
                end if

            end do
        else
            l_phot(1:N_spoints) = l_spl(1:N_spoints)
            do iL = 1, N_spoints
                rad = r_spl(iL)
                psi2 = psi_spl(iL) - psi0
                r2_hat = cos(psi2)*r0_hat + sin(psi2)*psi0_hat
                coord_phot(iL, :) = rad * r2_hat
                call cross(r2_hat, n_hat, psi2_hat)
                k_phot(iL, :) = dr_spl(iL)*r2_hat(:) + rdpsi_spl(iL)*psi2_hat(:)
                if (iL .eq. 1) then
                    rdpsi = f(rad)
                    df(1) = dfdr(rad)
                    dk0 = ((df(1)/(1 + rdpsi**2.d0) + rdpsi/rad)/(1+rdpsi**2.d0)) * (psi2_hat - rdpsi * r2_hat)
                else if (iL .eq. N_spoints) then
                    rdpsi = f(rad)
                    df(1) = dfdr(rad)
                    dkn = ((df(1)/(1 + rdpsi**2.d0) + rdpsi/rad)/(1+rdpsi**2.d0)) * (psi2_hat - rdpsi * r2_hat)
                end if
            end do
        end if
    end if
    !determine what l the spline ends at
    if ((cos(eta) .ge. 0) .or. (b > b_crit)) then
        call splint(r_spl, l_spl, l_spl2, N_spoints, Rmax, ln, abscissa_increasing)
    else
        call splint(r_spl, l_spl, l_spl2, N_spoints, 1.d0, ln, abscissa_increasing)
        ln = -ln
    end if
    return
    end subroutine traj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end module trajectory
