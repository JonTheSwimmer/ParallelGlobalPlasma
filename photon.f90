!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!:
    module photon
    use parameters_base
    use routines
    use spectrum
    use grid
    use OMP_LIB

    real*8 :: beta_var !storage for beta value across integral
    real*8 :: N_emission !photon emission rate, based on average energy of photons in initial spectrum 
    contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set N_emission for photon creation
! must be run before rate is determined! Otherwise everything blows up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine init_photon
    implicit none
    real*8 :: E_char
    E_char = omega_avg(omega_min, omega_max)
    print *, "E_char =", E_char
    N_emission = luminosity/E_char

    end subroutine init_photon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! mu_i = k_i * B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function mu_i(k, B)
    implicit none
    real*8, dimension(3), intent(in) :: k, B
    real*8 :: mu_i
    mu_i = sum(k*B)
    return
    end function mu_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! cross section for two photons in background magnetic field, including flux rate
! omega_i should be in units of ergs
! B is in units of B_Q
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function sigma_c(omega1, omega2, k1, k2, B)
    implicit none
    real*8, intent(in) :: omega1, omega2
    real*8, dimension(3), intent(in) :: k1, k2, B
    real*8 :: C1, mu1, mu2, mu12
    real*8, dimension(3) :: bhat
    real*8 :: sigma_c
    bhat = B / norm2(B)

    mu1 = mu_i(k1, Bhat)
    mu2 = mu_i(k2, Bhat)
    C1 = C_plus(omega1, omega2, mu1, mu2)
    if (norm2(B) .ge. 1.d0) then
        if (C1 > 4*m_e**2) then
            sigma_c = 32.d0 * pi * (r_elec**2.d0) * norm2(B) &
            * sqrt(C1 - 4*m_e**2) / (C1**0.5d0 * (omega1*omega2)**3.d0) &
            *(C1**2.d0 * (1-mu1**2) * (1-mu2**2)*m_e**6) &
            /((C1*(1-mu1**2)*(1-mu2**2) + 4*m_e**2*(mu1 - mu2)**2)**2.d0)
        else
            sigma_c = 0
        end if
    else 
        if (E_com(omega1, omega2, k1, k2) > m_e) then
            sigma_c = sigma0(omega1, omega2, k1, k2)
        else 
            sigma_c = 0.d0
        end if
    end if
    return
    end function sigma_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function C_plus(omega1, omega2, mu1, mu2)
    implicit none
    real*8, intent(in) :: omega1, omega2, mu1, mu2
    real*8 :: C_plus

    C_plus = (omega1 + omega2)**2.d0 - (mu1*omega1 + mu2*omega2)**2.d0
    return
    end function C_plus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function C_perp(omega, mu)
    implicit none
    real*8, intent(in) :: omega, mu
    real*8 :: C_perp

    C_perp = (1 - mu**2.d0) * omega**2.d0
    return
    end function C_perp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function C_cross(omega1, omega2, mu1, mu2)
    implicit none
    real*8, intent(in) :: omega1, omega2, mu1, mu2
    real*8 :: C_cross
    C_cross = omega1 * omega2 * (1 - mu1 * mu2)

    return
    end function C_cross

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! cross section for two photons without magnetic field
! Assumes photons are in parallel polarization (see Breit and Wheeler 1934)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function sigma0(omega1, omega2, k1, k2)
    implicit none
    real*8, intent(in) :: omega1, omega2
    real*8, dimension(3), intent(in) :: k1, k2
    real*8 :: sigma0, theta
    theta = acosh(E_com(omega1, omega2, k1, k2) / m_e)
    sigma0 = pi * (r_elec**2.d0)
    sigma0 = sigma0 * (2.d0*theta*(cosh(theta)**(-2) + cosh(theta)**(-4) - 0.75d0 * cosh(theta)**(-6)) &
                - sinh(theta) * (cosh(theta)**(-3) + 1.5d0 * cosh(theta)**(-5)))
    sigma0 = sigma0 * (1.d0 - sum(k1 * k2))
    return
    end function sigma0

    function E_com(omega1, omega2, k1, k2)
    implicit none
    real*8, intent(in) :: omega1, omega2
    real*8, dimension(3), intent(in) :: k1, k2
    real*8 :: E_com
    E_com = sqrt(omega1 * omega2 * (1.d0 - sum(k1 * k2)) / 2.d0)
    return
    end function E_com


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solid angle given bounding coordinates
! theta2 > theta1, phi2 > phi1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function s_angle(theta1, theta2, phi1, phi2)
    implicit none
    real*8, intent(in) :: theta1, theta2, phi1, phi2
    real*8 :: s_angle
    s_angle = (cos(theta1) - cos(theta2)) * (phi2 - phi1)
    return
    end function s_angle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! creation rate dR/dt for a given spacial cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rate(indR, indTheta, indPhi, N_rate, E_rate)
    implicit none
    !input variables
    integer, intent(in) :: indR, indTheta, indPhi
    !output variables
    real*8, intent(out) :: N_rate, E_rate
    !local variables
    integer :: iT1, iT2, iP1, iP2, iW1, iW2
    real*8 :: n1, n2, l_tot, n_tot, p_theta1, p_theta2, p_phi1, p_phi2, omega1, omega2
    real*8 :: Ndot
    real*8, dimension(3) :: k1, k2, B
    N_rate = 0.d0; E_rate = 0.d0
    l_tot = cell_time(indR, indTheta, indPhi)
    n_tot = (l_tot / (c_light * cell_volume(indR, indTheta, indPhi))) * (N_emission / N_prop)
    if (n_tot .eq. 0) then
        return
    end if 
    call cellB(indR, indTheta, indPhi, B)
    do iT1 = 1, Nk_theta
        do iP1 = 1, Nk_phi(iT1)
            do iT2 = 1, Nk_theta
                do iP2 = 1, Nk_phi(iT2)
                    p_theta1 = acos((cos(ktheta1(iT1)) + cos(ktheta2(iT1))) / 2)
                    p_theta2 = acos((cos(ktheta1(iT2)) + cos(ktheta2(iT2))) / 2)
                    p_phi1 = (kphi1(iT1, iP1) + kphi2(iT1, iP1)) / 2
                    p_phi2 = (kphi1(iT2, iP2) + kphi2(iT2, iP2)) / 2
                    call cartesian(1.d0, p_theta1, p_phi1, k1(1), k1(2), k1(3))
                    call cartesian(1.d0, p_theta2, p_phi2, k2(1), k2(2), k2(3))
                    do iW1 = 1, Nomega
                        omega1 = omega_mid(iW1)
                        n1 = n_tot * cell_bins(indR, indTheta, indPhi, iT1, iP1, iW1) / l_tot
                        do iW2 = 1, Nomega
                            omega2 = omega_mid(iW2)
                            n2 = n_tot * cell_bins(indR, indTheta, indPhi, iT2, iP2, iW2) / l_tot
                            if ((n1 .gt. 0) .and. (n2 .gt. 0)) then
                                Ndot = (n1 * n2 * sigma_c(omega1, omega2, k1, k2, B) * c_light)
                                N_rate = N_rate + Ndot 
                                E_rate = E_rate + Ndot*(omega1 + omega2)
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do
    N_rate = N_rate * cell_volume(indR, indTheta, indPhi)
    E_rate = E_rate * cell_volume(indR, indTheta, indPhi)
    return
    end subroutine rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! creation rate for all spacial cells 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rate_fill(Edot)
    implicit none
    !output variable
    real*8, intent(out) :: Edot !erg/s of photon annihilation to pairs
    real*8 :: Eloss
    integer :: indR, indT, indP
    Edot = 0.d0
    do indR = 1, Nr+1
        do indT = 1, 2*Ntheta
            do indP = 1, Nphi
                if (cell_volume(indR, indT, indP) > 0) then
                    call rate(indR, indT, indP, cell_create(indR, indT, indP), Eloss) 
                    Edot = Edot + Eloss
                end if
            end do
        end do
    end do
    end subroutine rate_fill

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! creation rate for all spacial cells, done in parallel 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rate_para(Edot)
    implicit none
    !output variable
    real*8, intent(out) :: Edot !erg/s of photon annihilation to pairs
    real*8 :: Ehold, Eloss
    integer :: indR, indT, indP
    Edot = 0.d0
    do indR = 1, Nr+1
        do indT = 1, 2*Ntheta
            !$OMP PARALLEL PRIVATE(indP, Eloss, Ehold)
            Eloss = 0.d0
            !$OMP DO
            do indP = 1, Nphi
                if (cell_volume(indR, indT, indP) > 0) then
                    call rate(indR, indT, indP, cell_create(indR, indT, indP), Ehold)
                    Eloss = Eloss + Ehold
                end if
            end do
            !$OMP END DO
            !$OMP CRITICAL
            Edot = Edot + Eloss
            !$OMP END CRITICAL
            !$OMP END PARALLEL
        end do
    end do
    end subroutine rate_para

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! export creation rate data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine output_rate(runnum)
    implicit none
    integer, intent(in) :: runnum
    integer :: iR, iT
    character(len=8) :: num ! format descriptor
    write (num, "(I0)") runnum
    call file_open(2, trim(num)//'_cell_create.dat')
101 format(*(g0, :, ", "))
    do iR = 1, Nr+1
        do iT = 1, 2*Ntheta
            write (2, 101) cell_create(iR, iT, :)
        end do
    end do
    close (2)
    return
    end subroutine output_rate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! surface annihilation rate for a given flux bundle
! assumes n_pair is in n/R_NS**3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function ann_surface(indR, indP, n_pair)
    implicit none
    integer :: indR, indP
    real*8 :: n_pair, x1, x2, dphi
    real*8 :: ann_surface
    if (indR .lt. Nr+1) then !non-polar cells
        dphi = twopi / Nphi
        x1 = cos(theta02(indR))
        x2 = cos(theta01(indR))
    else    !polar cap cells
        dphi = pi !2pi/2 because polar cap only has one annihilation surface
        x1 = 1.d0
        x2 = cos(theta02(Nr))
    end if
    ann_surface = dphi * (2.d0/3.d0) * n_pair * beta_surf * c_light * (sqrt(3*x1**2.d0 + 1) - sqrt(3*x2**2.d0 + 1)) 
    return
    end function ann_surface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! outflow loss rate for a given flux bundle
! x-ray energy assumed to be monochromatic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function ann_outflow(indR, indP, n_pair, beta)
    implicit none
    !input variables
    integer :: indR, indP
    real*8 :: n_pair, beta
    !local variables
    real*8 :: x1, x2, dphi
    real*8 :: ann_outflow
    real*8 :: R_out !radius at which resonant scattering occurs
    R_out = (m_e/omega_out)**(1.d0 / 3) * (B_dip * 2.d0)**(1.d0/3)
    
    if (indR .lt. Nr+1) then !non-polar cells
        if (rmax_mid(indR) > R_out) then
            dphi = twopi / Nphi
            x1 = cos(theta02(indR))
            x2 = cos(theta01(indR))
        else !bundle does not experience resonant scattering
            ann_outflow = 0
            return
        end if
    else    !polar cap cells will always have loss to outflow
        dphi = pi !2pi/2 because polar cap only has one annihilation surface
        x1 = 1.d0
        x2 = cos(theta02(Nr))
    end if
    ann_outflow = dphi * (2.d0/3.d0) * n_pair * beta * c_light * (sqrt(3*x1**2.d0 + 1) - sqrt(3*x2**2.d0 + 1))
    return
    end function ann_outflow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! average cross section for a spatial cell, weighted by angular momentum bins 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function sigma_avg(indR, indTheta, indPhi)
    implicit none
    integer, intent(in) :: indR, indTheta, indPhi
    integer :: iT1, iT2, iP1, iP2, iW1, iW2
    real*8 :: A1, A2, n1, n2, l_tot, n_tot, p_theta1, p_theta2, p_phi1, p_phi2, omega1, omega2
    real*8 :: sigma_avg
    real*8, dimension(3) :: k1, k2, B
    sigma_avg = 0.d0
    l_tot = cell_time(indR, indTheta, indPhi)
    if (l_tot .eq. 0) then
        return
    end if 
    call cellB(indR, indTheta, indPhi, B)
    do iT1 = 1, Nk_theta
        do iP1 = 1, Nk_phi(iT1)
            do iT2 = 1, Nk_theta
                do iP2 = 1, Nk_phi(iT2)
                    p_theta1 = acos((cos(ktheta1(iT1)) + cos(ktheta2(iT1))) / 2)
                    p_theta2 = acos((cos(ktheta1(iT2)) + cos(ktheta2(iT2))) / 2)
                    p_phi1 = (kphi1(iT1, iP1) + kphi2(iT1, iP1)) / 2
                    p_phi2 = (kphi1(iT2, iP2) + kphi2(iT2, iP2)) / 2
                    call cartesian(1.d0, p_theta1, p_phi1, k1(1), k1(2), k1(3))
                    call cartesian(1.d0, p_theta2, p_phi2, k2(1), k2(2), k2(3))
                    do iW1 = 1, Nomega
                        omega1 = omega_mid(iW1)
                        n1 = cell_bins(indR, indTheta, indPhi, iT1, iP1, iW1) / l_tot
                        do iW2 = 1, Nomega
                            omega2 = omega_mid(iW2)
                            n2 = cell_bins(indR, indTheta, indPhi, iT2, iP2, iW2) / l_tot
                            if ((n1 .gt. 0) .and. (n2 .gt. 0)) then
                                sigma_avg = sigma_avg + (n1 * n2 * sigma_c(omega1, omega2, k1, k2, B))
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do
    return
    end function sigma_avg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! annihilation cross section in each cell 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function sigma_ann(indR, indTheta, indP, beta)
    implicit none
    integer, intent(in) :: indR, indTheta, indP
    real*8 :: sigma_ann, beta, gam, sigmaB0
    real*8, dimension(3) :: B
    call cellB(indR, indTheta, indP, B)
    sigma_ann = min(sigma_annB(B, beta), sigma_ann0(beta))
    return
    end function sigma_ann

    function sigma_annB(B, beta)
    implicit none
    real*8 :: sigma_annB, beta, gam
    real*8, dimension(3) :: B
    sigma_annB = 8*pi*r_elec**2.d0 / (norm2(B))
    gam = (1 - beta**2.d0)**(-0.5d0)
    if (gam .le. 2.d0) then
        sigma_annB =sigma_annB*(beta/gam**2)*((5.d0/9)-(1.29116*beta**2)-(0.0886*beta**4) - &
                log(beta)*((4.d0/3)+(0.93634*beta**2.d0)+(4.53916*beta**4)))
    else
        sigma_annB = sigma_annB*(gam**-6.d0)*(1.d0 + 2.84/gam**2 + 82/gam**6)
    end if
    return
    end function sigma_annB

    function sigma_ann0(beta)
    implicit none
    real*8 :: sigma_ann0, beta, f
    f = (1-beta**2)*(4*beta)**(-1.d0)*(beta**(-1.d0)* (3-beta**4)*log((1+beta)/(1-beta)) - 2*(2-beta**2))
    sigma_ann0 = pi * r_elec**2.d0 * f
    return
    end function sigma_ann0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes the prefactor for N_ann = C*n(R_NS)**2
! n(R_NS) should be in units of 1/R_NS**3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function ann_average(indR, indP, beta)
    implicit none
    integer, intent(in)  :: indR, indP
    real*8, intent(in) :: beta
    real*8 :: sigma_const, B(3), dphi, N1(1), N2(1)
    real*8 :: ann_average, theta_max
    beta_var = beta
    call cellB(1, Ntheta, 1, B)
    !sigma_const = sigma_ann(1, Ntheta, 1, beta) * norm2(B) / B_dip !prefactor for sigma_ann without magnetic field
    sigma_const = sigma_annB(B, beta) * norm2(B)/B_dip !prefactor for sigma_ann without magnetic field
    N1(1) = 0; N2(1) = 0
    if (indR .eq. Nr+1) then !polar cap 
        dphi = twopi
        theta_max = asin(sqrt(Rmax*sin(theta02(Nr))**2.d0))
        call integrate_ode(1, N1, 0.d0, theta_max, deriv_capB)
        !call integrate_dlsode(1, N2, 0.d0, theta_max, deriv_cap0)
        call integrate_ode(1, N2, 1.d0, Rmax, deriv_cap0_r)
    else
        dphi = twopi/Nphi
        call integrate_ode(1, N1, rmax1(indR), rmax2(indR), deriv_annB)
        call integrate_ode(1, N2, rmax1(indR), rmax2(indR), deriv_ann0)
    end if 
    N1(1) = N1(1) * sigma_const
    N2(1) = N2(1) * sigma_ann0(beta)
    ann_average = (c_light * dphi * beta) * (N1(1) + N2(1))
    return
    end function ann_average

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes the prefactor for N_ann = C*n(R_NS)**2
! n(R_NS) should be in units of 1/R_NS**3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function ann_trapz(indR, indP, beta)
    implicit none
    integer, intent(in)  :: indR, indP
    real*8, intent(in) :: beta
    real*8 :: sigma_const, B(3), dphi
    real*8 :: ann_trapz, theta_max, R1, R2, N1, N2
    integer :: indT
    integer, parameter :: pts = 1000
    real*8, dimension(pts) :: dN1dx, dN2dx, X
    real*8, dimension(1):: Y, dYdx
    beta_var = beta
    call cellB(1, Ntheta, 1, B)
    sigma_const = sigma_ann(1, Ntheta, 1, beta) * norm2(B) / B_dip !prefactor for sigma_ann without magnetic field
    Y(1) = 1.d0 
    if (indR .eq. Nr+1) then !polar cap 
        dphi = twopi
        theta_max = asin(sqrt(Rmax*sin(theta02(Nr))**2.d0))
        do indT = 1, pts
            X(indT) = theta_max * (indT-1)/(pts-1)
            call deriv_capB(1, X(indT), Y, dYdx)
            dN1dx(indT) = dYdx(1)
            call deriv_cap0(1, X(indT), Y, dYdx)
            dN2dx(indT) = dYdx(1)
        end do
    else
        dphi = twopi/Nphi
        R1 = rmax1(indR)
        R2 = rmax2(indR)
        do indT = 1, pts
            X(indT) = R1 + (R2 - R1) * (indT - 1) / (pts-1)
            call deriv_annB(1, X(indT), Y, dYdx)
            dN1dx(indT) = dYdx(1)
            call deriv_ann0(1, X(indT), Y, dYdx)
            dN2dx(indT) = dYdx(1)
        end do
    end if
    N1 = 0; N2 = 0;
    do indT = 1, pts-1
        N1 = N1 + (dN1dx(indT+1)+dN1dx(indT))*(X(indT+1)-X(indT))/2.d0
        N2 = N2 + (dN2dx(indT+1)+dN2dx(indT))*(X(indT+1)-X(indT))/2.d0
    end do
    N1 = N1 * sigma_const
    N2 = N2 * sigma_ann0(beta)
    ann_trapz = (c_light * dphi * beta) * (N1 + N2)
    return
    end function ann_trapz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! integrand for r_max annihilation integral when magnetic field is dominant
! X should be dimension 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine deriv_annB(neq, R_in, X, dXdr)
    integer, intent(in) :: neq
    real*8, intent(in) :: R_in, X(neq)
    real*8, intent(out) :: dXdr(neq)
    real*8 :: t0, t1, x0, x1
    t0 = asin(sqrt(1.d0/R_in))
    if (R_in .gt. Rmax) then
        t1 = asin(sqrt(Rmax/R_in))
    else
        t1 = pi/2.d0
    end if
    t1 = min(t1, thetaB0(R_in, beta_var)) 
    x0 = cos(t0) ; x1 = cos(t1)
    if (t0 .lt. t1) then
        dXdr(1) = (4 * R_in - 3)**(-1.d0) * 2*(indef_theta(x0) - indef_theta(x1))
    else
        dXdr(1) = 0.d0
    end if
    return
    end subroutine deriv_annB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! integrand for r_max annihilation integral when magnetic field is negligible
! X should be dimension 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine deriv_ann0(neq, R_in, X, dXdr)
    integer, intent(in) :: neq
    real*8, intent(in) :: R_in, X(neq)
    real*8, intent(out) :: dXdr(neq)
    real*8 :: t0, t1, x0, x1
    t0 = asin(sqrt(1.d0/R_in))
    t0 = max(t0, thetaB0(R_in, beta_var))
    if (R_in .gt. Rmax) then
        t1 = asin(sqrt(Rmax/R_in))
    else
        t1 = pi/2.d0
    end if
    x0 = cos(t0) ; x1 = cos(t1)
    if (t0 .lt. t1) then
        dXdr(1) = (R_in**3.d0*(4*R_in - 3))**(-1.d0) * 2*(indef2_theta(x0) - indef2_theta(x1))
    else
        dXdr(1) = 0.d0
    end if
    return
    end subroutine deriv_ann0 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! indefinite integral over theta for the magnetic field volume integral
! x = cos(theta)
! integral of sqrt(1+3x**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function indef_theta(x)
    real*8 :: x, indef_theta
    real*8 :: y
    y = sqrt(3.d0) * x
    indef_theta = (3*x*sqrt(y**2.d0+1) + sqrt(3.d0)*asinh(y))/6.d0
    return
    end function indef_theta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! indefinite integral over theta for the magnetic field squared volume integral
! x = cos(theta)
! integral of sqrt(1+3x**2) / (1 - x**2)**3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function indef2_theta(x)
    real*8 :: x, indef2_theta
    indef2_theta = (x/(x**2.d0 - 1)**2.d0)
    return
    end function indef2_theta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! integrand for cap annihilation integral when magnetic field is dominant
! X should be dimension 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine deriv_capB(neq, T_in, X, dXdT)
    integer, intent(in) :: neq
    real*8, intent(in) :: T_in, X(neq)
    real*8, intent(out) :: dXdT(neq)
    real*8 :: r0, r1, B(3), u
    r0 = max(1.d0, (sin(T_in)/sin(theta02(Nr)))**2.d0)
    call Bfield(1.d0, T_in, 0.d0, B(1), B(2), B(3))
    r1 = (sigma_ann0(beta_var) / sigma_annB(B, beta_var))**(1.d0/3.d0)
    r1 = min(Rmax, r1)
    !print *, T_in, r0, r1
    if (r0 .lt. r1) then
        u = sin(T_in)
        dXdT(1) = u*sqrt(4.d0-3.d0*u**2.d0)*log((4.d0*r1-3*u**2.d0)/(4.d0*r0-3*u**2))/4.d0
        !print *, T_in, r0, r1, log((4.d0*r1-3*u**2.d0)/(4.d0*r0-3*u**2))/4.d0
    else
        dXdT(1) = 0.d0
    end if
    return
    end subroutine deriv_capB 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! integrand for cap annihilation integral when magnetic field is negligible
! X should be dimension 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine deriv_cap0(neq, T_in, X, dXdT)
    integer, intent(in) :: neq
    real*8, intent(in) :: T_in, X(neq)
    real*8, intent(out) :: dXdT(neq)
    real*8 :: r0, r1, B(3), u
    r0 = max(1.d0, (sin(T_in)/sin(theta02(Nr)))**2.d0)
    call Bfield(1.d0, T_in, 0.d0, B(1), B(2), B(3))
    r0 = max(r0, (sigma_ann0(beta_var) / sigma_annB(B, beta_var))**(1.d0/3.d0))
    r1 = Rmax
    !print *, T_in, r0, r1
    if (r0 .lt. r1) then
        u = sin(T_in)
        dXdT(1) = u*(4.d0-3.d0*u**2.d0)*(indef2_r(r1, u) - indef2_r(r0, u)) 
        !print *, T_in, r0, r1, dXdT(1)
    else
        dXdT(1) = 0.d0
    end if
    return
    end subroutine deriv_cap0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! integrand for cap annihilation integral when magnetic field is negligible
! X should be dimension 1
! integrating over r, then theta because the other way is numerically unstable >:(
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine deriv_cap0_r(neq, R_in, X, dXdT)
    integer, intent(in) :: neq
    real*8, intent(in) :: R_in, X(neq)
    real*8, intent(out) :: dXdT(neq)
    real*8 :: t0, t1, sigB, sig0
    real*8 :: u0, u1, A
    real*8, dimension(3) :: B
    call Bfield(R_in, 0.d0, 0.d0, B(1), B(2), B(3))
    sig0 = sigma_ann0(beta_var)
    sigB = sigma_annB(B, beta_var)
    if (sigB .ge. sig0) then
        !entire theta section is non-magnetic
        t0 = 0.d0
    else if (sigB .le. sig0/2.d0) then 
        !entire section is magnetic
        t0 = pi/2
    else
        t0 = acos(sqrt(((2.d0*sigB/sig0)**2.d0 - 1.d0)/3.d0))
    end if
    t1 = asin(sqrt(R_in * sin(theta02(Nr))**2.d0))
    !print *, t0, t1
    !print *, T_in, r0, r1
    if (t0 .lt. t1) then
        u0 = cos(t0)
        u1 = cos(t1)
        A = 4 * R_in - 3
        dXdT(1) = R_in**(-3.d0) * (u0 - u1 - (A-1.d0)/sqrt(3.d0*A) * &
                (atan(sqrt(3/A)*u0)- atan(sqrt(3/A)*u1)))
        !print *, T_in, r0, r1, dXdT(1)
    else
        dXdT(1) = 0.d0
    end if
    return
    end subroutine deriv_cap0_r


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! indefinite integral over r for the non-magnetic cap annihilation integral
! ST = sin(theta)
! integral of (1/x^3) * (4x - 3 sin^2(theta))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function indef2_r(x, ST)
    !input variables
    real*8, intent(in) :: x, ST
    !output variable
    real*8 :: indef2_r
    !local variables
    real*8 :: y
    y = 3*ST**2.d0
    if (y .eq. 0.d0) then
        indef2_r = (-12.d0 * x**3.d0)**(-1.d0)
    else 
        indef2_r = (3.2d1*x**2*(log(4.d0-y/x))+y*(8*x+y))/(2*y**3*x**2)
    end if 
    return
    end function indef2_r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Along a given flux bundle, return the theta such that the B=0 annihilation cross section is equal to
! the normal cross section with strong magnetic fields
! input parameters R = R_max of flux bundle, beta = electron beta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function thetaB0(Rmax, beta)
    implicit none
    real*8 :: thetaB0, Rmax, beta
    real*8 :: t_low, t_high, t_mid
    real*8 :: sig_diff, r, tol, sig0
    real*8, dimension(3) :: B
    integer :: loop_count
    logical :: cont 
    t_low = asin(sqrt(1.d0/Rmax))
    t_high = pi/2
    tol = 0.0001
    cont = .true.
    loop_count = 0
    sig0 = sigma_ann0(beta)
    ! check if B=0 annihilation cross section is larger or smaller along entire flux tube
    r = 1.d0
    call Bfield(r, t_low, 0.d0, B(1), B(2), B(3)) 
    sig_diff = sigma_annB(B, beta) - sig0
    if (sig_diff .ge. 0) then
        thetaB0 = t_low
        return
    else
        call Bfield(Rmax, t_high, 0.d0, B(1), B(2), B(3))
        sig_diff = sigma_annB(B, beta) - sig0
        if (sig_diff .le. 0) then
            thetaB0 = t_high
            return
        end if
    end if
    do while (cont)
        t_mid = (t_high + t_low) / 2.d0
        r = Rmax * sin(t_mid)**2.d0
        call Bfield(r, t_mid, 0.d0, B(1), B(2), B(3))
        sig_diff = sigma_annB(B, beta) - sig0
        if (sig_diff .lt. 0) then
            t_low = t_mid
        else if (sig_diff .gt. 0) then
            t_high = t_mid
        else !if sig_diff is exactly 0
            thetaB0 = t_mid
            return
        end if

        if ((loop_count .gt. 1000) .or. ((t_high - t_low)/t_high .lt. tol)) then 
            cont = .false.
        end if
        loop_count = loop_count + 1
    end do
    thetaB0 = t_mid
    return
    end function thetaB0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! redistribute photons in each energy bin according to initial spectrum
! does not conserve energy!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine dist_spec()
    integer :: iR, iT, iP, iKT, iKP
    real*8 :: phot_sum
    do iR = 1, Nr+1
        do iT = 1, 2*Ntheta
            do iP = 1, Nphi
                do iKT = 1, Nk_theta
                    do iKP = 1, Nk_phi(iKT)
                        phot_sum = sum(cell_bins(iR, iT, iP, iKT, iKP, :))
                        cell_bins(iR, iT, iP, iKT, iKP, :) = phot_sum * spec_weights(:)
                    end do
                end do
            end do
        end do
    end do
    return
    end subroutine dist_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! percent of collisions above threshold for pair conversion
! for magnetic cells, threshold is lower
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function threshold(indR, indTheta, indPhi)
    implicit none
    integer, intent(in) :: indR, indTheta, indPhi
    integer :: iT1, iT2, iP1, iP2, iW1, iW2
    real*8 :: n1, n2, l_tot, n_tot, p_theta1, p_theta2, p_phi1, p_phi2, omega1, omega2
    real*8 :: threshold, mu1, mu2, mu12, B_norm
    real*8, dimension(3) :: k1, k2, B
    threshold = 0.d0
    l_tot = cell_time(indR, indTheta, indPhi)
    if (l_tot .eq. 0) then
        return
    end if 
    call cellB(indR, indTheta, indPhi, B)
    B_norm = norm2(B)
    do iT1 = 1, Nk_theta
        do iP1 = 1, Nk_phi(iT1)
            do iT2 = 1, Nk_theta
                do iP2 = 1, Nk_phi(iT2)
                    p_theta1 = acos((cos(ktheta1(iT1)) + cos(ktheta2(iT1))) / 2)
                    p_theta2 = acos((cos(ktheta1(iT2)) + cos(ktheta2(iT2))) / 2)
                    p_phi1 = (kphi1(iT1, iP1) + kphi2(iT1, iP1)) / 2
                    p_phi2 = (kphi1(iT2, iP2) + kphi2(iT2, iP2)) / 2
                    call cartesian(1.d0, p_theta1, p_phi1, k1(1), k1(2), k1(3))
                    call cartesian(1.d0, p_theta2, p_phi2, k2(1), k2(2), k2(3))
                    mu1 = mu_i(k1, B)/B_norm
                    mu2 = mu_i(k2, B)/B_norm
                    mu12 = mu_i(k1, k2)
                    do iW1 = 1, Nomega
                        omega1 = omega_mid(iW1)
                        n1 = cell_bins(indR, indTheta, indPhi, iT1, iP1, iW1) / l_tot
                        do iW2 = 1, Nomega
                            omega2 = omega_mid(iW2)
                            n2 = cell_bins(indR, indTheta, indPhi, iT2, iP2, iW2) / l_tot
                            if (B_norm .ge. 1) then
                                if (((omega1+omega2)**2 - (mu1*omega1 + mu2*omega2)**2) .ge. 4.d0*m_e**2) then
                                    threshold = threshold + n1*n2
                                end if
                            else
                                !if ((2*omega1*omega2*(1.d0-mu12)) .ge. 4.d0*m_e**2) then
                                if (E_com(omega1, omega2, k1, k2) .ge. m_e) then
                                    threshold = threshold + n1*n2
                                end if
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do

    return
    end function threshold
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! determine surface flux for each flux bundle for photon energy
! final units in erg/(s cm^2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function photon_flux(indR, indPhi)
    real*8 :: photon_flux
    integer, intent(in) :: indR, indPhi
    real*8 :: T01, T02, P1, P2, vol, KT, KP, mu
    real*8, dimension(Nk_theta, Nk_phimax, Nomega) :: phot_data
    integer :: indTS, indKT, indKP
    real*8, dimension(3) :: n_hat, k_hat
    logical :: first
    phot_data = 0; vol = 0; photon_flux = 0
    !determine what theta cells are in the bundle touch the surface of the star
    if (indR .lt. Nr+1) then
        T01 = theta01(indR); T02 = theta02(indR); P1 = phi1(indPhi); P2 = phi2(indPhi)
        call cartesian(1.d0, (T01 + T02)/2.d0, (P1 + P2)/2.d0, n_hat(1), n_hat(2), n_hat(3)) 
    else
        T01 = theta02(Nr); T02 = 0.d0; P1 = 0.d0; P2 = 0.d0;
        n_hat = 0
        n_hat(3) = 1.d0 
    end if
    first = .true.
    do indTS = 1, Ntheta       
        if (cell_volume(indR, indTS, indPhi) .gt. 0) then
            if ((theta_star1(indTS) .le. T01) .or. first) then
                phot_data(:, :, :) = phot_data(:, :, :) + cell_bins(indR, indTS, indPhi, :, :, :)
                vol = vol + cell_volume(indR, indTS, indPhi)*R_NS**3.d0
                first = .false.
            end if
        end if
    end do
    if (sum(phot_data) .eq. 0) then
        photon_flux = 0.d0
        return
    end if
    !convert to photon density in cm^-3
    phot_data = (phot_data / (c_light * vol)) * (N_emission / N_prop)
    do indKT = 1, Nk_theta
        KT = (ktheta1(indKT) + ktheta2(indKT))/2.d0
        do indKP = 1, Nk_phi(indKT)
            KP = (kphi1(indKT, indKP) + kphi2(indKT, indKP))/2.d0
            call cartesian(1.d0, KT, KP, k_hat(1), k_hat(2), k_hat(3))
            mu = sum(n_hat*k_hat)
            if (mu .lt. 0) then
                photon_flux = photon_flux + abs(mu) * sum(phot_data(indKT, indKP, :) * omega_mid(:))
            end if
        end do
    end do
    photon_flux = photon_flux * 0.5 * c_light*R_NS
    return
    end function photon_flux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! determine surface temp from deposited energy
! final units in K
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function flux2temp(indR, indNS, indPhi)
    integer, intent(in):: indR, indNS, indPhi
    real*8 :: flux2temp
    real*8 :: area, time, dphi, dtheta
    if (indR .eq. Nr + 1) then
        dphi = twopi
        dtheta = 1.d0 - cos(theta02(Nr)) 
    else
        dphi = twopi/Nphi
        dtheta = cos(theta02(indR)) - cos(theta01(indR))
    end if
    area = R_NS**2.d0 * dphi * dtheta
    time = N_prop / N_emission
    flux2temp = (surf_eff * surf_flux(indR, indNS, indPhi) / (area * time* sigma_SB))**(0.25d0)

    end function flux2temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! zero out bins where the photon energy exceeds the single photon conversion criteria
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine single_photon_filter()
    integer :: iR, iT, iP, iKT, iKP, iW
    real*8, dimension(3)  :: B, Bhat, khat
    real*8 :: KT, KP, sin_kb
    do iR = 1, Nr + 1
        do iT = 1, 2*Ntheta
            do iP = 1, Nphi
                if (cell_volume(iR, iT, iP) .gt. 0) then
                call cellB(iR, iT, iP, B)
                Bhat = B/norm2(B)
                do iKT = 1, Nk_theta
                    KT = (ktheta1(iKT) + ktheta2(iKT))/2.d0
                    do iKP = 1, Nk_phi(iKT)
                        KP = (kphi1(iKT, iKP) + kphi2(iKT, iKP))/2.d0
                        call cartesian(1.d0, KT, KP, khat(1), khat(2), khat(3))
                        sin_kb = (1.d0 - sum(khat*Bhat)**2.d0)**(0.5d0)
                        do iW = 1, Nomega
                            if (omega_mid(iW) * sin_kb .ge. 2*m_e) then
                                cell_bins(iR, iT, iP, iKT, iKP, iW) = 0.d0
                            end if
                        end do
                    end do
                end do
                end if
            end do
        end do
    end do
    return
    end subroutine single_photon_filter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! return pair production rate for each photon energy bin combo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rate_omega_array(indR, indTheta, indPhi) result(rate_arr)
    implicit none
    integer, intent(in) :: indR, indTheta, indPhi
    integer :: iT1, iT2, iP1, iP2, iW1, iW2
    real*8, dimension(Nomega, Nomega) :: rate_arr
    real*8 :: n1, n2, l_tot, n_tot, p_theta1, p_theta2, p_phi1, p_phi2, omega1, omega2
    real*8, dimension(3) :: k1, k2, B
    real*8 :: vol, rate
    rate_arr = 0
    vol = cell_volume(indR, indTheta, indPhi)
    l_tot = cell_time(indR, indTheta, indPhi)
    n_tot = (l_tot / (c_light * cell_volume(indR, indTheta, indPhi))) * (N_emission / N_prop)
    if (n_tot .eq. 0) then
        return
    end if
    call cellB(indR, indTheta, indPhi, B)
    do iT1 = 1, Nk_theta
        do iP1 = 1, Nk_phi(iT1)
            do iT2 = 1, Nk_theta
                do iP2 = 1, Nk_phi(iT2)
                    !A1 = s_angle(ktheta1(iT1), ktheta2(iT1), kphi1(iT1, iP1), kphi2(iT1, iP1)) / (4 * pi)
                    !A2 = s_angle(ktheta1(iT2), ktheta2(iT2), kphi1(iT2, iP2), kphi2(iT2, iP2)) / (4 * pi)
                    p_theta1 = acos((cos(ktheta1(iT1)) + cos(ktheta2(iT1))) / 2)
                    p_theta2 = acos((cos(ktheta1(iT2)) + cos(ktheta2(iT2))) / 2)
                    p_phi1 = (kphi1(iT1, iP1) + kphi2(iT1, iP1)) / 2
                    p_phi2 = (kphi1(iT2, iP2) + kphi2(iT2, iP2)) / 2
                    call cartesian(1.d0, p_theta1, p_phi1, k1(1), k1(2), k1(3))
                    call cartesian(1.d0, p_theta2, p_phi2, k2(1), k2(2), k2(3))
                    do iW1 = 1, Nomega
                        omega1 = omega_mid(iW1)
                        n1 = n_tot * cell_bins(indR, indTheta, indPhi, iT1, iP1, iW1) / l_tot
                        do iW2 = 1, Nomega
                            omega2 = omega_mid(iW2)
                            n2 = n_tot * cell_bins(indR, indTheta, indPhi, iT2, iP2, iW2) / l_tot
                            if ((n1 .gt. 0) .and. (n2 .gt. 0)) then
                                rate = (n1 * n2 * sigma_c(omega1, omega2, k1, k2, B) * c_light) * vol  
                                rate_arr(iW1, iW2) = rate_arr(iW1, iW2) + rate 
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do
    return
    end function rate_omega_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! return total luminosity of produced pairs as a function of the higher photon energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function pair_lum_dist(indR, indTheta, indPhi) result(lum_arr)
    implicit none
    !input variables to select cell
    integer, intent(in) :: indR, indTheta, indPhi
    !output variable
    real*8, dimension(Nomega) :: lum_arr
    !local variables
    integer :: iT1, iT2, iP1, iP2, iW1, iW2, iL
    real*8 :: n1, n2, l_tot, n_tot, p_theta1, p_theta2, p_phi1, p_phi2, omega1, omega2
    real*8, dimension(3) :: k1, k2, B
    real*8 :: vol, rate, lum
    lum_arr = 0
    vol = cell_volume(indR, indTheta, indPhi)
    l_tot = cell_time(indR, indTheta, indPhi)
    n_tot = (l_tot / (c_light * cell_volume(indR, indTheta, indPhi))) * (N_emission / N_prop)
    if (n_tot .eq. 0) then
        return
    end if
    call cellB(indR, indTheta, indPhi, B)
    do iT1 = 1, Nk_theta
        do iP1 = 1, Nk_phi(iT1)
            do iT2 = 1, Nk_theta
                do iP2 = 1, Nk_phi(iT2)
                    p_theta1 = acos((cos(ktheta1(iT1)) + cos(ktheta2(iT1))) / 2)
                    p_theta2 = acos((cos(ktheta1(iT2)) + cos(ktheta2(iT2))) / 2)
                    p_phi1 = (kphi1(iT1, iP1) + kphi2(iT1, iP1)) / 2
                    p_phi2 = (kphi1(iT2, iP2) + kphi2(iT2, iP2)) / 2
                    call cartesian(1.d0, p_theta1, p_phi1, k1(1), k1(2), k1(3))
                    call cartesian(1.d0, p_theta2, p_phi2, k2(1), k2(2), k2(3))
                    do iW1 = 1, Nomega
                        omega1 = omega_mid(iW1)
                        n1 = n_tot * cell_bins(indR, indTheta, indPhi, iT1, iP1, iW1) / l_tot
                        do iW2 = 1, Nomega
                            omega2 = omega_mid(iW2)
                            n2 = n_tot * cell_bins(indR, indTheta, indPhi, iT2, iP2, iW2) / l_tot
                            if ((n1 .gt. 0) .and. (n2 .gt. 0)) then
                                rate = (n1 * n2 * sigma_c(omega1, omega2, k1, k2, B) * c_light) * vol  
                                lum = rate * (omega1 + omega2)
                                !pick higher energy photon
                                iL = max(iW1, iW2)
                                lum_arr(iL) = lum_arr(iL) + lum
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do
    return
    end function pair_lum_dist


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! creation rate for all spacial cells, done in parallel 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pair_lum_para(lum_arr)
    implicit none
    !output variable
    real*8, dimension(Nomega), intent(out) :: lum_arr 
    real*8, dimension(Nomega) :: lum_hold, lum_loss
    integer :: indR, indT, indP
    lum_arr = 0.d0
    do indR = 1, Nr+1
        do indT = 1, 2*Ntheta
            !$OMP PARALLEL PRIVATE(indP, lum_loss)
            lum_loss = 0.d0
            !$OMP DO
            do indP = 1, Nphi
                if (cell_volume(indR, indT, indP) > 0) then
                    lum_loss = lum_loss + pair_lum_dist(indR, indT, indP)
                end if
            end do
            !$OMP END DO
            !$OMP CRITICAL
            lum_arr = lum_arr + lum_loss
            !$OMP END CRITICAL
            !$OMP END PARALLEL
        end do
    end do
    end subroutine pair_lum_para
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end module photon

