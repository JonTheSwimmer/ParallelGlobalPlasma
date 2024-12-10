!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Functions for the photon spectrum 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module spectrum
    use parameters_base
    use routines
    implicit none

    real*8 :: spec_norm !normalization constant for spectrum
    real*8 :: spec_weights(Nomega) !integrated weight for each bin during initial spectrum
    real*8, dimension(Nomega), save :: omega01, omega02, omega_mid
    real*8, dimension(N_ospline), save :: spec_omega_spl, spec_num_spl, spec_omega2_spl
    contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dL/dE without normalization 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    function base_spec(omega)
    real*8 :: omega, base_spec
    base_spec = exp(-((omega/omega_cut)**(N_spec)))
    return
    end function base_spec

    subroutine deriv_spec(neq, w_in, f, dfdw)
    integer, intent(in) :: neq
    real*8, intent(in) :: w_in, f(neq)
    real*8, intent(out) :: dfdw(neq)
    dfdw(1) = base_spec(w_in)
    return 
    end subroutine deriv_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! find normalization for spectrum, and set up bins
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine init_spectrum
    real*8 :: f(1), A(1)
    real*8 :: omega_step, weight_norm, omega, dN0(1), dN1(1)
    integer :: iomega
    f(1) = 0.d0
    call integrate_dlsode(1, f, omega_min, omega_max, deriv_spec)
    spec_norm = 1.d0/f(1)
    ! set photon energy boundaries for bins
    omega_step = (log10(omega_binmax) - log10(omega_binmin))/Nomega
    do iomega = 1, Nomega
        omega01(iomega) = 1.d1**(log10(omega_binmin) + (iomega-1) * omega_step) 
        omega02(iomega) = 1.d1**(log10(omega_binmin) + (iomega)*omega_step)
        omega_mid(iomega) = omega_avg(omega01(iomega), omega02(iomega))
    end do
    A(1) = 0.d0
    call integrate_dlsode(1, A, omega_min, omega_max, deriv_weight)
    weight_norm = A(1)
    do iomega = 1, Nomega
        A(1) = 0.d0
        call integrate_ode(1, A, omega01(iomega), omega02(iomega), deriv_weight)
        spec_weights(iomega) = A(1) / weight_norm
    end do
    !set up spline for converting a uniform distribution in (0, 1) to photon spectrum
    omega_step = (log10(omega_max) - log10(omega_min)) / (N_ospline - 1)
    spec_omega_spl(1) = omega_min
    spec_num_spl(1) = 0.d0
    A(1) = 0.d0
    do iomega = 2, N_ospline
        omega = 10**(log10(omega_min) + omega_step * (iomega - 1))
        spec_omega_spl(iomega) = omega
        call integrate_ode(1, A, spec_omega_spl(iomega-1), spec_omega_spl(iomega), deriv_weight)
        spec_num_spl(iomega) = A(1)
    end do
    spec_num_spl = spec_num_spl / A(1)
    call deriv_weight(1, omega_min, A, dN0)
    call deriv_weight(1, omega_max, A, dN1)
    dN0 = dN0 / A
    dN1 = dN1 / A
    call spline(spec_num_spl, spec_omega_spl, N_ospline, 1/dN0(1), 1/dN1(1), spec_omega2_spl)
    return
    end subroutine init_spectrum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! normalized dL/dE, assuming init_spectrum has already been run 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function photon_spec(omega)
    real*8 :: photon_spec, omega
    photon_spec = spec_norm * base_spec(omega)
    return
    end function photon_spec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! return weighted average energy between omega1 and omega2 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine deriv_weight(neq, w_in, f, dfdw)
    integer, intent(in) :: neq
    real*8, intent(in) :: w_in, f(neq)
    real*8, intent(out) :: dfdw(neq)
    dfdw(1) = base_spec(w_in)/w_in
    return
    end subroutine deriv_weight

    function omega_avg(omega1, omega2)
    real*8, intent(in) :: omega1, omega2
    real*8 :: omega_avg, o_step, num, denom
    real*8, dimension(1000) :: o_arr, w_arr, f_arr
    integer:: iW
    !f(1) = 0.d0
    !g(1) = 0.d0 
    
    !call integrate_dlsode(1, g, omega1, omega2, deriv_spec)
    !call integrate_dlsode(1, f, omega1, omega2, deriv_weight)
    !omega_avg = f(1) / g(1)
    o_step = (omega2 - omega1) / (999)
    do iW = 1, 1000
        o_arr(iW) = omega1 + o_step * (iW - 1)
        w_arr(iW) = base_spec(o_arr(iW))
        f_arr(iW) = base_spec(o_arr(iW))/o_arr(iW)
    end do
    num = 0.d0; denom = 0.d0
    do iW = 1, 999
        num = num + (w_arr(iW) + w_arr(iW+1)) * (o_arr(iW+1) - o_arr(iW)) / 2
        denom = denom + (f_arr(iW) + f_arr(iW+1)) * (o_arr(iW+1) - o_arr(iW))/2
    end do
    omega_avg = num/denom
    return
    end function omega_avg
    
    function omega_avg2(omega1, omega2)
    real*8, intent(in) :: omega1, omega2
    real*8 :: omega_avg2, o_step, num, denom
    real*8, dimension(1) :: f, g
    f(1) = 0.d0
    g(1) = 0.d0 
    
    call integrate_ode(1, g, omega1, omega2, deriv_spec)
    call integrate_ode(1, f, omega1, omega2, deriv_weight)
    omega_avg2 = g(1) / f(1)
    return
    end function omega_avg2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! return total luminosity between omega1 and omega2 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function lum_total(omega1, omega2) result(lum)
    !input variables
    real*8, intent(in) :: omega1, omega2
    !output variable
    real*8 :: lum
    !local variables
    real*8, dimension(1) :: lum_int
    lum_int = 0.d0
    call integrate_ode(1, lum_int, omega1, omega2, deriv_spec)
    lum = lum_int(1) * spec_norm  * luminosity
    return
    end function lum_total

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! return total photon emission rate between omega1 and omega2
! int_omega1^omega2 1/E dL/dE dE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function rate_total(omega1, omega2) result(rate)
    !input variables 
    real*8, intent(in) :: omega1, omega2
    !output variable
    real*8 :: rate
    !local variables
    real*8, dimension(1) :: rate_int
    rate_int = 0.d0
    call integrate_ode(1, rate_int, omega1, omega2, deriv_weight)
    rate = rate_int(1) * spec_norm  * luminosity
    end function rate_total

    end module spectrum
