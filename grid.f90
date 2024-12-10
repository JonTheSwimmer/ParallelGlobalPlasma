!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Propagate and collide gamma rays in magnetar magnetospheres
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module grid
    use parameters_base
    use routines
    use spectrum
    implicit none

    real*8, dimension(Nr), save :: rmax_mid, rmax1, rmax2, theta0_mid, theta01, theta02 
    integer, dimension(Nr), save :: cell_start
    real*8, dimension(Nphi), save :: phi1, phi2
    real*8, dimension(2*Ntheta), save :: theta_star1, theta_star2, theta_starmid
    real*8, dimension(Nr + 1, 2*Ntheta, Nphi), save :: cell_time, cell_volume, cell_length
    real*8, dimension(Nr + 1, 2*Ntheta, Nphi), save :: cell_create
    real*8, dimension(Nr + 1, 2*Ntheta, Nphi, Nomega), save :: cell_create_split, cell_loss_split
    real*8, dimension(Nr + 1, 2*Ntheta, Nphi, Nk_theta, Nk_phimax, Nomega), save :: cell_bins
    ! array for surface energy, with flux bundle, north/south, and phi section
    real*8, dimension(Nr+1, 2, Nphi), save :: surf_flux 
    integer(8), save :: N_prop, N_col
    integer(8), dimension(4), save:: stat_emit, stat_surf, stat_conv
    real*8, dimension(Nk_theta), save :: ktheta1, ktheta2
    integer, dimension(Nk_theta), save :: Nk_phi !number of phi divisions in each theta cell 
    real*8, dimension(Nk_theta, Nk_phimax), save :: kphi1, kphi2
    real*8 :: TS1, TS2 !dummy variables to update during cell volume integration



    contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Length along magnetic dipole field line, measured from theta0 to theta < pi/2-theta0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
    function L(theta0, theta)
    implicit none
    real*8, intent(in) :: theta0, theta
    real*8 :: L, x, x0

    x0 = sqrt(3.d0) * cos(theta0);  x = sqrt(3.d0) * cos(theta)
    L = ( log( (x0 + sqrt(1.d0 + x0**2))/(x + sqrt(1.d0 + x**2)) ) + &
           x0*sqrt(1.d0+x0**2) - x*sqrt(1.d0+x**2) ) / (2.d0*sqrt(3.d0)*(1 - (x0**2)/3.d0))

    return
    end function L


!!!!i!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! r_pole given theta_*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
    function r_pole(theta)
    implicit none
    real*8, intent(in) :: theta
    real*8 :: r_pole, x
    x = theta
    if (x .eq. pi/2) then
        r_pole = 1.d6
        return
    end if
    if (x .gt. pi/2) then
        x = pi - x
    end if
    r_pole = exp((1/cos(x)) - 1) 
    return 
    end function r_pole

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! intersection theta of a flux line and equipotential line, given thetaS (eq line) and R (flux line)
! assumes that flux and equipotential lines are in 0 < theta < theta_max < pi/2
! theta_max = theta_star1(Ntheta) < pi/2
! intersection function is (1-x^2)exp((1-x)/x) = r_pole / R and is decreasing on [0, 1]
! for x = cos(theta)
! must be used after theta lines are set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    function theta_int(thetaS, R)
    implicit none
    real*8 :: thetaS, R
    real*8 :: theta_int, x, x_high, x_low, y, y_test, tol
    logical :: flip
    tol = 1.d-12
    if (thetaS == pi/2.d0) then
        theta_int = pi/2.d0
        return
    end if
    flip = .false.
    if (thetaS > pi/2.d0) then
        flip = .true.
        thetaS = pi - thetaS
    end if
    y = r_pole(thetaS) / R
    x_low = 0; x_high = 1; x = (x_low + x_high) / 2.d0 !x = cos(theta)
    do while (((x_high - x_low) / x_low) > tol)
        y_test = (1 - x**2.d0) * exp((1-x)/x) 
        if (y_test < y) then
            x_high = x
        else if (y_test > y) then
            x_low = x
        else !if it somehow lands to be exactly equal
            theta_int = acos(x)
            return
        end if
        x = (x_high + x_low) / 2.d0
    end do
    theta_int = acos(x)
    if (flip) then
        theta_int = pi - theta_int
        thetaS = pi - thetaS
    end if
    return
    end function theta_int


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convert between grid coordinates (curly R, curly theta) -> spherical (R, theta)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine grid2sphere(curlyR, thetaS, r, theta)
    implicit none
    real*8, intent(in) :: curlyR, thetaS
    real*8, intent(out) :: r, theta
    theta = theta_int(thetaS, curlyR)
    r = curlyR * sin(theta)**2.d0
    return
    end subroutine grid2sphere

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! magnetic field in cartesian 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine Bfield(r, theta, phi, Bx, By, Bz)
    implicit none
    real*8, intent(in) :: r, theta, phi
    real*8, intent(out) :: Bx, By, Bz
    real*8 :: x, y, z
    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)
    Bx = B_dip * 3.d0*x*z / r**5.d0
    By = B_dip * 3.d0*y*z / r**5.d0
    Bz = B_dip * ((3.d0*z**2 / r**5.d0) - (1.d0 / r**3.d0))
    return
    end subroutine Bfield

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dA/dR to find volume of each cell
! dphi * int(dA/dR) = volume
! dependence on equipotential lines is inserted via rp1, rp2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine deriv_area(neq, R_in, X, dXdr)
    integer, intent(in) :: neq
    real*8, intent(in) :: R_in, X(neq)
    real*8, intent(out) :: dXdr(neq)
    real*8 :: theta1, theta2, theta_min, theta_max
    theta_min = asin(R_in**(-0.5d0))
    theta_max = asin((Rmax / R_in)**(0.5d0))
    !filter theta bounds to grid limits
    theta1 = max(theta_min, min(theta_max, theta_int(TS1, R_in)))
    theta2 = max(theta_min, min(theta_max, theta_int(TS2, R_in)))
    !print *, theta1, theta2
    dXdr(1) = (R_in**(2.d0) / 2240) * (1225*(cos(theta1) - cos(theta2)) &
                - 245*(cos(3*theta1) - cos(3*theta2)) + 49* (cos(5*theta1) - cos(5*theta2)) &
                - 5*(cos(7*theta1) - cos(7*theta2)))

    return
    end subroutine deriv_area

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! volume calculation for polar caps
! integrate from the bounding R line to the minimum of theta02 (deriv_cap1)
! then use deriv_cap2 to integrate the rest as a conical-ish slice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine deriv_cap1(neq, R_in, X, dXdr)
    integer, intent(in) :: neq
    real*8, intent(in) :: R_in, X(neq)
    real*8, intent(out) :: dXdr(neq)
    real*8 :: theta1, theta2, theta_min, theta_max
    theta_min = theta02(Nr)
    theta_max = asin((Rmax / R_in)**(0.5d0))
    
    theta1 = max(theta_min, min(theta_max, theta_int(TS1, R_in)))
    theta2 = max(theta_min, min(theta_max, theta_int(TS2, R_in)))
    !print *, theta1, theta2
    dXdr(1) = (R_in**(2.d0) / 2240) * (1225*(cos(theta1) - cos(theta2)) &
                - 245*(cos(3*theta1) - cos(3*theta2)) + 49* (cos(5*theta1) - cos(5*theta2)) &
                - 5*(cos(7*theta1) - cos(7*theta2)))
    return
    end subroutine deriv_cap1
    
    subroutine deriv_cap2(neq, T_in, X, dXdT)
    integer, intent(in) :: neq
    real*8, intent(in) :: T_in, X(neq)
    real*8, intent(out) :: dXdT(neq)
    real*8 :: r1, r2
    r1 = r_pole(TS1) * exp(1.d0 - (1/cos(T_in)))
    r2 = r_pole(TS2) * exp(1.d0 - (1/cos(T_in)))
    !filter rp1, rp2 to grid limits
    r1 = max(1.d0, min(r1, Rmax))
    r2 = max(1.d0, min(r2, Rmax))

    dXdT = (sin(T_in) / 3.d0) * (r2**3.d0 - r1**3.d0)
    
    return
    end subroutine deriv_cap2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! r_mid == mid-radius of magnetic flux bundle at theta = pi/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine create_grid
    implicit none
    integer :: ir, itheta, iphi, neq, iomega
    real*8 :: dlogrmax_mid, factor, dtheta, dphi, theta2_step
    logical :: array_inc
    real*8 :: A(1), R1, R2, rp1, rp2, TM !theta_min
    real*8 :: lPhi, lR, lTS, theta, r
    real*8 :: theta_step, phi_step, omega_step
    integer:: stat
    !Two 
    dlogrmax_mid = (log10(Rmax-1.d0) - log10(dRmax)) / (Nr-Nextra - 1)
    factor = 1.d1**(0.5d0*dlogrmax_mid)
! Radius of equatorial mid-point of flux bundle, in units of R_star
    do ir = 1, Nr
         rmax_mid(ir) = 1.d0 + 1.d1**(log10(dRmax) + dlogrmax_mid*(ir-1))
         theta0_mid(ir) = asin(rmax_mid(ir)**(-0.5))
    end do    
! r1 = inner radial boundary;  r2 = outer radial boundary 
    rmax1(1) = 1.d0;  rmax2(1) = 1.d0 + (rmax_mid(1)-1.d0)*factor
    theta01(1) = asin(rmax1(1)**(-0.5)); theta02(1) = asin(rmax2(1)**(-0.5))
    do ir = 2, Nr
        rmax1(ir) = rmax2(ir-1)
        rmax2(ir) = 1.d0 + (rmax_mid(ir)-1.d0) * factor
        theta01(ir) = asin(rmax1(ir)**(-0.5));  theta02(ir) = asin(rmax2(ir)**(-0.5))
    end do
    
!   cos(theta_*) splits in grid
    theta2_step = ((pi / 2.d0)**2.d0) / Ntheta
    do itheta = 1, Ntheta
        theta_star1(itheta) = ((itheta - 1) * theta2_step)**(0.5d0)
        theta_starmid(itheta) = ((itheta-0.5) * theta2_step)**(0.5d0)
        theta_star2(itheta) = (itheta * theta2_step)**(0.5d0)
    end do
    do itheta = Ntheta + 1, 2*Ntheta
        theta_star1(itheta) = theta_star2(itheta-1)
        theta_starmid(itheta) = pi - theta_starmid(2*Ntheta+1-itheta)
        theta_star2(itheta) = pi - theta_star1(2*Ntheta + 1 - itheta)
    end do

!   find first equipotential line that bounds a given flux line
    array_inc = .true.
    do ir = 1, Nr
        cell_start(ir) = bisect(theta_star1, theta02(ir), 2*Ntheta, array_inc)
    end do
!   phi splits in grid
    dphi = twopi / Nphi
    do iphi = 1, Nphi

        phi1(iphi) = (iphi-1)*dphi
        phi2(iphi) = phi1(iphi) + dphi

    end do
!   integrate cell volume for all cells
!   non-polar cells first because they're easier
    cell_volume(:, :, :) = 0
    do ir = 1, Nr
        R1 = rmax1(ir)
        R2 = rmax2(ir)
        do itheta = cell_start(ir), Ntheta
            TS1 = theta_star1(itheta) 
            TS2 = theta_star2(itheta)
            A(1) = 0
            call integrate_ode(1, A, R1, R2, deriv_area)
            cell_volume(ir, itheta, :) = A(1) * dphi
        end do
        do itheta = Ntheta+1, 2*Ntheta
            cell_volume(ir, itheta, :) = cell_volume(ir, 2*Ntheta+1 - itheta, :)
        end do
    end do
!   integrate volume for polar cells
    TM = theta02(Nr)
    do itheta = 1, Ntheta
        TS1 = theta_star1(itheta)
        TS2 = theta_star2(itheta)
        if (TS1 < acos(cos(TM) / (cos(TM) *log(Rmax) + 1))) then
            R1 = rmax2(Nr)
            R2 = r_pole(TS2) * exp(1 - 1/cos(theta02(Nr))) / sin(theta02(Nr))**(2.d0)
            A(1) = 0
            call integrate_ode(1, A, R1, R2, deriv_cap1)
            call integrate_ode(1, A, 0.d0, theta02(Nr), deriv_cap2)
            cell_volume(Nr+1, itheta, 1) = A(1) * twopi
            cell_volume(Nr+1, 2*Ntheta+1 - itheta, 1) = cell_volume(Nr+1, itheta, 1)
        end if
    end do
!   determine characteristic cell length for each cell
    do ir = 1, Nr
        R1 = rmax1(ir)
        R2 = rmax2(ir)
        do itheta = cell_start(ir), Ntheta
            TS1 = theta_star1(itheta)
            TS2 = theta_star2(itheta)
            theta = theta_int(TS2, R2)
            r = R2 * sin(theta)**2
            lPhi = r * sin(theta) * (twopi/Nphi)
            lR = R2 - R1
            lTS = r * (TS2 - TS1)
            cell_length(ir, itheta, :) = (lPhi**2.d0 + lTS**2.d0 + lR**2.d0)**0.5d0
        end do
        do itheta = Ntheta+1, 2*Ntheta
            cell_length(ir, itheta, :) = cell_length(ir, 2*Ntheta+1 - itheta, :)
        end do
    end do
    do itheta = 1, Ntheta
        TS1 = theta_star1(itheta)
        TS2 = theta_star2(itheta)
        if (TS1 < acos(cos(TM) / (cos(TM) *log(Rmax) + 1))) then
            lTS = r_pole(TS2) - r_pole(TS1)
            lR = r_pole(TS2)*theta02(Nr)*2
            cell_length(Nr+1, itheta, 1) = (lR**2.d0 + lTS**2.d0)**0.5d0
            cell_length(Nr+1, 2*Ntheta+1-itheta, 1) = cell_length(Nr+1, itheta, 1)
        end if
    end do

!   Set angular cell boundaries for photon trajectories
    theta_step = pi / (Nk_theta)
    do itheta = 1, Nk_theta
        ktheta1(itheta) = theta_step * (itheta-1) 
        ktheta2(itheta) = theta_step * (itheta) 
        Nk_phi(itheta) = min(Nk_phimax, floor(Nk_phimax * sin(0.5d0 * (ktheta1(itheta) + ktheta2(itheta))) + 1))
        phi_step = twopi / Nk_phi(itheta)
        do iphi = 1, Nk_phi(itheta)
            kphi1(itheta, iphi) = (iphi - 1) * phi_step
            kphi2(itheta, iphi) = iphi * phi_step
        end do
    end do
    
    cell_time= 0
    cell_bins = 0
    cell_create = 0
    N_prop = 0
    N_col = 0
    surf_flux = 0
    stat_emit = 0; stat_surf = 0; stat_conv = 0
    call init_spectrum()
    return
    end subroutine create_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Output relevant parameters for the grid    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine output_grid(runnum)
    implicit none
    integer :: iR, iT
    integer, intent(in) :: runnum
    character(len=8) :: num ! format descriptor
    write (num, "(I0)") runnum
    call file_open(2, trim(num)//'_parameters.dat')
101 format(*(g0, :, ", "))
    write (2, 101) Nr, Nphi, Ntheta, Nk_theta, Nk_phimax, Nomega, Rmax
    close (2)
    call file_open(3, trim(num)//'_theta01.dat')
    call file_open(4, trim(num)//'_theta02.dat')
    call file_open(5, trim(num)//'_theta0_mid.dat')
    write(3, 101) theta01
    write(4, 101) theta02
    write(5, 101) theta0_mid
    close(3)
    close(4)
    close(5)
    call file_open(8, trim(num)//'_theta_star1.dat')
    call file_open(9, trim(num)//'_theta_star2.dat')
    call file_open(10, trim(num)//'_theta_starmid.dat')
    write(8, 101) theta_star1
    write(9, 101) theta_star2
    write(10, 101) theta_starmid
    close(8)
    close(9)
    close(10)
    call file_open(11, trim(num)//'_cell_vol.dat')
    do iR = 1, Nr+1
        do iT = 1, 2*Ntheta
            write(11, 101) cell_volume(iR, iT, :)
        end do
    end do
    close(11)
    call file_open(12, trim(num)//"_Nk_phi.dat")
    call file_open(13, trim(num)//"_ktheta1.dat")
    call file_open(14, trim(num)//"_ktheta2.dat")
    write (12, 101) Nk_phi
    write (13, 101) ktheta1
    write (14, 101) ktheta2
    close (12)
    close (13)
    close (14)
    call file_open(15, trim(num)//"_kphi1.dat")
    call file_open(16, trim(num)//"_kphi2.dat")
    do iR = 1, Nk_theta
        write (15, 101) kphi1(iR, :)
        write (16, 101) kphi2(iR, :)
    end do
    close (15)
    close (16)
    call file_open(18, trim(num)//"_omega01.dat")
    call file_open(19, trim(num)//"_omega02.dat")
    write (18, 101) omega01
    write (19, 101) omega02
    close(18); close(19)
    return
    end subroutine output_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Output relevant parameters for the grid    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine output_data(runnum)
    implicit none
    integer, intent(in) :: runnum
    integer :: iR, iT, iP, iK, iW
    character(len=8) :: num ! format descriptor
    write (num, "(I0)") runnum
    call file_open(2, 'cell_time'//trim(num)//'.dat')
    call file_open(3, 'cell_bins'//trim(num)//'.dat')
    call file_open(5, 'surf_flux'//trim(num)//'.dat')
101 format(*(g0, :, ", "))
    do iR = 1, Nr+1
        write(5, 101) surf_flux(iR, 1, :)
        write(5, 101) surf_flux(iR, 2, :)
        do iT = 1, 2*Ntheta
            write (2, 101) cell_time(iR, iT, :)
            do iP = 1, Nphi
                do iK = 1, Nk_theta
                    do iW = 1, Nk_phimax
                        write (3, 101) cell_bins(iR, iT, iP, iK, iW, :)
                    end do
                end do
            end do
        end do
    end do
    close (2)
    close (3)
    close (5)
    call file_open(4, 'N_prop'//trim(num)//'.dat')
    write (4, 101) N_prop
    close (4)
    call file_open(20, 'emit_stats'//trim(num)//'.dat')
    write (20, 101) stat_emit
    write (20, 101) stat_surf
    write (20, 101) stat_conv
    close (20)
    return
    end subroutine output_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Load data for the grid    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine load_data(runnum)
    implicit none
    integer, intent(in) :: runnum
    character(len=8) :: num
    integer :: iR, iT, iP, iK, iW
    real*8, dimension(Nphi) :: time_data, flux_data
    real*8, dimension(Nomega) :: bin_data
    real*8 :: N_p
    integer, dimension(4):: emit_data, surf_data, conv_data
    write (num, "(I0)") runnum
    open(2, file = 'data/cell_time'//trim(num)//'.dat')
    open(3, file = 'data/cell_bins'//trim(num)//'.dat')
    open(5, file = 'data/surf_flux'//trim(num)//'.dat')
    do iR = 1, (Nr + 1)
        read(5, *) flux_data
        surf_flux(iR, 1, :) = surf_flux(iR, 1, :) + flux_data
        read(5, *) flux_data
        surf_flux(iR, 2, :) = surf_flux(iR, 2, :) + flux_data
        do iT = 1, 2*Ntheta
            read(2, *) time_data(:)
            cell_time(iR, iT, :) = cell_time(iR, iT, :) + time_data(:)
            do iP = 1, Nphi
                do iK = 1, Nk_theta
                    do iW = 1, Nk_phimax
                        read(3, *) bin_data(:)
                        cell_bins(iR, iT, iP, iK, iW, :) = cell_bins(iR, iT, iP, iK, iW, :) + bin_data(:)
                    end do
                end do
            end do
        end do
    end do
    close(2)
    close(3)
    close(5)
    open(4, file = 'data/N_prop'//trim(num)//'.dat')
    read(4, *) N_p
    N_prop = N_prop + N_p
    close(4)
    open(20, file = 'data/emit_stats'//trim(num)//'.dat')
    read (20, *) emit_data
    stat_emit = emit_data
    read (20, *) surf_data
    stat_surf = surf_data
    read (20, *) conv_data
    stat_conv = conv_data
    close (20)
    end subroutine load_data
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Load data to private grid variables   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine priv_load(runnum, priv_flux, priv_time, priv_bins, priv_N)
    implicit none
    !input variable
    integer, intent(in) :: runnum
    !output variables
    real*8, dimension(Nr+1, 2, Nphi) :: priv_flux
    real*8, dimension(Nr+1, 2*Ntheta, Nphi) :: priv_time
    real*8, dimension(Nr+1, 2*Ntheta, Nphi, Nk_theta, Nk_phimax, Nomega) :: priv_bins
    integer :: priv_N
    !local variables
    character(len=8) :: num
    integer :: iR, iT, iP, iK, iW
    real*8, dimension(Nphi) :: time_data, flux_data
    real*8, dimension(Nomega) :: bin_data
    real*8 :: N_p
    write (num, "(I0)") runnum
    open(2, file = 'data/cell_time'//trim(num)//'.dat')
    open(3, file = 'data/cell_bins'//trim(num)//'.dat')
    open(5, file = 'data/surf_flux'//trim(num)//'.dat')
    do iR = 1, (Nr + 1)
        read(5, *) flux_data
        priv_flux(iR, 1, :) = priv_flux(iR, 1, :) + flux_data
        read(5, *) flux_data
        priv_flux(iR, 2, :) = priv_flux(iR, 2, :) + flux_data
        do iT = 1, 2*Ntheta
            read(2, *) time_data(:)
            priv_time(iR, iT, :) = priv_time(iR, iT, :) + time_data(:)
            do iP = 1, Nphi
                do iK = 1, Nk_theta
                    do iW = 1, Nk_phimax
                        read(3, *) bin_data(:)
                        priv_bins(iR, iT, iP, iK, iW, :) = priv_bins(iR, iT, iP, iK, iW, :)+bin_data(:)
                    end do
                end do
            end do
        end do
    end do
    close(2)
    close(3)
    close(5)
    open(4, file = 'data/N_prop'//trim(num)//'.dat')
    read(4, *) N_p
    priv_N = priv_N + N_p
    close(4)
    end subroutine priv_load

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Load data for the grid, if data had no omega data 
!   Loads data as if all photons were monochromatic with energy ~omega_mono
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine load_monodata(runnum)
    implicit none
    integer, intent(in) :: runnum
    character(len=8) :: num
    integer :: iR, iT, iP, iK, iOmega 
    real*8, dimension(Nphi) :: time_data, flux_data
    real*8, dimension(Nk_phimax) :: bin_data
    real*8 :: N_p
    iOmega = bisect(omega01, omega_mono, Nomega, .true.)
    !omega_mid(iOmega) = omega_mono
    write (num, "(I0)") runnum
    open(2, file = 'data/cell_time'//trim(num)//'.dat')
    open(3, file = 'data/cell_bins'//trim(num)//'.dat')
    open(5, file = 'data/surf_flux'//trim(num)//'.dat')
    do iR = 1, (Nr + 1)
        read(5, *) flux_data
        surf_flux(iR, 1, :) = surf_flux(iR, 1, :) + flux_data
        read(5, *) flux_data
        surf_flux(iR, 2, :) = surf_flux(iR, 2, :) + flux_data
        do iT = 1, 2*Ntheta
            read(2, *) time_data(:)
            cell_time(iR, iT, :) = cell_time(iR, iT, :) + time_data(:)
            do iP = 1, Nphi
                do iK = 1, Nk_theta
                    read(3, *) bin_data(:)
                    cell_bins(iR, iT, iP, iK, :, iOmega) = cell_bins(iR, iT, iP, iK, :, iOmega) + bin_data(:)
                end do
            end do
        end do
    end do
    close(2)
    close(3)
    close(5)
    open(4, file = 'data/N_prop'//trim(num)//'.dat')
    read(4, *) N_p
    N_prop = N_prop + N_p
    close(4)
    end subroutine load_monodata
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine which grid cell a coordinate is in
! indTheta0 is Nr+1 for polar cells, as they extend past the furthest theta0 line by definition
! indPhi will always be 1 for polar cells
! indThetaS is determined by the equipotential cell grid divisions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine place(R, theta, phi, indTheta0, indThetaS, indPhi)
    implicit none
    real*8, intent(in) :: R, theta, phi
    integer, intent(out) :: indThetaS, indTheta0, indPhi
    real*8 :: theta0, thetaS, thetaflip
    logical :: array_increasing
    integer :: Nmax

    if ((R < 1) .or. (R > Rmax)) then
        indThetaS = -1;  indTheta0 = -1; indPhi = -1;
        return
    endif
!   Assumes that phi divisions are uniformly split
    indPhi = floor(mod(phi, twopi) * Nphi / twopi) + 1
    theta0 = asin((sin(theta)**2 / R)**0.5)
!   check if theta0 is in either of the polar caps
    if (theta0 < theta02(Nr)) then
        indTheta0 = Nr + 1; indPhi = 1 
    else
        array_increasing = .false.
        indTheta0 = bisect(theta01, theta0, Nr, array_increasing)
    end if
!   determine which equipotential line bounds the point, i.e indThetaS    
    if (theta < pi/2.d0) then
        thetaS = acos(cos(theta) / (1.d0 + cos(theta) * log(R)))
    else
        thetaflip = pi - theta
        thetaS = pi - acos(cos(thetaflip) / (1.d0 + cos(thetaflip) * log(R)))
    end if
    array_increasing = .true.
    indThetaS = bisect(theta_star1, thetaS, Ntheta*2, array_increasing) 
    return
    end subroutine place
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine which surface cell a photon hits 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine place_surf(theta, phi, indTheta0, indNS, indPhi)
    implicit none
    real*8, intent(in) :: theta, phi
    integer, intent(out) :: indTheta0, indNS, indPhi
    real*8 :: thetaflip
    logical :: array_inc
    indNS = 1
    thetaflip = theta
    if (theta .gt. pi/2) then
        indNS = 2 !southern hemisphere
        thetaflip = pi - theta
    end if
    if (thetaflip < theta02(Nr)) then
        indTheta0 = Nr+1; indPhi = 1
    else
        array_inc = .false.
        indTheta0 = bisect(theta01, thetaflip, Nr, array_inc)
        indPhi = floor(mod(phi, twopi) * Nphi / twopi) + 1
    end if

    return
    end subroutine place_surf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine which angular cell and frequency bin a photon is in
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine place_k(k, omega, indTheta, indPhi, indOmega)
    implicit none
    real*8, intent(in) :: k(3), omega !cartesian direction of photon momentum
    integer, intent(out) :: indTheta, indPhi, indOmega
    real*8 :: R, theta, phi
    logical :: array_increasing
    call spherical(k(1), k(2), k(3), R, theta, phi)
    indTheta = bisect(ktheta1, theta, Nk_theta, .true.)
    associate (N => Nk_phi(indTheta))
    indPhi = bisect(kphi1(indTheta, 1:N), phi, N, .true.)
    end associate
    indOmega = bisect(omega01, omega, Nomega, .true.)
    return
    end subroutine place_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check that all indices are within bounds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function inBounds(indT0, indTS, indPhi, indKT, indKP)
    implicit none
    logical :: inBounds, c1, c2, c3, c4, c5
    integer, intent(in) :: indT0, indTS, indPhi, indKT, indKP
    !check that indT0 is in bounds
    c1 = ((indT0 .ge. 1) .and. (indT0 .le. (Nr+1)))
    c2 = ((indTS .ge. 1) .and. (indTS .le. 2*Ntheta))
    c3 = ((indPhi .ge. 1) .and. (indPhi .le. Nphi))
    c4 = ((indKT .ge. 1) .and. (indKT .le. Nk_theta))
    c5 = ((indKP .ge. 1) .and. (indKP .le. Nk_phimax))
    inBounds = (c1 .and. c2 .and. c3 .and. c4 .and. c5)
    return
    end function inBounds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! B-field for a given cell 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cellB(indR, indTheta, indP, B)
    integer :: indR, indTheta, indP
    real*8, dimension(3) :: B
    real*8 :: curlyR, curlyT, cellR, cellT, cellP
    if (indR .lt. Nr+1) then
        curlyR = (rmax1(indR) + rmax2(indR))/2.d0
        curlyT = (theta_star1(indTheta) + theta_star2(indTheta))/2.d0
        cellP = (phi1(indP) + phi2(indP))/2.d0
        call grid2sphere(curlyR, curlyT, cellR, cellT)
    else
        if (indTheta .le. Ntheta) then
            cellT = 0.d0
        else
            cellT = pi
        end if
        cellP = 0.d0
        cellR = (r_pole(theta_star1(indTheta)) + r_pole(theta_star2(indTheta))) / 2.d0
    end if
    call Bfield(cellR, cellT, cellP, B(1), B(2), B(3))
    end subroutine cellB
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end module grid

