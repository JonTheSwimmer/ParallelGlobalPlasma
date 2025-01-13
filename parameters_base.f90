!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module parameters_base
    implicit none

    real*8, parameter :: pi = 3.1415926535897932d0,  twopi = 2.d0 * pi

! Parameters for magnetospheric grid
! Maximum, minimum radius of equatorial mid-point of cell
    real*8, parameter :: Rmax = 2.d0**5
    real*8, parameter :: dRmax = 0.1d0
! Distribute magnetic flux bundles uniformly in lg(r-Rstar)
! Number of bundles
    integer, parameter :: Nextra = 4 !number of bundles that exceed Rmax
    integer, parameter :: Nr = 16 + Nextra
! Number of grid cells in phi
    integer, parameter :: Nphi = 2**6
! Number of grid cells in theta
    integer, parameter :: Ntheta = 2**5
! Number of total possible  cells
    integer, parameter :: Ncells = (Nphi+1)*2*Ntheta*Nr
! Number of theta cells for tracking photon direction
    integer, parameter :: Nk_theta = 2**3
! Max number of phi cells in one theta cell for photon direction
    integer, parameter :: Nk_phimax = 2**3

! Magnetar parameters
    real*8, parameter :: M_NS = 1.4d0 * 2.d33 ! mass in grams
    real*8, parameter :: R_NS = 1.d6  ! radius in cm
    real*8, parameter :: rg = 2.d0 * 7.4213758738776d-29 * M_NS / R_NS !gravitational radius divided by magnetar radius
    real*8, parameter :: b_crit = (1.d0- rg)**(-0.5) ! critical impact parameter
! Physical parameters
    real*8, parameter :: c_light = 2.99792458d10 / R_NS !speed of light in units of R_NS/second
    real*8, parameter :: r_elec = 2.8179d-13 / R_NS !radius of electron in units of R_NS
    real*8, parameter :: m_e = 8.1871057769d-7 !electron mass in ergs/c^2
    real*8, parameter :: keV = 1.60218d-9 !keV to ergs
    real*8, parameter :: z_max = sqrt(1 - rg/Rmax) / sqrt(1 - rg) !max gravitational redshift, > 1
    real*8, parameter :: sigma_SB = 5.670374419d-5 !stefan-boltzmann constant in cgs
! photon and plasma parameters
    real*8, parameter :: luminosity = 1.d35 !luminosity of emission structure in ergs/s
    real*8, parameter :: omega_mono = 2.d0 * m_e !energy for monochromatic photon, in units of ergs
    !real*8, parameter :: B_dip = 5.d0 !dipole magnetic moment in units of B_q
    real*8, save :: B_dip 
    real*8, parameter :: beta_e = 0.7 !velocity of pair plasma in flux bundles
    real*8, parameter :: beta_surf = 0.2 !velocity of pair plasma after surface screening
    real*8, parameter :: omega_min = 1.5d1*keV !minimum energy in spectrum
    real*8, parameter :: omega_max = 2.5d0*m_e !max energy in spectrum
    real*8, parameter :: omega_binmin = omega_min / z_max !min energy in bins, due to gravitational redshift
    real*8, parameter :: omega_binmax = 3.d0*m_e !omega_max * z_max !max energy in bins, due to gravitational blueshift
    real*8, parameter :: omega_cut = 1.5d0 * m_e !energy for cutoff in spectrum
    real*8, parameter :: N_spec = 4.d0 !power law in exponential for spectrum
    integer, parameter :: Nomega = 2**7
    integer, parameter :: N_ospline = 1024
    real*8, parameter :: omega_out = 1 * keV !energy of x-rays that resonantly scatter plasma in outer grid
    real*8, parameter :: surf_eff = 0.5d0 !efficiency of photon surface absorption to determine heating
! Parameters for trajectory
    integer, parameter :: N_spline = 2**11 ! number of splines to interpolate with
    integer, parameter :: N_spoints = 2**10 ! number of points for each spline
    real*8, parameter  :: b_max = Rmax * 1.1 ! maximum impact parameter to sample
    real*8, parameter :: logb_min = -3.d0 ! minimum impact parameter
    real*8, parameter  :: R_limit = Rmax*1.5d0 ! maximum radius of trajectory to travel
    integer, parameter :: N_smax = N_spoints*3 ! max number of points in a photon trajectory
    real*8, parameter :: R_inner = 0.7d0 !minimum radius to consider in photon trajectory
    real*8, parameter :: b_inner = R_inner * (1 - rg/R_inner)**(-0.5d0) !impact parameter for trajectory skimming R_inner
    real*8, parameter :: dpsi = 0.05 !angle of second point to avoid singularity at R_min

! Integration tolerance
    real*8, parameter :: tolerance = 1.d-12

! Parameters for parallel propagation
    integer, parameter :: N_para_list = 2**11
    integer, save :: N_threads
    integer, parameter :: N_threadmax = 128 !max number of threads

    end module parameters_base


