!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    module parameters
    use parameters_base
    implicit none


! Parameters for emission
    real*8, parameter :: theta_out = pi/4 !theta0 for emission surface assuming isotropic emission on dipole field bundle
    real*8, parameter :: theta_in = 0.9088084183150976d0 !inner theta0 for bundle
    real*8, parameter :: E_phi0 =-1.d0*pi/4!-1.d0*pi/4! !phi coord of right side of emission structure
    real*8, parameter :: E_phi1 = 1.d0*pi/4!1.d0*pi/4 !phi coord of left side of emission structure
    integer, parameter :: N_espline = 2500 !number of points in l-theta spline on the emisison surface


    end module parameters


