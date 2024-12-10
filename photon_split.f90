!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!new photon module to determine the pair creation rate and photon consumption rate
!the old photon module was getting too long....
!pair creation rate and photon consumption rate are now split into Nomega bins
!each bin is the PC/PC rate determined by the higher photon energy for cuts

module photon_split
use parameters_base
use routines
use spectrum
use grid
use OMP_LIB
use photon

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! pair creation rate for a given spacial cell
! also determines photon consumption rate which is basically the same, except doubled when:
!    iT1 = iT2, iP1 = iP2, and iW1 = iW2
! N_pair, N_consumption, and E_rate are labelled into Nomega bins by higher photon energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rate_split(indR, indTheta, indPhi, N_pair, N_consumption, E_rate)
    implicit none
    !input variables
    integer, intent(in) :: indR, indTheta, indPhi
    !output variables
    real*8, intent(out), dimension(Nomega) :: N_pair, N_consumption, E_rate
    !local variables
    integer :: iT1, iT2, iP1, iP2, iW1, iW2, iWmax
    real*8 :: n1, n2, l_tot, n_tot, p_theta1, p_theta2, p_phi1, p_phi2, omega1, omega2
    real*8 :: Ndot
    real*8, dimension(3) :: k1, k2, B
    N_pair = 0.d0; N_consumption = 0.d0; E_rate = 0.d0
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
                                iWmax = max(iW1, iW2) !find the higher energy photon
                                N_pair(iWmax) = N_pair(iWmax) + Ndot 
                                E_rate(iWmax) = E_rate(iWmax) + Ndot*(omega1 + omega2)
                                if ((iW1 .eq. iW2) .and. (iT1 .eq. iT2) .and. (iP1 .eq. iP2)) then
                                    N_consumption(iWmax) = N_consumption(iWmax) + Ndot*2.d0
                                else 
                                    N_consumption(iWmax) = N_consumption(iWmax) + Ndot
                                end if
                            end if
                        end do
                    end do
                end do
            end do
        end do
    end do
    N_pair = N_pair * cell_volume(indR, indTheta, indPhi)
    N_consumption = N_consumption * cell_volume(indR, indTheta, indPhi)
    E_rate = E_rate * cell_volume(indR, indTheta, indPhi)
    return
    end subroutine rate_split

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! creation rate for all spacial cells, done in parallel
! now split over by higher energy photon label
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rate_split_para(Edot)
    implicit none
    !output variable
    real*8, intent(out) :: Edot(Nomega) !erg/s of photon annihilation to pairs
    real*8, dimension(Nomega) :: Ehold, Eloss
    integer :: indR, indT, indP
    Edot = 0.d0;
    do indR = 1, Nr+1
        do indT = 1, 2*Ntheta
            !$OMP PARALLEL PRIVATE(indP, Eloss, Ehold)
            Eloss = 0.d0
            !$OMP DO
            do indP = 1, Nphi
                if (cell_volume(indR, indT, indP) > 0) then
                    call rate_split(indR, indT, indP, cell_create_split(indR, indT, indP, :),&
                         cell_loss_split(indR, indT, indP, :), Ehold)
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
    end subroutine rate_split_para

end module photon_split
