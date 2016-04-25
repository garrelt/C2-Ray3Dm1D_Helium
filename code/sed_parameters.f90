!>
!! \brief This module contains SED parameters
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! Note: this file contains parameters for spectral energy distributions
!!       (SEDs) for different types of SEDs:
!!       - Black body
!!       - General power law 
!!       - Quasar power law

module sed_parameters

  use precision, only: dp
  use cgsconstants, only: ev2fr
  use cgsphotoconstants, only: ion_freq_HeII
  use astroconstants, only: YEAR
  use sizes, only: mesh

  implicit none

  logical,parameter :: ask_for_sed = .false.

  !> Parameters for nominal SED (Black Body)
  !> Effective temperature (K)
  real(kind=dp),parameter :: T_eff_nominal=5.0e4
  !> Number of ionizing photons / second
  real(kind=dp),parameter :: bb_S_star_nominal=1e48_dp

  !> Parameters for nominal SED (General Power law)
#ifdef PL
  !> nominal Eddington efficiency
  real(kind=dp),parameter :: pl_EddLeff_nominal=1.0_dp
  !> nominal power law index (for photon number)
  real(kind=dp),parameter :: pl_index_nominal=2.5_dp
  !> nominal black hole mass for Eddington luminosity (M0)
  real(kind=dp),parameter :: pl_mass_nominal=1.0e6_dp
  !> Eddington luminosity per mass_nominal solar mass (erg/s)
  real(kind=dp),parameter :: pl_EddLum=1.38e38*pl_mass_nominal
  !> Number of ionizing photons / second
  real(kind=dp),parameter :: pl_S_star_nominal=1e48_dp
  !> nominal minimum and maximum frequency for power law source
  real(kind=dp),parameter :: pl_MinFreq_nominal=0.3*1e3*ev2fr
  real(kind=dp),parameter :: pl_MaxFreq_nominal=ion_freq_HeII * 100.00_dp
#endif

  !> Parameters for nominal SED (Quasar Power law)
#ifdef QUASARS
  !> nominal quasar Eddington efficiency
  real(kind=dp),parameter :: qpl_EddLeff_nominal=1.0_dp
  !> nominal quasar index (for photon number)
  real(kind=dp),parameter :: qpl_index_nominal=1.8_dp
  !> nominal quasar black hole mass for Eddington luminosity (M0)
  real(kind=dp),parameter :: qpl_mass_nominal=1.0e6_dp
  !> Eddington luminosity per qmass_nominal solar mass (erg/s)
  real(kind=dp),parameter :: qpl_EddLum=1.38e38*qpl_mass_nominal
  !> Number of ionizing photons / second for quasars
  real(kind=dp),parameter :: qpl_S_star_nominal=1e48_dp
  !> nominal minimum and maximum frequency for quasar source
  real(kind=dp),parameter :: qpl_MinFreq_nominal=0.3*1e3*ev2fr
  real(kind=dp),parameter :: qpl_MaxFreq_nominal=ion_freq_HeII * 100.00_dp
#endif

end module sed_parameters
