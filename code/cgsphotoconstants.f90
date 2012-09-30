!>
!! \brief This module contains radiation-related physical constants
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date:
!!
!! \b Version: cgs units
!!

module cgsphotoconstants

  ! A collection of physical constants and conversion factors for 
  ! photo-ionization calculations
  ! Units: cgs
  
  use precision, only: dp
  use cgsconstants

  !	!> Helium ionization potentials (eV)
  !	real(kind=dp), dimension(0:1),parameter :: ethe=(/24.587,54.416/)
  !> Hydrogen cross section
  real(kind=dp), parameter :: sigh=6.346e-18
  !> Helium cross section
  real(kind=dp), parameter :: sighe0=7.43e-18
  !> He+ cross section
  real(kind=dp), parameter :: sighe1=1.589e-18
  !> H ionization energy in frequency
  real(kind=dp), parameter :: frth0=ev2fr*eth0
  !> He ionization energy in frequency
  real(kind=dp), parameter :: frthe0=ev2fr*ethe(0)
  !> He+ ionization energy in frequency
  real(kind=dp), parameter :: frthe1=ev2fr*ethe(1)
  !> Frequency dependence of H cross section parameter
  real(kind=dp),parameter :: betah0=1.0
  !> Frequency dependence of He cross section parameter
  real(kind=dp), dimension(0:1),parameter :: betahe=(/1.0,1.0/)

  !> Frequency dependence of H cross section parameter
  !real(kind=dp), parameter :: sh0=2.8
  !> Frequency dependence of He cross section parameter
  !real(kind=dp), parameter :: she0=1.9
  !> Frequency dependence of He+ cross section parameter
  !real(kind=dp), parameter :: she1=2.8

  !> maximum T_eff for black body
  real(kind=dp), parameter :: thigh=200000.0
  !> minimum T_eff for black body
  real(kind=dp), parameter :: tlow=2000.0
  !real(kind=dp), parameter :: frtop1=700.0*tlow/47979.72484*1e15

  !> Upper limits for integrals: due to arithmetic precision,
  !!     exp(700) exceeds double precision limit
  real(kind=dp), parameter :: frtop1=700.0*tlow*kb/hplanck
  !> Upper limits for integrals: due to the form of the planck
  !!     curve: take 10 times the frequency of maximum intensity
  real(kind=dp), parameter :: frtop2=5.88e-05*thigh*1e15
  
  !> H optical depth fit parameter for frequency range 2 & 3 (in 1 =1)
!  real(kind=dp) :: tfh(1:8) 
  
  !> He^0 optical depth fit parameter for frequency range 3 (in 2 =1) 
!  real(kind=dp) :: tfhe0(1:8) 

  !> He^0 optical depth fit parameter for frequency range 1,2 & 3,  =0 & =1 
!  real(kind=dp) :: tfhe1(1:8) 



  real(kind=dp) :: s_H_heth
  real(kind=dp) :: s_H_heLya
  real(kind=dp) :: s_He_heLya
  real(kind=dp) :: s_He_he2
  real(kind=dp) :: s_H_he2
  integer :: xix
  ! 
contains

  !> This subroutine initializes the factors needed
  !! to combine the three optical depths (H, He0, He+) 
  !! in the spectral regions where they all occur.
  !! This is from Tenorio-Tagle et al. (1983)
  !! It also contains the optical depths for doric
  subroutine ini_He_factors ( )

    !These values are taken from Verner et al. 1996
    !> opt depth of H at he0 ion threshold
    !s_H_heth  = sigh*(frthe0/frth0)**(-sh0)       
     s_H_heth  = 1.238e-18_dp     
   
    !> opt depth of He0 at he0 ion threshold
    !s_He_heth = sighe0 !not needed because it is only sighe0      
                             
    !> opt depth of H  at he+Lya  (this is 6.21e16 Hz
    !s_H_heLya = sigh*(40.817_dp*ev2fr/frth0)**(-sh0)      
     s_H_heLya = 9.907e-22_dp

    !> opt depth of He at he+Lya 
    !s_He_heLya= sighe0*(40.817_dp*ev2fr/frthe0)**(-she0)   
    s_He_heLya= 1.301e-20_dp
 
    ! He0 cross-section at He+ ionization threshold
    s_He_he2=1.690780687052975e-18_dp
    
    ! H0 cross-section at He+ ionization threshold 
    s_H_he2=1.230695924714239e-19_dp

  end subroutine ini_He_factors

end module cgsphotoconstants




