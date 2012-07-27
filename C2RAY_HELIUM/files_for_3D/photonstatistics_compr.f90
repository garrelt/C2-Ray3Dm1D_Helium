!>
!! \brief This module contains data and routines for calculating the photon statistics
!!
!! Module for C2-Ray
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date:  26-Feb-2008
!!
!! \b Version: 

module photonstatistics
  
  ! This module handles the calculation of the photon statistics
  ! For: C2-Ray

  ! Author: Garrelt Mellema

  ! Date: 26-Feb-2008

  ! Photon statistics
  ! photon_loss is a sum over all sources, summing is done in evolve0d
  ! and evolve3d (in case of parallelization).
 
  use precision, only: dp
  use cgsconstants, only: albpow,bh00,colh0,temph0
  use sizes, only: mesh
  use grid, only: vol
  use material, only: ndens, xh_compr, temper, clumping, neutral_from_compr, ionized_from_compr
  use tped, only: electrondens

  !> true if checking photonstatistics
  logical,parameter :: do_photonstatistics=.true. 
  !> Total number of recombinations
  real(kind=dp) :: totrec
  !> Total number of collisional ionizations
  real(kind=dp) :: totcollisions
  !> Change in number of neutral H atoms
  real(kind=dp) :: dh0
  !> Total number of ionizing photons used
  real(kind=dp) :: total_ion
  !> Grand total number of ionizing photons used
  real(kind=dp) :: grtotal_ion
  !> Number of photons leaving the grid
  real(kind=dp) :: photon_loss

  real(kind=dp),private :: h0_before !< number of H atoms at start of time step
  real(kind=dp),private :: h0_after !< number of H atoms at end of time step
  real(kind=dp),private :: h1_before !< number of H ions at start of time step
  real(kind=dp),private :: h1_after !< number of H ions at end of time step

  integer,private :: i !< mesh loop index (x)
  integer,private :: j !< mesh loop index (y)
  integer,private :: k !< mesh loop index (z)

contains
  
  !----------------------------------------------------------------------------

  !> Initialize the photon statistics
  subroutine initialize_photonstatistics ()

    ! set total number of ionizing photons used to zero
    grtotal_ion=0.0

  end subroutine initialize_photonstatistics

  !----------------------------------------------------------------------------

  !> Call the individual routines needed for photon statistics calculation
  subroutine calculate_photon_statistics (dt)

    real(kind=dp),intent(in) :: dt

    ! Call the individual routines needed for this calculation

    call state_after () ! number of neutrals after integration
    call total_rates (dt) ! total photons used in balancing recombinations etc.
    call total_ionizations () ! final statistics
    
  end subroutine calculate_photon_statistics

  !----------------------------------------------------------------------------

  !> Calculate the state at the start of the time step
  subroutine state_before ()

    ! Photon statistics: calculate the number of neutrals before integration
    h0_before=0.0
    h1_before=0.0
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             h0_before=h0_before+ndens(i,j,k)*neutral_from_compr(xh_compr(i,j,k))
             h1_before=h1_before+ndens(i,j,k)*ionized_from_compr(xh_compr(i,j,k))
          enddo
       enddo
    enddo

    h0_before=h0_before*vol
    h1_before=h1_before*vol

  end subroutine state_before

  !----------------------------------------------------------------------------

  !> Calculate (sum) all the rates
  subroutine total_rates(dt)

    real(kind=dp),intent(in) :: dt

    real(kind=dp),dimension(0:1) :: yh
    real(kind=dp) :: ndens_p

    ! Photon statistics: Determine total number of recombinations/collisions
    ! Should match the code in doric_module

    totrec=0.0
    totcollisions=0.0
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             yh(0)=neutral_from_compr(xh_compr(i,j,k))
             yh(1)=ionized_from_compr(xh_compr(i,j,k))
             ndens_p=ndens(i,j,k)
             totrec=totrec+ndens_p*yh(1)*    &
                  electrondens(ndens_p,yh)*  &
                  clumping*bh00*(temper/1e4)**albpow
             totcollisions=totcollisions+ndens_p*   &
                  yh(0)*electrondens(ndens_p,yh)* &
                  colh0*sqrt(temper)*exp(-temph0/temper)
          enddo
       enddo
    enddo

    totrec=totrec*vol*dt
    totcollisions=totcollisions*vol*dt

  end subroutine total_rates
  
  !----------------------------------------------------------------------------

  !> Calculate the state at the end of the time step
  subroutine state_after()
    
    ! Photon statistics: Calculate the number of neutrals after the integration
    h0_after=0.0
    h1_after=0.0
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             h0_after=h0_after+ndens(i,j,k)*neutral_from_compr(xh_compr(i,j,k))
             h1_after=h1_after+ndens(i,j,k)*ionized_from_compr(xh_compr(i,j,k))
          enddo
       enddo
    enddo
    h0_after=h0_after*vol
    h1_after=h1_after*vol
    
  end subroutine state_after
  
  !----------------------------------------------------------------------------

  !> Calculate the total number of ionizing photons used
  subroutine total_ionizations ()
    
    ! Photon statistics: Total number of new ionizations
    dh0=(h0_before-h0_after)
    total_ion=totrec+dh0
    
  end subroutine total_ionizations

end module photonstatistics
