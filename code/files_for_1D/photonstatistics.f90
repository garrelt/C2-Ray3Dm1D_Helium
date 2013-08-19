!>
!! \brief This module contains data and routines for calculating the photon statistics
!!
!! Module for C2-Ray
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 05-Mar-2010
!!
!! \b Version: 1D

module photonstatistics
  
  ! Photon statistics
  ! photon_loss is a sum over all sources, summing is done in evolve0d.
  ! For parallelization consider a reduction, or may also be dropped
  ! entirely.
  ! Photon_losst contains the threaded values

  use precision, only: dp
  use my_mpi, only: rank
  use file_admin, only: logf
  use cgsconstants, only: albpow,bh00,colh0,temph0
  use cgsconstants, only: alcpow,bhe00,bhe10,colhe,temphe
  use sizes, only: mesh
  use grid, only: vol
  use abundances, only: abu_he
  use material, only: ndens, xh, xhe, temper, clumping
  use tped, only: electrondens
  use radiation, only: S_star, NumFreqBnd

  implicit none

  !> true if checking photonstatistics
  logical,parameter :: do_photonstatistics=.true.
  !> Total number of recombinations
  real(kind=dp) :: totrec
  !> Total number of collisional ionizations
  real(kind=dp) :: totcollisions
  !> Change in number of neutral H, He atoms
  real(kind=dp) :: dh0,dhe0,dhe1,dhe2, dh1
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
  real(kind=dp),private :: he0_before,he0_after,he1_before,he1_after
  real(kind=dp),private :: he2_before,he2_after
  integer,private :: i,j,k !< mesh loop index

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

    real(kind=dp),intent(in) :: dt !< time step

    ! Call the individual routines needed for this calculation

    call state_after () ! number of neutrals after integration
    call total_rates (dt) ! total photons used in balancing recombinations etc.
    call total_ionizations () ! final statistics
    
  end subroutine calculate_photon_statistics

  !----------------------------------------------------------------------------

  !> Calculates the number of neutrals and ions at the start of the time step
  subroutine state_before ()

    ! Photon statistics: calculate the number of neutrals before integration
    h0_before=0.0
    h1_before=0.0
    he0_before=0.0
    he1_before=0.0
    he2_before=0.0
    do i=1,mesh(1)
       h0_before=h0_before+vol(i)*ndens(i,1,1)*xh(i,0)*(1.0_dp-abu_he)
       h1_before=h1_before+vol(i)*ndens(i,1,1)*xh(i,1)*(1.0_dp-abu_he)
       he0_before=he0_before+vol(i)*ndens(i,1,1)*xhe(i,0)*abu_he
       he1_before=he1_before+vol(i)*ndens(i,1,1)*xhe(i,1)*abu_he
       he2_before=he2_before+vol(i)*ndens(i,1,1)*xhe(i,2)*abu_he

    enddo
    
  end subroutine state_before

  !----------------------------------------------------------------------------

  !> Calculates total number of recombinations and collisions
  subroutine total_rates(dt)

    real(kind=dp),intent(in) :: dt !< time step


    real(kind=dp),dimension(0:1) :: yh
    real(kind=dp),dimension(0:2) :: yhe
    real(kind=dp) :: ndens_p ! needed because ndens may be single precision

    ! Photon statistics: Determine total number of recombinations/collisions
    ! Should match the code in doric_module

    totrec=0.0
    totcollisions=0.0
    do i=1,mesh(1)
       yh(0)=xh(i,0)
       yh(1)=xh(i,1)
       yhe(0)=xhe(i,0)
       yhe(1)=xhe(i,1)
       yhe(2)=xhe(i,2)
       ndens_p=ndens(i,1,1)
       totrec=totrec+vol(i)*ndens_p* electrondens(ndens_p,yh,yhe)*clumping* &
            (xh(i,1)*(1.0_dp-abu_he)*  & 
            1.0_dp/(1.0_dp/(bh00*(temper(i)/1e4)**albpow)+1.0_dp/(bh00*5.0_dp*(temper(i)/1e4)**(1.95_dp*albpow))) +&
            xhe(i,1)*abu_he*  &
            (bhe00*(temper(i)/1e4)**alcpow) + &
            xhe(i,2)*abu_he*  &
            1.0_dp/(1.0_dp/(bhe10*(temper(i)/1e4)**(0.95_dp*albpow))+1.0_dp/(bhe10*11.0_dp*(temper(i)/1e4)**(albpow*1.95_dp))))


       totcollisions=totcollisions+vol(i)*ndens_p*(1.0_dp-abu_he)* &
            xh(i,0)*electrondens(ndens_p,yh,yhe)* &
            colh0*sqrt(temper(i))*exp(-temph0/temper(i))+ &
            vol(i)*ndens_p*abu_he*   &
            xhe(i,0)*electrondens(ndens_p,yh,yhe)* &
            colhe(0)*sqrt(temper(i))*exp(-temphe(0)/temper(i))+ &
            vol(i)*ndens_p*abu_he*   &
            xhe(i,1)*electrondens(ndens_p,yh,yhe)* &
            colhe(1)*sqrt(temper(i))*exp(-temphe(1)/temper(i))
    enddo

    totrec=totrec*dt
    totcollisions=totcollisions*dt

  end subroutine total_rates
  
  !----------------------------------------------------------------------------

  !> Calculates the number of neutrals and ions at the end of the time step
  subroutine state_after()

    ! Photon statistics: Calculate the number of neutrals after the integration
    h0_after=0.0
    h1_after=0.0
    he0_after=0.0
    he1_after=0.0
    he2_after=0.0

    do i=1,mesh(1)
       h0_after=h0_after+vol(i)*ndens(i,1,1)*xh(i,0)*(1.0_dp-abu_he)
       h1_after=h1_after+vol(i)*ndens(i,1,1)*xh(i,1)*(1.0_dp-abu_he)
       he0_after=he0_after+vol(i)*ndens(i,1,1)*xhe(i,0)*abu_he
       he1_after=he1_after+vol(i)*ndens(i,1,1)*xhe(i,1)*abu_he
       he2_after=he2_after+vol(i)*ndens(i,1,1)*xhe(i,2)*abu_he
    enddo
    
  end subroutine state_after
  
  !----------------------------------------------------------------------------

  !> Calculate the total number of ionizing photons used
  subroutine total_ionizations ()

    ! Photon statistics: Total number of new ionizations
    dh0=(h0_before-h0_after)
    dh1=(h1_before-h1_after)
    dhe0=(he0_before-he0_after)
    dhe1=(he1_before-he1_after)
    dhe2=(he2_before-he2_after)
    total_ion=totrec+dh0+dhe0+dhe1
    
  end subroutine total_ionizations

  !----------------------------------------------------------------------------

  !> Calculate the total number of ionizing photons used
  subroutine report_photonstatistics (dt)

    real(kind=dp),intent(in) :: dt

    real(kind=dp) :: totalsrc,photcons,total_photon_loss

    !total_photon_loss=sum(photon_loss)*dt* &
    total_photon_loss=photon_loss*dt* &
         real(mesh(1))
    totalsrc=s_star*dt
    photcons=(total_ion-totcollisions)/totalsrc
    if (rank == 0) then
       write(logf,"(7(1pe10.3))") &
            total_ion, totalsrc, &
            photcons, &
            dh0/total_ion, &
            totrec/total_ion, &
            total_photon_loss/totalsrc, &
            totcollisions/total_ion
    endif
    
  end subroutine report_photonstatistics

end module photonstatistics
