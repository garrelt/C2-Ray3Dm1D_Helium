!>
!! \brief This module contains data and routines for calculating the photon statistics
!!
!! Module for C2-Ray
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date:  26-Feb-2008
!!
!! \b Version: 3D

module photonstatistics
  
  ! Photon statistics
  ! photon_loss: this is kept here, but calculated in the evolve module.

  use precision, only: dp
  use my_mpi, only: rank
  use file_admin, only: logf
  use cgsconstants, only: acolh0, acolhe0, acolhe1,brech0, breche0,breche1
  use cgsphotoconstants, only: sigh  
  use sizes, only: mesh
  use grid, only: vol
  use material, only: ndens, clumping, clumping_point
  use tped, only: electrondens
  use sourceprops, only: NormFlux, NumSrc
  use radiation, only: S_star, NumFreqBnd
  use c2ray_parameters, only: type_of_clumping
  use abundances, only: abu_he

  implicit none

  !> true if checking photonstatistics
  logical,parameter :: do_photonstatistics=.true.
  !> Total number of recombinations
  real(kind=dp) :: totrec
  ! >Total number of recombination ionizations
  real(kind=dp) :: recomions
  !> Total number of collisional ionizations
  real(kind=dp) :: totcollisions
  !> Change in number of neutral H atoms
  real(kind=dp) :: dh0
  !> Change in number of neutral He atoms
  real(kind=dp) :: dhe0
  !> Change in number of double ionized He+ ions
  real(kind=dp) :: dhe2
  !> Total number of ionizing photons used
  real(kind=dp) :: total_ion
  !> Total number of photons lost in LLSs
  real(kind=dp) :: LLS_loss  
  !> Grand total number of ionizing photons used
  real(kind=dp) :: grtotal_ion
  !> Grand total number of ionizing photons produced
  real(kind=dp) :: grtotal_src
  !> Number of photons leaving the grid in each frequency band
  real(kind=dp) :: photon_loss(NumFreqBnd)

  real(kind=dp),private :: h0_before !< number of H atoms at start of time step
  real(kind=dp),private :: h0_after !< number of H atoms at end of time step
  real(kind=dp),private :: h1_before !< number of H ions at start of time step
  real(kind=dp),private :: h1_after !< number of H ions at end of time step
  real(kind=dp),private :: he0_before !< number of He atoms at start of time step
  real(kind=dp),private :: he0_after !< number of He atoms at end of time step
  real(kind=dp),private :: he1_before !< number of He+ ions at start of time step
  real(kind=dp),private :: he1_after !< number of He+ ions at end of time step
  real(kind=dp),private :: he2_before !< number of He++ ions at start of time step
  real(kind=dp),private :: he2_after !< number of He++ ions at end of time step


  integer,private :: i !< mesh loop index (x)
  integer,private :: j !< mesh loop index (y)
  integer,private :: k !< mesh loop index (z)

contains

  !----------------------------------------------------------------------------

  !> Initialize the photon statistics
  subroutine initialize_photonstatistics ()

    ! set total number of ionizing photons used to zero
    grtotal_ion=0.0
    grtotal_src=0.0

  end subroutine initialize_photonstatistics

  !----------------------------------------------------------------------------

  !> Call the individual routines needed for photon statistics calculation
  subroutine calculate_photon_statistics (dt,xh_l,xh_r,xhe_l,xhe_r)

    real(kind=dp),intent(in) :: dt !< time step
    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:1),intent(in) :: xh_l
    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:2),intent(in) :: xhe_l
    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:1),intent(in) :: xh_r
    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:2),intent(in) :: xhe_r
    ! Call the individual routines needed for this calculation

    call state_after (xh_l,xhe_l) ! number of neutrals after integration
    call total_rates (dt,xh_r,xhe_r) ! total photons used in balancing recombinations etc.
    call total_ionizations () ! final statistics
    
  end subroutine calculate_photon_statistics

  !----------------------------------------------------------------------------

  !> Calculates the number of neutrals and ions at the start of the time step
  subroutine state_before (xh_l,xhe_l)

    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:1),intent(in) :: xh_l
    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:2),intent(in) :: xhe_l
    ! Photon statistics: calculate the number of neutrals before integration
    h0_before=0.0
    h1_before=0.0
    he0_before=0.0
    he1_before=0.0
    he2_before=0.0
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             h0_before=h0_before+ndens(i,j,k)*xh_l(i,j,k,0)
             h1_before=h1_before+ndens(i,j,k)*xh_l(i,j,k,1)
             he0_before=he0_before+ndens(i,j,k)*xhe_l(i,j,k,0)
             he1_before=he1_before+ndens(i,j,k)*xhe_l(i,j,k,1)
             he2_before=he2_before+ndens(i,j,k)*xhe_l(i,j,k,2)
          enddo
       enddo
    enddo

    h0_before=h0_before*vol*(1.0_dp-abu_he)
    h1_before=h1_before*vol*(1.0_dp-abu_he)
    he0_before=he0_before*vol*abu_he
    he1_before=he1_before*vol*abu_he
    he2_before=he2_before*vol*abu_he
  end subroutine state_before

  !----------------------------------------------------------------------------

  !> Calculates total number of recombinations and collisions
!*** this is generally colled with the average values
  subroutine total_rates(dt,xh_l, xhe_l)

    real(kind=dp),intent(in) :: dt !< time step
    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:1),intent(in) :: xh_l
    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:2),intent(in) :: xhe_l

    real(kind=dp),dimension(0:1) :: yh
    real(kind=dp),dimension(0:2) :: yhe
    real(kind=dp) :: ndens_p ! needed because ndens may be single precision

    ! Photon statistics: Determine total number of recombinations/collisions
    ! Should match the code in doric_module

    totrec=0.0
    totcollisions=0.0
    recomions=0.0
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             yh(0)=xh_l(i,j,k,0)
             yh(1)=xh_l(i,j,k,1)
             yhe(0)=xhe_l(i,j,k,0)
             yhe(1)=xhe_l(i,j,k,1)
             yhe(2)=xhe_l(i,j,k,2)
             ndens_p=ndens(i,j,k)
             ! Set clumping to local value if we have a clumping grid
             if (type_of_clumping == 5) &
                  call clumping_point (i,j,k)
             
             totrec=totrec+ndens_p* &
                  (xh_l(i,j,k,1)*brech0*(1.0_dp-abu_he) +   &
                  xhe_l(i,j,k,1)*breche0*abu_he*0.04_dp    &
                  !   xhe_l(i,j,k,1)*(breche0*abu_he)         &
                  )*electrondens(ndens_p,yh,yhe)*clumping
             
             totcollisions=totcollisions+ &
                  ndens_p*electrondens(ndens_p,yh,yhe)* ( &
                  xh_l(i,j,k,0)*acolh0   +                            &
                  xhe_l(i,j,k,0)*acolhe0 +                            & 
                  xhe_l(i,j,k,1)*acolhe1)

             recomions= recomions+ &
                  ndens_p*abu_he*clumping*(xhe_l(i,j,k,2)*1.121_dp*breche1+  &   ! This is: 1+vfrac*1.425-vfrac)
                  xhe_l(i,j,k,1)*breche0*0.96_dp)* &
                  abu_he*electrondens(ndens_p,yh,yhe)
          enddo
       enddo
    enddo
   
    totrec=totrec*vol*dt
    totcollisions=totcollisions*vol*dt
    recomions=recomions*vol*dt

  end subroutine total_rates
  
  !----------------------------------------------------------------------------

  !> Calculates the number of neutrals and ions at the end of the time step
  subroutine state_after(xh_l,xhe_l)
    
    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:1),intent(in) :: xh_l
    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),0:2),intent(in) :: xhe_l
    ! Photon statistics: Calculate the number of neutrals after the integration
    h0_after=0.0
    h1_after=0.0
    he0_after=0.0
    he1_after=0.0
    he2_after=0.0
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             h0_after=h0_after+ndens(i,j,k)*xh_l(i,j,k,0)
             h1_after=h1_after+ndens(i,j,k)*xh_l(i,j,k,1)
             he0_after=he0_after+ndens(i,j,k)*xhe_l(i,j,k,0)
             he1_after=he1_after+ndens(i,j,k)*xhe_l(i,j,k,1)
             he2_after=he2_after+ndens(i,j,k)*xhe_l(i,j,k,2)
          enddo
       enddo
    enddo
    h0_after=h0_after*vol*(1.0_dp-abu_he)
    h1_after=h1_after*vol*(1.0_dp-abu_he)
    he0_after=he0_after*vol*abu_he   
    he1_after=he1_after*vol*abu_he   
    he2_after=he2_after*vol*abu_he   
  end subroutine state_after
  
  !----------------------------------------------------------------------------

  !> Calculate the total number of ionizing photons used
  subroutine total_ionizations ()
    
    ! Photon statistics: Total number of new ionizations
    dh0=(h0_before-h0_after)    ! count the old neutrals - new neutrals --> ionizations
    dhe0=(he0_before-he0_after) ! count the old neutrals - new neutrals --> ionizations
    dhe2=(he2_after-he2_before) ! count the new ions     - old ions     --> ionizations 
    total_ion=dh0+dhe0+dhe2  !add them up
    
  end subroutine total_ionizations
  
  !----------------------------------------------------------------------------
  !> Calculate the total number of photons lost in LLSs
  subroutine total_LLS_loss (phi_out,coldensh_LLS)

    real(kind=dp),intent(in) :: phi_out !< Photon flux out off cell
    real(kind=dp),intent(in) :: coldensh_LLS !< local (cell) column density of LLSs

    ! Note LLS_loss is set to zero at the start of a time step in
    ! evolve:pass_all_sources
    real(kind=dp) :: tau_LLS

    tau_LLS=sigh*coldensh_LLS
    ! Note: This expression is only valid for grey opacities!!
    ! GM (110131): Why do we divide by exp(-tau_LLS)? The numbers
    ! generated by this expression seem to be wrong. Comment out
    ! exp(-tau_LLS)
    LLS_loss = LLS_loss + phi_out*(1.0-exp(-tau_LLS))!/exp(-tau_LLS)

  end subroutine total_LLS_loss

  !----------------------------------------------------------------------------

  !> Calculate the total number of ionizing photons used
  subroutine report_photonstatistics (dt)

    real(kind=dp),intent(in) :: dt

    real(kind=dp) :: totalsrc,photcons,total_photon_loss,total_LLS_loss

    total_photon_loss=sum(photon_loss)*dt* &
         real(mesh(1))*real(mesh(2))*real(mesh(3))
    total_LLS_loss = LLS_loss*dt         
    totalsrc=sum(NormFlux(1:NumSrc))*s_star*dt
    photcons=(total_ion-totcollisions-recomions)/(totalsrc)
    if (rank == 0) then
       write(90,"(9(1pe10.3))") &
            total_ion, totalsrc, &          ! # of new ionizations (after-before) | # of emitted photons
            recomions, total_photon_loss, &                  ! # number of ionizations due to He recombination ! # photon loss  
            totrec, totcollisions, &                       ! recombinations that don't ionize;  
            totrec/total_ion, &             !  fraction of recombinations 
            total_photon_loss/totalsrc, &   ! lost photon fraction
            totcollisions/total_ion         ! fraction of collisions
       write(logf,*) h1_before,h1_after
    endif
    
  end subroutine report_photonstatistics

  !----------------------------------------------------------------------------

  !> Calculate the total number of ionizing photons produced
  subroutine update_grandtotal_photonstatistics (dt)

    real(kind=dp),intent(in) :: dt !< time step

    grtotal_src=grtotal_src+sum(NormFlux(1:NumSrc))*s_star*dt
    grtotal_ion=grtotal_ion+total_ion-totcollisions

  end subroutine update_grandtotal_photonstatistics

end module photonstatistics
