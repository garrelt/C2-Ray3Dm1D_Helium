!>
!! \brief This module contains routines having to do with the calculation of
!! the thermal evolution of a single point/cell. 
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 
!!

module thermalevolution

  ! This file contains routines having to do with the calculation of
  ! the thermal evolution of a single point/cell.
  ! It can used for Yguazu-a or non-hydro photo-ionization calculations.

  ! - thermal : time dependent solution of the radiative heating/cooling

  use precision, only: dp
  use c2ray_parameters, only: minitemp,relative_denergy
  use radiative_cooling, only: coolin
  use tped, only: temper2pressr,pressr2temper,electrondens
  use cgsconstants
  use atomic
  use cosmology
  use radiation, only: photrates
  use material, only: ionstates

  implicit none

contains

  !=======================================================================


  subroutine thermal (dt,temper,avg_temper,rhe,rhh,ion,phi)


    ! calculates the thermal evolution

    ! Author: Garrelt Mellema
    ! Date: 11-May-2005 (30 July 2004, 1 June 2004)

    ! Version:
    ! Simplified version for testing photon conservation.
    ! - sub-timesteps
    ! - H only, one frequency range
    ! - Cooling curve

    ! Changes
    ! 30 Jul 2004: the initial value for the internal energy is
    !              now calculated with the initial ionization
    !              fractions (xh0). The final with the final (xh), 
    !              and the intermediate temperatures with the average
    !              value (xh_av).
    ! 29 Jul 2004: call to electrondens was wrong. This function needs
    !              density and ionization fraction.

    !real(kind=dp),parameter :: minitemp=1.0_dp ! minimum temperature
    ! fraction of the cooling time step below which no iteration is done
    !real(kind=dp),parameter :: relative_denergy=0.1_dp    

    real(kind=dp),intent(in) :: dt
    real(kind=dp),intent(inout) :: temper
    real(kind=dp),intent(out) :: avg_temper
    real(kind=dp),intent(in) :: rhe,rhh !*,xh(0:1),xh_av(0:1),xh0(0:1)
!*    real(kind=dp),intent(in) :: xhe(0:2),xhe_av(0:2),xhe0(0:2)

    type(photrates),intent(in) :: phi
!*TEST
    type(ionstates), intent(in) :: ion
!*****

    real(kind=dp) :: temper0
    real(kind=dp) :: dt0,dt1,timestep,e_int,dt_thermal,eplus
    real(kind=dp) :: emin,eplush,eplushe0,eplushe1,thermalrt,cosmo_cool_rate
    integer :: nstep,nstepmax

    ! Photo-ionization heating
    ! GM/121115: This is the total heating rate including H, He0, He1 and
    ! secondary heating. See radiation module,  subroutine  lookuptable
    ! where the total heating variable fheat is defined and then put into
    ! phi%hv_h. We should probably change this, or the variable names.
    eplush=phi%hv_h

    !eplushe0=phi%hv_he(0)
    !eplushe1=phi%hv_he(1)

    eplus=eplush ! +eplushe0+eplushe1  (3.4e-13)*(temp0/1e4).^-0.60

    ! Find initial internal energy
!Change to the old here!
    e_int=temper2pressr(temper,rhh,electrondens(rhh,ion%h_old, ion%he_old))/ &
         (gamma1)

     ! write(*,*) e_int,eplus,phi%h

!    if (cosmological) then
       ! Disabled for testing
!       cosmo_cool_rate=cosmo_cool(e_int)
!    else
       cosmo_cool_rate=0.0
!    endif

    ! Do nothing if temperature is below minitemp

    if (temper.gt.minitemp) then
       !
       ! Do the cooling/heating.
       ! First figure out the cooling/heating rate (thermalrt) 
       ! and corresponding time step (dt_thermal). Then take
       ! a fraction relative_denergy of this time step if the real time step
       ! is larger than this. Loop through this until we reach
       ! the full time step (dt).
       ! We are effectively following the cooling curve here.
       ! Along the way we collect the temperatures passed so
       ! that an average temperature (avg_temper) can be calculated.
       !
       dt0=dt          
       dt1=dt
       timestep=0.0    ! stores how much of time step is done
       nstep=0         ! counter
       avg_temper=0.0  ! initialize time averaged temperature
       temper0=temper  ! initial temperature
       do
          ! update counter              
          nstep=nstep+1 
         
          ! Find cooling rate (using average ionization fraction)
          !emin=min(1d-50,1e-8*eplus)!

          emin=coolin(rhh,rhe,ion%h_av, ion%he_av, temper)

          emin=emin+cosmo_cool_rate

          ! Find total energy change rate
          thermalrt=max(1d-50,abs(emin-eplus))

          ! Calculate thermal time
          dt_thermal=e_int/abs(thermalrt)
          ! Calculate time step needed to limit energy change
          ! to a fraction relative_denergy
          dt1=relative_denergy*dt_thermal
          ! Time step to large, change it to dt1
          ! Make sure we do not integrate for longer than the
          ! total time step dtcgs
          dt0=min(dt1,dt-timestep)
          ! Find new internal energy density
          e_int=e_int+dt0*(eplus-emin)
          ! Update avg_temper sum (first part of dt1 sub time step)
          avg_temper=avg_temper+0.5*temper*dt0
          ! Find new temperature from the internal energy density
          temper=pressr2temper(e_int*gamma1,rhh,electrondens(rhh,ion%h_av, ion%he_av))
          ! Update avg_temper sum (second part of dt1 sub time step)
          avg_temper=avg_temper+0.5*temper*dt0
          
          
          ! Take measures if temperature drops below minitemp
          if (temper.lt.minitemp) then
             e_int=temper2pressr(minitemp,rhh,electrondens(rhh,ion%h_av,ion%he_av))
             temper=minitemp
          endif
                    
          ! Update fractional timestep
          timestep=timestep+dt0
  
          ! Exit if we reach dt
          ! Mind fp precision here, so check for nearness
          if (timestep.ge.dt.or.abs(timestep-dt).lt.1e-6*dt) exit
          	if (nstep.gt.10000)  then
          	write(nstep) 
          	exit
          	endif
          	
     	
       enddo
       

       
       ! Calculate time averaged temperature
       if (dt.gt.0.0) then
          avg_temper=avg_temper/dt
       else
          avg_temper=temper0
       endif
       
       !write(*,*) '--'
       !write(*,*) e_int,rhh,gamma1
       !pause
       ! Calculate temperature with final ionization fractions
       temper=pressr2temper(e_int*gamma1,rhh,electrondens(rhh,ion%h,ion%he))
       
    endif
    
  end subroutine thermal
  
end module thermalevolution
