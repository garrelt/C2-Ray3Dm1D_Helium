! This module contains routines having to do with the calculation of
! the thermal evolution of a single point/cell. 

module thermalevolution

  use precision, only: dp
  use c2ray_parameters, only: minitemp, relative_denergy
  use cooling, only: coolin
  use tped, only: temper2pressr, pressr2temper, electrondens
  use cgsconstants
  use atomic
  use cosmology
  use type_definition, only: photrates,ionstates

  implicit none

contains

  ! calculates the thermal evolution of one grid point
  subroutine thermal (dt,end_temper,avg_temper,ndens_electron,ndens_atom,&
                      begin_xHI, begin_xHII, begin_xHeI, begin_xHeII, begin_xHeIII, &
                      avg_xHI, avg_xHII, avg_xHeI, avg_xHeII, avg_xHeIII, &
                      end_xHI, end_xHII, end_xHeI, end_xHeII, end_xHeIII, & 
                      heating_rate,pos)

    implicit none

    ! The time step
    real(kind=dp), intent(in) :: dt
    ! end time temperature of the cell
    real(kind=dp), intent(inout) :: end_temper
    ! average temperature of the cell
    real(kind=dp), intent(out) :: avg_temper
    ! Electron density of the cell
    real(kind=dp), intent(in) :: ndens_electron
    ! Number density of atoms of the cell
    real(kind=dp), intent(in) :: ndens_atom
    real(kind=dp),intent(in) :: begin_xHI
    real(kind=dp),intent(in) :: begin_xHII
    real(kind=dp),intent(in) :: begin_xHeI
    real(kind=dp),intent(in) :: begin_xHeII
    real(kind=dp),intent(in) :: begin_xHeIII
    real(kind=dp),intent(in) :: avg_xHI
    real(kind=dp),intent(in) :: avg_xHII
    real(kind=dp),intent(in) :: avg_xHeI
    real(kind=dp),intent(in) :: avg_xHeII
    real(kind=dp),intent(in) :: avg_xHeIII
    real(kind=dp),intent(in) :: end_xHI
    real(kind=dp),intent(in) :: end_xHII
    real(kind=dp),intent(in) :: end_xHeI
    real(kind=dp),intent(in) :: end_xHeII
    real(kind=dp),intent(in) :: end_xHeIII
    real(kind=dp),intent(in) :: heating_rate
 integer,dimension(1:3),intent(in) :: pos

    ! initial temperature
    real(kind=dp) :: initial_temp
    ! timestep taken to solve the ODE
    real(kind=dp) :: dt_ODE
    ! timestep related to thermal timescale
    real(kind=dp) :: dt_thermal
    ! record the time elapsed
    real(kind=dp) :: timestep
    ! internal energy of the cell
    real(kind=dp) :: internal_energy
    ! thermal timescale, used to calculate the thermal timestep
    real(kind=dp) :: thermal_timescale
    ! heating rate
    real(kind=dp) :: heating
    ! cooling rate
    real(kind=dp) :: cooling
    ! difference of heating and cooling rate
    real(kind=dp) :: thermal_rate
    ! cosmological cooling rate
    real(kind=dp) :: cosmo_cool_rate
    ! Counter of number of thermal timesteps taken
    integer :: i_heating

    ! heating rate
    heating = heating_rate

    ! Find initial internal energy
    internal_energy = temper2pressr(end_temper,ndens_atom, &
             electrondens(ndens_atom,begin_xHII,begin_xHeII,begin_xHeIII))/(gamma1)

    ! Set the cosmological cooling rate
    if (cosmological) then
       ! Disabled for testing
       cosmo_cool_rate=cosmo_cool(internal_energy)
    else
       cosmo_cool_rate=0.0
    endif

    ! Thermal process is only done if the temperature of the cell 
    ! is larger than the minimum temperature requirement

    if (end_temper.gt.minitemp) then

       ! stores the time elapsed is done
       timestep = 0.0 
   
       ! initialize the counter
       i_heating = 0

       ! initialize time averaged temperature
       avg_temper = 0.0 

       ! initial temperature
       initial_temp = end_temper  

       ! thermal process begins
       do

          ! update counter              
          i_heating = i_heating+1 
         
          ! update cooling rate from cooling tables
          cooling = coolin(ndens_atom,ndens_electron,avg_xHI,avg_xHII,avg_xHeI,avg_xHeII,avg_xHeIII,end_temper)+ &
                    cosmo_cool_rate
!cooling = 0.0

          ! Find total energy change rate
          thermal_rate = max(1d-50,abs(cooling-heating))

          ! Calculate thermal time scale
          thermal_timescale = internal_energy/abs(thermal_rate)

          ! Calculate time step needed to limit energy change
          ! to a fraction relative_denergy
          dt_thermal = relative_denergy*thermal_timescale

          ! Time step to large, change it to dt_thermal
          ! Make sure we do not integrate for longer than the
          ! total time step
          dt_ODE = min(dt_thermal,dt-timestep)

          ! Find new internal energy density
          internal_energy = internal_energy+dt_ODE*(heating-cooling)

          ! Update avg_temper sum (first part of dt_thermal sub time step)
          avg_temper = avg_temper+0.5*end_temper*dt_ODE

          ! Find new temperature from the internal energy density
          end_temper = pressr2temper(internal_energy*gamma1,ndens_atom, &
               electrondens(ndens_atom,avg_xHII,avg_xHeII,avg_xHeIII))

          ! Update avg_temper sum (second part of dt_thermal sub time step)
          avg_temper = avg_temper+0.5*end_temper*dt_ODE
                    
          ! Take measures if temperature drops below minitemp
          if (end_temper.lt.minitemp) then
             internal_energy = temper2pressr(minitemp,ndens_atom, &
                  electrondens(ndens_atom,avg_xHII,avg_xHeII,avg_xHeIII))
             end_temper = minitemp
          endif
                    
          ! Update fractional timestep
          timestep = timestep+dt_ODE
  
          ! Exit if we reach dt
          if (timestep.ge.dt.or.abs(timestep-dt).lt.1e-6*dt) exit
        	     	
       enddo
              
       ! Calculate the averaged temperature
       if (dt.gt.0.0) then

          avg_temper = avg_temper/dt

       else
          avg_temper = initial_temp
       endif
       
       ! Calculate the final temperature 
       end_temper = pressr2temper(internal_energy*gamma1,ndens_atom, &
            electrondens(ndens_atom,end_xHII,end_xHeII,end_xHeIII))

       
    endif
    
  end subroutine thermal
  
end module thermalevolution
