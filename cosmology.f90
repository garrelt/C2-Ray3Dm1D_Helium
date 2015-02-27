!>
!! \brief This module contains data and routines for cosmological problems
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 
!!
!! This module keeps track of the current redshift
!!

module cosmology

  ! This file contains routines having to do with the cosmological 
  ! evolution
  
  ! - cosmology_initialization: initializes cosmological time and sets lengths
  !            and volumes from comoving to proper scaling.
  ! - time_to_redshift: convert time to redshift
  ! - redshift_to_time: convert redshift to time
  ! - redshift_evolution: calculate redshift from time, and zfactor between the
  !    current and previous time
  ! - cosmological_evolution: cosmological evolution of space, density
  ! - cosmo_cool: cosmological adiabatic cooling rate
  ! - compton_cool: Compton cooling wrt the CMB.

  use precision, only: dp
  use my_mpi
  use cosmology_parameters, only: H0,Omega0,cmbtemp,cosmo_id
  use file_admin, only: logf
  use cgsconstants, only: c
  use c2ray_parameters, only: cosmological
  use array, only: number_density_array
  use redshift, only: initial_redshift
  use grid, only: x_coordinate, y_coordinate, z_coordinate
  use parameter, only: cell_size, cell_volume
  use parameter, only: current_redshift    !< current redshift

  implicit none

  real(kind=dp) :: zred_t0 !< initial redshift
  real(kind=dp) :: t0      !< time of initial redshift
  real(kind=dp) :: Hz      !< Hubble constant at current redshift
  real(kind=dp) :: relative_scale_factor !< scaling factor between two redshifts

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! initializes cosmological time and sets lengths and cell_volumes from 
  ! comoving to proper scaling.
  subroutine cosmology_initialization()
    
    implicit none

    real(kind=dp) :: relative_scale_factor3

    ! Cosmological time corresponding to (initial) redshift initial_redshift
    ! NOTE: Good only for high-z!!!
    t0 = 2.*(1.+initial_redshift)**(-1.5)/(3.*H0*sqrt(Omega0))
    
    current_redshift = initial_redshift ! needs to be zero, so comoving will be changed to proper

    ! Take the average relative_scale_factor between previous_redshift and current_redshift
    relative_scale_factor = 1.0/(1.+current_redshift) ! a = 1/(1+z)

    ! Hubble constant at current redshift (cgs)
    Hz = H0*(1.+current_redshift)**(1.5)*sqrt(Omega0) 

    relative_scale_factor3 = relative_scale_factor*relative_scale_factor*relative_scale_factor

    ! Change the grid coordinates
    x_coordinate(:) = x_coordinate(:)*relative_scale_factor
    y_coordinate(:) = y_coordinate(:)*relative_scale_factor
    z_coordinate(:) = z_coordinate(:)*relative_scale_factor
    cell_size = cell_size*relative_scale_factor
    cell_volume = cell_volume*relative_scale_factor3

  end subroutine cosmology_initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculates the cosmological redshift from time
  ! and the scale factor relative_scale_factor for use in cosmological_evolution
  subroutine redshift_evolution(time)

    implicit none

    real(kind=dp),intent(in) :: time
    real(kind=dp) :: previous_redshift

    ! Calculate the change since the previous redshift.
    ! Note: the initial redshift should be ZERO since
    ! the variables are initialized as comoving!
    ! NOTE: Good only for high-z!!!
    previous_redshift = current_redshift
    current_redshift = -1+(1.+initial_redshift)*((t0+time)/t0)**(-2./3.)

    ! Take the average relative_scale_factor between previous_redshift and current_redshift
    !relative_scale_factor = (1.0+previous_redshift)/(1.+current_redshift) ! a1/a2 = (1+z2)/(1+z1)

    relative_scale_factor = 1.0 ! Test 1

    ! Hubble constant at current redshift (cgs)
    Hz = H0*(1.+current_redshift)**(1.5)*sqrt(Omega0) 

    if (rank .eq. 0) then
      write(logf,*) "Redshift = ",current_redshift  
    endif

  end subroutine redshift_evolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculates the cosmological evolution of space and densities\n
  ! Changes variables: number_density_array, r, cell_size, cell_vol
  subroutine cosmological_evolution()

    ! Calculates the cosmological evolution of space and densities

    ! Author: Garrelt Mellema
    ! Date: 04-Mar-2006
    ! Version: F90 first version

    ! History:
    ! - 19-Nov-2004: first version f77

    !use sizes
    use grid, only: x_coordinate, y_coordinate, z_coordinate, cell_size, cell_volume
    
    real(kind=dp) :: relative_scale_factor3

    relative_scale_factor3 = relative_scale_factor*relative_scale_factor*relative_scale_factor

    ! Change the grid coordinates
    x_coordinate(:) = x_coordinate(:)*relative_scale_factor
    y_coordinate(:) = y_coordinate(:)*relative_scale_factor
    z_coordinate(:) = z_coordinate(:)*relative_scale_factor
    !cell_size = cell_size*relative_scale_factor  ! in test4, no need do this
    !cell_volume = cell_volume*relative_scale_factor3  ! in test4, no need do this   

    ! Change the densities
    !number_density_array(:,:,:) = number_density_array(:,:,:)/relative_scale_factor3  ! in test4, no need do this

  end subroutine cosmological_evolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculates the cosmological redshift for a given time
  real(kind=dp) function time_to_redshift(time)

    ! Calculates the cosmological redshift for a given time

    ! Author: Garrelt Mellema

    ! Date: 20-Aug-2006 (f77: 21-May-2005)
    
    ! Version: f90

    ! History: - 20-Aug-2006, conversion to f90

    real(kind=dp),intent(in) :: time

    ! Calculate the redshift
    ! NOTE: Good only for high-z!!!
    time_to_redshift = -1+(1.+initial_redshift)*(t0/(t0+time))**(2./3.)

  end function time_to_redshift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculates the time for a given cosmological redshift
  real(kind=dp) function redshift_to_time(zred1)

    real(kind=dp),intent(in) :: zred1

    ! Calculate the redshift
    ! NOTE: Good only for high-z!!!
    redshift_to_time = t0*( ((1.0+initial_redshift)/(1.0+zred1))**1.5 - 1.0 )

  end function redshift_to_time
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculates the cosmological adiabatic cooling
  real(kind=dp) function cosmo_cool (e_int)

    ! Calculates the cosmological adiabatic cooling

    ! Author: Garrelt Mellema
    ! Date: 04-Mar-2006
    ! Version: f90

    ! History:
    ! 19-Nov-2004: first version f77

    real(kind=dp),intent(in) :: e_int

    real(kind=dp) :: dzdt

    ! Cosmological cooling rate:
    ! 2*(da/dt)/a*e
    ! or
    ! 2*(dz/dt)/(1+z)*e
    ! with a the cosmological scale factor

    ! dz/dt (for flat LambdaCDM)
    dzdt=H0*(1.+current_redshift)*sqrt(Omega0*(1.+current_redshift)**3+1.-Omega0)

    !Cooling rate
    cosmo_cool=e_int*2.0/(1.0+current_redshift)*dzdt

  end function cosmo_cool

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculates the (cosmological) Compton cooling rate (against the CMB)
  real(kind=dp) function compton_cool (temper,eldens)
    
    ! Calculates the (cosmological) Compton cooling rate
    
    ! Author: Garrelt Mellema
    ! Date: 04-Mar-2006
    ! Version: first version
    
    ! History:
    ! 16-May-2005: first version f77
    

    ! parameter reference?

    real(kind=dp),intent(in) :: temper ! temperature
    real(kind=dp),intent(in) :: eldens ! electron density
    
    !Cooling rate
    compton_cool=5.65e-36*eldens*(1.0+current_redshift)**4* &
         (temper-cmbtemp*(1.0+current_redshift))
    
  end function compton_cool
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module cosmology
