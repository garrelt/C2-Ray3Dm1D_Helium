!>
!! \brief This module handles the time variables
!!
!! Module for C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 20-Aug-2006
!<
module time

  ! This module handles the time variables

  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use astroconstants, only: YEAR
  use cosmology, only: t0, zred_t0, redshift_to_time
  use parameter, only: number_timesteps,number_outputs

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Sets time steps for calculation and outputting
  subroutine set_timesteps (begin_redshift,end_redshift,end_time,dt,output_dt)

    real(kind=dp),intent(in) :: begin_redshift !< starting redshift
    real(kind=dp),intent(in) :: end_redshift !< ending redshift
    real(kind=dp),intent(out) :: end_time !< end time
    real(kind=dp),intent(out) :: dt !< time step
    real(kind=dp),intent(out) :: output_dt !< output time step
    real(kind=dp) :: current_time

    ! Convert to time (in seconds)
    current_time = redshift_to_time(begin_redshift)
    end_time = redshift_to_time(end_redshift) 

    ! Set value of time step
    dt = (end_time-current_time)/real(number_timesteps)

    ! Convert to time 
    output_dt = (end_time-current_time)/real(number_outputs)

    if (rank .eq. 0) then
      write(logf,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(logf,*) "In the loop of redshift ",begin_redshift
      write(logf,*) "The evolution time-step = ",dt/YEAR,"years"
      write(logf,*) "The output time-step = ",output_dt/YEAR,"years"
      flush(logf) 
    endif

  end subroutine set_timesteps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module time
