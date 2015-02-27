module clumping

  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use parameter, only: current_redshift 
  use parameter, only: type_of_clumping, constant_clumping_factor

  implicit none

  real :: clumping_factor

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine clumping_evolution()
   
    implicit none

    select case (type_of_clumping)
    case(1)
       clumping_factor = constant_clumping_factor
    case(2) 
       clumping_factor = 27.466*exp(-0.114*current_redshift+0.001328*current_redshift*current_redshift)
    case(3)
       clumping_factor = 26.2917*exp(-0.1822*current_redshift+0.003505*current_redshift*current_redshift)
    case(4)
       clumping_factor = 17.57*exp(-0.101*current_redshift+0.0011*current_redshift*current_redshift)
    end select
   
    if (rank .eq. 0) then
      write(logf,*) "Clumping factor = ", clumping_factor
      flush(logf) 
    endif

  end subroutine clumping_evolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module clumping

