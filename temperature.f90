module temperature

  use precision, only: dp
  use my_mpi
  use array, only: temperature_array
  use parameter, only: input_IGM_temperature

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine temperature_initialization()
 
    implicit none
    
    temperature_array(:,:,:) = input_IGM_temperature

  end subroutine temperature_initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module temperature
