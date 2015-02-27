module evolve_array

  use precision, only: dp
  use my_mpi
  use parameter, only: mesh
  use parameter, only: NumFreqBnd
  use array, only: intermediate_average_xHI,intermediate_average_xHII
  use array, only: intermediate_end_xHI,intermediate_end_xHII
  use array, only: intermediate_average_xHeI,intermediate_average_xHeII,intermediate_average_xHeIII
  use array, only: intermediate_end_xHeI,intermediate_end_xHeII,intermediate_end_xHeIII
  use array, only: intermediate_end_temperature,intermediate_average_temperature
  use array, only: photoionization_HI_array
  use array, only: photoionization_HeI_array
  use array, only: photoionization_HeII_array
  use array, only: heating_array
  use array, only: coldens_in_HI
  use array, only: coldens_out_HI
  use array, only: coldens_in_HeI
  use array, only: coldens_out_HeI
  use array, only: coldens_in_HeII
  use array, only: coldens_out_HeII
  use array, only: ring_vol
  use array, only: buffer

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine evolve_array_initialization()
   
    implicit none
 
    allocate(photoionization_HI_array(1:mesh,1:mesh,1:mesh))
    photoionization_HI_array = 0.0 
    allocate(photoionization_HeI_array(1:mesh,1:mesh,1:mesh))
    photoionization_HeI_array = 0.0 
    allocate(photoionization_HeII_array(1:mesh,1:mesh,1:mesh))
    photoionization_HeII_array = 0.0 
    allocate(heating_array(1:mesh,1:mesh,1:mesh))
    heating_array=0.0
    
    allocate(intermediate_average_xHI(1:mesh,1:mesh,1:mesh))
    allocate(intermediate_average_xHII(1:mesh,1:mesh,1:mesh))

    allocate(intermediate_average_xHeI(1:mesh,1:mesh,1:mesh))
    allocate(intermediate_average_xHeII(1:mesh,1:mesh,1:mesh))
    allocate(intermediate_average_xHeIII(1:mesh,1:mesh,1:mesh))

    allocate(intermediate_end_xHI(1:mesh,1:mesh,1:mesh))
    allocate(intermediate_end_xHII(1:mesh,1:mesh,1:mesh))

    allocate(intermediate_end_xHeI(1:mesh,1:mesh,1:mesh))
    allocate(intermediate_end_xHeII(1:mesh,1:mesh,1:mesh))
    allocate(intermediate_end_xHeIII(1:mesh,1:mesh,1:mesh))

    allocate(intermediate_end_temperature(1:mesh,1:mesh,1:mesh))
    allocate(intermediate_average_temperature(1:mesh,1:mesh,1:mesh))

    allocate(coldens_in_HI(1:mesh,1:mesh,1:mesh))
    allocate(coldens_out_HI(1:mesh,1:mesh,1:mesh))
    allocate(coldens_in_HeI(1:mesh,1:mesh,1:mesh))
    allocate(coldens_out_HeI(1:mesh,1:mesh,1:mesh))
    allocate(coldens_in_HeII(1:mesh,1:mesh,1:mesh))
    allocate(coldens_out_HeII(1:mesh,1:mesh,1:mesh))

    allocate(ring_vol(1:mesh,1:mesh,1:mesh))
    allocate(buffer(1:mesh,1:mesh,1:mesh))

  end subroutine evolve_array_initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module evolve_array
