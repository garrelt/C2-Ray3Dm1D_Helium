module material_array

  use parameter, only: mesh
  use array, only: number_density_array, temperature_array, xHI_array, xHII_array, &
                   xHeI_array, xHeII_array, xHeIII_array

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine material_array_initialization()
    
    implicit none

      allocate(number_density_array(1:mesh,1:mesh,1:mesh))
      allocate(temperature_array(1:mesh,1:mesh,1:mesh))
      allocate(xHI_array(1:mesh,1:mesh,1:mesh))
      allocate(xHII_array(1:mesh,1:mesh,1:mesh))
      allocate(xHeI_array(1:mesh,1:mesh,1:mesh))
      allocate(xHeII_array(1:mesh,1:mesh,1:mesh))
      allocate(xHeIII_array(1:mesh,1:mesh,1:mesh))

  end subroutine material_array_initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module material_array
