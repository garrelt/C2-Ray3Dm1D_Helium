module xfrac

  use precision, only: dp
  use c2ray_parameters, only: epsilon_value
  use array, only: xHI_array,xHII_array
  use array, only: xHeI_array,xHeII_array,xHeIII_array

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine xfrac_initialization

    implicit none

    xHI_array(:,:,:) = 1.0-epsilon_value
    xHII_array(:,:,:) = epsilon_value
    xHeI_array(:,:,:) = 1.0-2.0_dp*epsilon_value
    xHeII_array(:,:,:) = epsilon_value	
    xHeIII_array(:,:,:) = epsilon_value 

    !xHI_array(:,:,:) = 1.0e-3
    !xHII_array(:,:,:) = 1.0-1-0e-3
    !xHeI_array(:,:,:) = 1.0e-3
    !xHeII_array(:,:,:) = 1.0-1.0e-3-epsilon_value	
    !xHeIII_array(:,:,:) = epsilon_value 

  end subroutine xfrac_initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module xfrac
