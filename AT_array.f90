module AT_array

  use precision, only: dp
  use parameter, only: mesh, time_partition, number_of_source, do_periodic, &
                       length_pos_x, length_neg_x, length_pos_y, &
                       length_neg_y, length_pos_z, length_neg_z, &
                       real_length_pos_x, real_length_neg_x, real_length_pos_y, &
                       real_length_neg_y, real_length_pos_z, real_length_neg_z
  use sourceprops, only: srcpos
  use array, only: xfinal, dt_array, ionized_array, ifront_plus_ionized_array, &
                   global_Ifront_array, source_Ifront_array, &
                   non_screened_Ifront_array, &
                   global_photoionization_HI_array, &
                   global_photoionization_HeI_array, global_photoionization_HeII_array, &
                   global_timestep_array, &
                   pos_x_pos_y_pos_z_octant, pos_x_pos_y_neg_z_octant, &
                   pos_x_neg_y_pos_z_octant, pos_x_neg_y_neg_z_octant, &
                   neg_x_pos_y_pos_z_octant, neg_x_pos_y_neg_z_octant, &
                   neg_x_neg_y_pos_z_octant, neg_x_neg_y_neg_z_octant, &
                   pos_x_pos_y_pos_z_dir_x_triangular_array, &
                   pos_x_pos_y_pos_z_dir_x_triangular_x_array, &
                   pos_x_pos_y_pos_z_dir_x_triangular_y_array, &
                   pos_x_pos_y_pos_z_dir_x_triangular_z_array, &
                   pos_x_pos_y_pos_z_dir_y_triangular_array, &
                   pos_x_pos_y_pos_z_dir_y_triangular_x_array, &
                   pos_x_pos_y_pos_z_dir_y_triangular_y_array, &
                   pos_x_pos_y_pos_z_dir_y_triangular_z_array, &
                   pos_x_pos_y_pos_z_dir_z_triangular_array, &
                   pos_x_pos_y_pos_z_dir_z_triangular_x_array, &
                   pos_x_pos_y_pos_z_dir_z_triangular_y_array, &
                   pos_x_pos_y_pos_z_dir_z_triangular_z_array, &
                   pos_x_pos_y_neg_z_dir_x_triangular_array, &
                   pos_x_pos_y_neg_z_dir_x_triangular_x_array, &
                   pos_x_pos_y_neg_z_dir_x_triangular_y_array, &
                   pos_x_pos_y_neg_z_dir_x_triangular_z_array, &
                   pos_x_pos_y_neg_z_dir_y_triangular_array, &
                   pos_x_pos_y_neg_z_dir_y_triangular_x_array, &
                   pos_x_pos_y_neg_z_dir_y_triangular_y_array, &
                   pos_x_pos_y_neg_z_dir_y_triangular_z_array, &
                   pos_x_pos_y_neg_z_dir_z_triangular_array, &
                   pos_x_pos_y_neg_z_dir_z_triangular_x_array, &
                   pos_x_pos_y_neg_z_dir_z_triangular_y_array, &
                   pos_x_pos_y_neg_z_dir_z_triangular_z_array, &
                   pos_x_neg_y_pos_z_dir_x_triangular_array, &
                   pos_x_neg_y_pos_z_dir_x_triangular_x_array, &
                   pos_x_neg_y_pos_z_dir_x_triangular_y_array, &
                   pos_x_neg_y_pos_z_dir_x_triangular_z_array, &
                   pos_x_neg_y_pos_z_dir_y_triangular_array, &
                   pos_x_neg_y_pos_z_dir_y_triangular_x_array, &
                   pos_x_neg_y_pos_z_dir_y_triangular_y_array, &
                   pos_x_neg_y_pos_z_dir_y_triangular_z_array, &
                   pos_x_neg_y_pos_z_dir_z_triangular_array, &
                   pos_x_neg_y_pos_z_dir_z_triangular_x_array, &
                   pos_x_neg_y_pos_z_dir_z_triangular_y_array, &
                   pos_x_neg_y_pos_z_dir_z_triangular_z_array, &
                   pos_x_neg_y_neg_z_dir_x_triangular_array, &
                   pos_x_neg_y_neg_z_dir_x_triangular_x_array, &
                   pos_x_neg_y_neg_z_dir_x_triangular_y_array, &
                   pos_x_neg_y_neg_z_dir_x_triangular_z_array, &
                   pos_x_neg_y_neg_z_dir_y_triangular_array, &
                   pos_x_neg_y_neg_z_dir_y_triangular_x_array, &
                   pos_x_neg_y_neg_z_dir_y_triangular_y_array, &
                   pos_x_neg_y_neg_z_dir_y_triangular_z_array, &
                   pos_x_neg_y_neg_z_dir_z_triangular_array, &
                   pos_x_neg_y_neg_z_dir_z_triangular_x_array, &
                   pos_x_neg_y_neg_z_dir_z_triangular_y_array, &
                   pos_x_neg_y_neg_z_dir_z_triangular_z_array, &
                   neg_x_pos_y_pos_z_dir_x_triangular_array, &
                   neg_x_pos_y_pos_z_dir_x_triangular_x_array, &
                   neg_x_pos_y_pos_z_dir_x_triangular_y_array, &
                   neg_x_pos_y_pos_z_dir_x_triangular_z_array, &
                   neg_x_pos_y_pos_z_dir_y_triangular_array, &
                   neg_x_pos_y_pos_z_dir_y_triangular_x_array, &
                   neg_x_pos_y_pos_z_dir_y_triangular_y_array, &
                   neg_x_pos_y_pos_z_dir_y_triangular_z_array, &
                   neg_x_pos_y_pos_z_dir_z_triangular_array, &
                   neg_x_pos_y_pos_z_dir_z_triangular_x_array, &
                   neg_x_pos_y_pos_z_dir_z_triangular_y_array, &
                   neg_x_pos_y_pos_z_dir_z_triangular_z_array, &
                   neg_x_pos_y_neg_z_dir_x_triangular_array, &
                   neg_x_pos_y_neg_z_dir_x_triangular_x_array, &
                   neg_x_pos_y_neg_z_dir_x_triangular_y_array, &
                   neg_x_pos_y_neg_z_dir_x_triangular_z_array, &
                   neg_x_pos_y_neg_z_dir_y_triangular_array, &
                   neg_x_pos_y_neg_z_dir_y_triangular_x_array, &
                   neg_x_pos_y_neg_z_dir_y_triangular_y_array, &
                   neg_x_pos_y_neg_z_dir_y_triangular_z_array, &
                   neg_x_pos_y_neg_z_dir_z_triangular_array, &
                   neg_x_pos_y_neg_z_dir_z_triangular_x_array, &
                   neg_x_pos_y_neg_z_dir_z_triangular_y_array, &
                   neg_x_pos_y_neg_z_dir_z_triangular_z_array, &
                   neg_x_neg_y_pos_z_dir_x_triangular_array, &
                   neg_x_neg_y_pos_z_dir_x_triangular_x_array, &
                   neg_x_neg_y_pos_z_dir_x_triangular_y_array, &
                   neg_x_neg_y_pos_z_dir_x_triangular_z_array, &
                   neg_x_neg_y_pos_z_dir_y_triangular_array, &
                   neg_x_neg_y_pos_z_dir_y_triangular_x_array, &
                   neg_x_neg_y_pos_z_dir_y_triangular_y_array, &
                   neg_x_neg_y_pos_z_dir_y_triangular_z_array, &
                   neg_x_neg_y_pos_z_dir_z_triangular_array, &
                   neg_x_neg_y_pos_z_dir_z_triangular_x_array, &
                   neg_x_neg_y_pos_z_dir_z_triangular_y_array, &
                   neg_x_neg_y_pos_z_dir_z_triangular_z_array, &
                   neg_x_neg_y_neg_z_dir_x_triangular_array, &
                   neg_x_neg_y_neg_z_dir_x_triangular_x_array, &
                   neg_x_neg_y_neg_z_dir_x_triangular_y_array, &
                   neg_x_neg_y_neg_z_dir_x_triangular_z_array, &
                   neg_x_neg_y_neg_z_dir_y_triangular_array, &
                   neg_x_neg_y_neg_z_dir_y_triangular_x_array, &
                   neg_x_neg_y_neg_z_dir_y_triangular_y_array, &
                   neg_x_neg_y_neg_z_dir_y_triangular_z_array, &
                   neg_x_neg_y_neg_z_dir_z_triangular_array, &
                   neg_x_neg_y_neg_z_dir_z_triangular_x_array, &
                   neg_x_neg_y_neg_z_dir_z_triangular_y_array, &
                   neg_x_neg_y_neg_z_dir_z_triangular_z_array
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AT_array_initialization_1()

    implicit none

    ! dummy variable of time partition
    integer :: i_partition  
                             
    do i_partition = 1,time_partition
      xfinal(i_partition) = real(i_partition)/real(time_partition+1)
    enddo

    allocate(dt_array(1:number_of_source))

  end subroutine AT_array_initialization_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AT_array_initialization_2(ns)

    implicit none

    integer, intent(in) :: ns

    if (do_periodic) then
      length_pos_x = mesh/2
      length_neg_x = -mesh/2+1
      length_pos_y = mesh/2
      length_neg_y = -mesh/2+1
      length_pos_z = mesh/2
      length_neg_z = -mesh/2+1
    else
      length_pos_x = mesh-srcpos(1,ns)
      length_neg_x = 1-srcpos(1,ns)
      length_pos_y = mesh-srcpos(2,ns)
      length_neg_y = 1-srcpos(2,ns)
      length_pos_z = mesh-srcpos(3,ns)
      length_neg_z = 1-srcpos(3,ns)
    endif

    real_length_pos_x = real(length_pos_x)
    real_length_neg_x = real(length_neg_x)
    real_length_pos_y = real(length_pos_y)
    real_length_neg_y = real(length_neg_y)
    real_length_pos_z = real(length_pos_z)
    real_length_neg_z = real(length_neg_z)

    allocate(source_Ifront_array(length_neg_x:length_pos_x,length_neg_y:length_pos_y,length_neg_z:length_pos_z))
    allocate(non_screened_Ifront_array(length_neg_x:length_pos_x,length_neg_y:length_pos_y,length_neg_z:length_pos_z))

    allocate(ionized_array(1:mesh,1:mesh,1:mesh))
    allocate(ifront_plus_ionized_array(1:mesh,1:mesh,1:mesh))

    allocate(global_photoionization_HI_array(1:mesh,1:mesh,1:mesh))
    allocate(global_photoionization_HeI_array(1:mesh,1:mesh,1:mesh))
    allocate(global_photoionization_HeII_array(1:mesh,1:mesh,1:mesh))
    allocate(global_timestep_array(1:mesh,1:mesh,1:mesh))

    allocate(pos_x_pos_y_pos_z_octant(0:length_pos_x,0:length_pos_y,0:length_pos_z))
    allocate(pos_x_pos_y_neg_z_octant(0:length_pos_x,0:length_pos_y,length_neg_z:0))
    allocate(pos_x_neg_y_pos_z_octant(0:length_pos_x,length_neg_y:0,0:length_pos_z))
    allocate(pos_x_neg_y_neg_z_octant(0:length_pos_x,length_neg_y:0,length_neg_z:0))
    allocate(neg_x_pos_y_pos_z_octant(length_neg_x:0,0:length_pos_y,0:length_pos_z))
    allocate(neg_x_pos_y_neg_z_octant(length_neg_x:0,0:length_pos_y,length_neg_z:0))
    allocate(neg_x_neg_y_pos_z_octant(length_neg_x:0,length_neg_y:0,0:length_pos_z))
    allocate(neg_x_neg_y_neg_z_octant(length_neg_x:0,length_neg_y:0,length_neg_z:0))

    global_Ifront_array = 0 
    source_Ifront_array = 0
    ionized_array = 0
    ifront_plus_ionized_array = 0
    non_screened_Ifront_array = 0
    global_photoionization_HI_array = 0
    global_photoionization_HeI_array = 0
    global_photoionization_HeII_array = 0
    global_timestep_array = -1.0
    pos_x_pos_y_pos_z_octant = 0
    pos_x_pos_y_neg_z_octant = 0
    pos_x_neg_y_pos_z_octant = 0
    pos_x_neg_y_neg_z_octant = 0
    neg_x_pos_y_pos_z_octant = 0
    neg_x_pos_y_neg_z_octant = 0
    neg_x_neg_y_pos_z_octant = 0
    neg_x_neg_y_neg_z_octant = 0
		
  end subroutine AT_array_initialization_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AT_dt_array_destruction()

    implicit none

    deallocate(dt_array)
  
  end subroutine AT_dt_array_destruction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AT_triangular_array_destruction()

    deallocate(pos_x_pos_y_pos_z_dir_x_triangular_array)
    deallocate(pos_x_pos_y_pos_z_dir_x_triangular_x_array)
    deallocate(pos_x_pos_y_pos_z_dir_x_triangular_y_array)
    deallocate(pos_x_pos_y_pos_z_dir_x_triangular_z_array)
    deallocate(pos_x_pos_y_pos_z_dir_y_triangular_array)
    deallocate(pos_x_pos_y_pos_z_dir_y_triangular_x_array)
    deallocate(pos_x_pos_y_pos_z_dir_y_triangular_y_array)
    deallocate(pos_x_pos_y_pos_z_dir_y_triangular_z_array)
    deallocate(pos_x_pos_y_pos_z_dir_z_triangular_array)
    deallocate(pos_x_pos_y_pos_z_dir_z_triangular_x_array)
    deallocate(pos_x_pos_y_pos_z_dir_z_triangular_y_array)
    deallocate(pos_x_pos_y_pos_z_dir_z_triangular_z_array)
    deallocate(pos_x_pos_y_neg_z_dir_x_triangular_array)
    deallocate(pos_x_pos_y_neg_z_dir_x_triangular_x_array)
    deallocate(pos_x_pos_y_neg_z_dir_x_triangular_y_array)
    deallocate(pos_x_pos_y_neg_z_dir_x_triangular_z_array)
    deallocate(pos_x_pos_y_neg_z_dir_y_triangular_array)
    deallocate(pos_x_pos_y_neg_z_dir_y_triangular_x_array)
    deallocate(pos_x_pos_y_neg_z_dir_y_triangular_y_array)
    deallocate(pos_x_pos_y_neg_z_dir_y_triangular_z_array)
    deallocate(pos_x_pos_y_neg_z_dir_z_triangular_array)
    deallocate(pos_x_pos_y_neg_z_dir_z_triangular_x_array)
    deallocate(pos_x_pos_y_neg_z_dir_z_triangular_y_array)
    deallocate(pos_x_pos_y_neg_z_dir_z_triangular_z_array)
    deallocate(pos_x_neg_y_pos_z_dir_x_triangular_array)
    deallocate(pos_x_neg_y_pos_z_dir_x_triangular_x_array)
    deallocate(pos_x_neg_y_pos_z_dir_x_triangular_y_array)
    deallocate(pos_x_neg_y_pos_z_dir_x_triangular_z_array)
    deallocate(pos_x_neg_y_pos_z_dir_y_triangular_array)
    deallocate(pos_x_neg_y_pos_z_dir_y_triangular_x_array)
    deallocate(pos_x_neg_y_pos_z_dir_y_triangular_y_array)
    deallocate(pos_x_neg_y_pos_z_dir_y_triangular_z_array)
    deallocate(pos_x_neg_y_pos_z_dir_z_triangular_array)
    deallocate(pos_x_neg_y_pos_z_dir_z_triangular_x_array)
    deallocate(pos_x_neg_y_pos_z_dir_z_triangular_y_array)
    deallocate(pos_x_neg_y_pos_z_dir_z_triangular_z_array)
    deallocate(pos_x_neg_y_neg_z_dir_x_triangular_array)
    deallocate(pos_x_neg_y_neg_z_dir_x_triangular_x_array)
    deallocate(pos_x_neg_y_neg_z_dir_x_triangular_y_array)
    deallocate(pos_x_neg_y_neg_z_dir_x_triangular_z_array)
    deallocate(pos_x_neg_y_neg_z_dir_y_triangular_array)
    deallocate(pos_x_neg_y_neg_z_dir_y_triangular_x_array)
    deallocate(pos_x_neg_y_neg_z_dir_y_triangular_y_array)
    deallocate(pos_x_neg_y_neg_z_dir_y_triangular_z_array)
    deallocate(pos_x_neg_y_neg_z_dir_z_triangular_array)
    deallocate(pos_x_neg_y_neg_z_dir_z_triangular_x_array)
    deallocate(pos_x_neg_y_neg_z_dir_z_triangular_y_array)
    deallocate(pos_x_neg_y_neg_z_dir_z_triangular_z_array)
    deallocate(neg_x_pos_y_pos_z_dir_x_triangular_array)
    deallocate(neg_x_pos_y_pos_z_dir_x_triangular_x_array)
    deallocate(neg_x_pos_y_pos_z_dir_x_triangular_y_array)
    deallocate(neg_x_pos_y_pos_z_dir_x_triangular_z_array)
    deallocate(neg_x_pos_y_pos_z_dir_y_triangular_array)
    deallocate(neg_x_pos_y_pos_z_dir_y_triangular_x_array)
    deallocate(neg_x_pos_y_pos_z_dir_y_triangular_y_array)
    deallocate(neg_x_pos_y_pos_z_dir_y_triangular_z_array)
    deallocate(neg_x_pos_y_pos_z_dir_z_triangular_array)
    deallocate(neg_x_pos_y_pos_z_dir_z_triangular_x_array)
    deallocate(neg_x_pos_y_pos_z_dir_z_triangular_y_array)
    deallocate(neg_x_pos_y_pos_z_dir_z_triangular_z_array)
    deallocate(neg_x_pos_y_neg_z_dir_x_triangular_array)
    deallocate(neg_x_pos_y_neg_z_dir_x_triangular_x_array)
    deallocate(neg_x_pos_y_neg_z_dir_x_triangular_y_array)
    deallocate(neg_x_pos_y_neg_z_dir_x_triangular_z_array)
    deallocate(neg_x_pos_y_neg_z_dir_y_triangular_array)
    deallocate(neg_x_pos_y_neg_z_dir_y_triangular_x_array)
    deallocate(neg_x_pos_y_neg_z_dir_y_triangular_y_array)
    deallocate(neg_x_pos_y_neg_z_dir_y_triangular_z_array)
    deallocate(neg_x_pos_y_neg_z_dir_z_triangular_array)
    deallocate(neg_x_pos_y_neg_z_dir_z_triangular_x_array)
    deallocate(neg_x_pos_y_neg_z_dir_z_triangular_y_array)
    deallocate(neg_x_pos_y_neg_z_dir_z_triangular_z_array)
    deallocate(neg_x_neg_y_pos_z_dir_x_triangular_array)
    deallocate(neg_x_neg_y_pos_z_dir_x_triangular_x_array)
    deallocate(neg_x_neg_y_pos_z_dir_x_triangular_y_array)
    deallocate(neg_x_neg_y_pos_z_dir_x_triangular_z_array)
    deallocate(neg_x_neg_y_pos_z_dir_y_triangular_array)
    deallocate(neg_x_neg_y_pos_z_dir_y_triangular_x_array)
    deallocate(neg_x_neg_y_pos_z_dir_y_triangular_y_array)
    deallocate(neg_x_neg_y_pos_z_dir_y_triangular_z_array)
    deallocate(neg_x_neg_y_pos_z_dir_z_triangular_array)
    deallocate(neg_x_neg_y_pos_z_dir_z_triangular_x_array)
    deallocate(neg_x_neg_y_pos_z_dir_z_triangular_y_array)
    deallocate(neg_x_neg_y_pos_z_dir_z_triangular_z_array)
    deallocate(neg_x_neg_y_neg_z_dir_x_triangular_array)
    deallocate(neg_x_neg_y_neg_z_dir_x_triangular_x_array)
    deallocate(neg_x_neg_y_neg_z_dir_x_triangular_y_array)
    deallocate(neg_x_neg_y_neg_z_dir_x_triangular_z_array)
    deallocate(neg_x_neg_y_neg_z_dir_y_triangular_array)
    deallocate(neg_x_neg_y_neg_z_dir_y_triangular_x_array)
    deallocate(neg_x_neg_y_neg_z_dir_y_triangular_y_array)
    deallocate(neg_x_neg_y_neg_z_dir_y_triangular_z_array)
    deallocate(neg_x_neg_y_neg_z_dir_z_triangular_array)
    deallocate(neg_x_neg_y_neg_z_dir_z_triangular_x_array)
    deallocate(neg_x_neg_y_neg_z_dir_z_triangular_y_array)
    deallocate(neg_x_neg_y_neg_z_dir_z_triangular_z_array)

  end subroutine AT_triangular_array_destruction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AT_general_array_destruction()

    implicit none

    deallocate(source_Ifront_array)
    deallocate(ionized_array)
    deallocate(ifront_plus_ionized_array)
    deallocate(non_screened_Ifront_array)
    deallocate(global_photoionization_HI_array)
    deallocate(global_photoionization_HeI_array)
    deallocate(global_photoionization_HeII_array)
    deallocate(global_timestep_array)
    deallocate(pos_x_pos_y_pos_z_octant)
    deallocate(pos_x_pos_y_neg_z_octant)
    deallocate(pos_x_neg_y_pos_z_octant)
    deallocate(pos_x_neg_y_neg_z_octant)
    deallocate(neg_x_pos_y_pos_z_octant)
    deallocate(neg_x_pos_y_neg_z_octant)
    deallocate(neg_x_neg_y_pos_z_octant)
    deallocate(neg_x_neg_y_neg_z_octant)
    
  end subroutine AT_general_array_destruction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module AT_array
