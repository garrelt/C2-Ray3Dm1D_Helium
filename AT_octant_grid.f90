module AT_octant_grid

  use parameter, only: mesh, &
                       length_pos_x, length_neg_x, length_pos_y, &
                       length_neg_y, length_pos_z, length_neg_z
  use my_mpi
  use array, only: source_Ifront_array, &
                   pos_x_pos_y_pos_z_octant, pos_x_pos_y_neg_z_octant, &
                   pos_x_neg_y_pos_z_octant, pos_x_neg_y_neg_z_octant, &
                   neg_x_pos_y_pos_z_octant, neg_x_pos_y_neg_z_octant, &
                   neg_x_neg_y_pos_z_octant, neg_x_neg_y_neg_z_octant

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AT_octant_grid_creation
	
    implicit none
	
    integer :: i, j, k

!$OMP PARALLEL
	
!$OMP SECTIONS PRIVATE(i,j,k)

!$OMP SECTION

    do i = 0,length_pos_x 
      do j = 0,length_pos_y
        do k = 0,length_pos_z
          pos_x_pos_y_pos_z_octant(i,j,k) = source_Ifront_array(i,j,k)
        enddo
      enddo
    enddo

!$OMP SECTION		
    do i = 0,length_pos_x
      do j = 0,length_pos_y
        do k = 0,length_neg_z,-1
          pos_x_pos_y_neg_z_octant(i,j,k) = source_Ifront_array(i,j,k)
        enddo
      enddo
    enddo
!$OMP SECTION		
    do i = 0,length_pos_x
      do j = 0,length_neg_y,-1
        do k = 0,length_pos_z
          pos_x_neg_y_pos_z_octant(i,j,k) = source_Ifront_array(i,j,k)
        enddo
      enddo
    enddo
!$OMP SECTION		
    do i = 0,length_pos_x
      do j = 0,length_neg_y,-1
        do k = 0,length_neg_z,-1
          pos_x_neg_y_neg_z_octant(i,j,k) = source_Ifront_array(i,j,k)
        enddo
      enddo
    enddo
!$OMP SECTION		
    do i = 0,length_neg_x,-1
      do j = 0,length_pos_y
        do k = 0,length_pos_z
          neg_x_pos_y_pos_z_octant(i,j,k) = source_Ifront_array(i,j,k)
        enddo
      enddo
    enddo
!$OMP SECTION		
    do i = 0,length_neg_x,-1
      do j = 0,length_pos_y
        do k = 0,length_neg_z,-1
          neg_x_pos_y_neg_z_octant(i,j,k) = source_Ifront_array(i,j,k)
        enddo
      enddo
    enddo
!$OMP SECTION		
    do i = 0,length_neg_x,-1
      do j = 0,length_neg_y,-1
        do k = 0,length_pos_z
          neg_x_neg_y_pos_z_octant(i,j,k) = source_Ifront_array(i,j,k)
        enddo
      enddo
    enddo
!$OMP SECTION		
    do i = 0,length_neg_x,-1
      do j = 0,length_neg_y,-1
        do k = 0,length_neg_z,-1
          neg_x_neg_y_neg_z_octant(i,j,k) = source_Ifront_array(i,j,k)
        enddo
      enddo
    enddo
!$OMP END SECTIONS
	
!$OMP END PARALLEL
	
  end subroutine AT_octant_grid_creation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module AT_octant_grid
