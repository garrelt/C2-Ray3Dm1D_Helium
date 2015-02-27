module AT_ray
  
  use my_mpi
  use file_admin, only: logf
  use parameter, only: mesh, half, &
                       real_length_pos_x, real_length_neg_x, &
                       real_length_pos_y, real_length_neg_y, &
                       real_length_pos_z, real_length_neg_z
  use array, only: non_screened_Ifront_array, &
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
                   neg_x_neg_y_neg_z_dir_z_triangular_z_array, &
                   source_Ifront_array

  implicit none 
	
  contains
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AT_ray_pass_through_octant(i_source)
	
    implicit none

    integer, intent(in) :: i_source

!$OMP PARALLEL

!$OMP SECTIONS 
!$OMP SECTION
    call ray_pass_through_pos_x_pos_y_pos_z_octant(i_source)
!$OMP SECTION
    call ray_pass_through_pos_x_pos_y_neg_z_octant(i_source)
!$OMP SECTION
    call ray_pass_through_pos_x_neg_y_pos_z_octant(i_source)
!$OMP SECTION
    call ray_pass_through_pos_x_neg_y_neg_z_octant(i_source)
!$OMP SECTION
    call ray_pass_through_neg_x_pos_y_pos_z_octant(i_source)
!$OMP SECTION
    call ray_pass_through_neg_x_pos_y_neg_z_octant(i_source)
!$OMP SECTION
    call ray_pass_through_neg_x_neg_y_pos_z_octant(i_source)
!$OMP SECTION
    call ray_pass_through_neg_x_neg_y_neg_z_octant(i_source)
!$OMP END SECTIONS

!$OMP END PARALLEL

    call update_non_screened_Ifront_array()

  end subroutine AT_ray_pass_through_octant

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
  subroutine ray_pass_through_pos_x_pos_y_pos_z_octant(i_source)

    implicit none
		
    integer, intent(in) :: i_source

    real :: ii,jj,kk
    logical :: pass
    integer :: passed_cell

    !direction x triangular array
    ii = real_length_pos_x
    jj = 0.0
    kk = 0.0	
	
    do while (abs(jj).le.abs(real_length_pos_x)) 
      kk = 0.0	
       do while (abs(kk).le.abs(real_length_pos_x))
        call find_non_screened_Ifront_cells('+','x','+','y','+','z',ii,kk,jj,'x', &
                                               size(pos_x_pos_y_pos_z_dir_x_triangular_array))
      enddo
      call update_triangular_array('+','x','+','y','+','z',jj,'x',size(pos_x_pos_y_pos_z_dir_x_triangular_array))
      jj = jj + 1.0
    enddo
	
    !direction y triangular array
    ii = 0.0
    jj = real_length_pos_y
    kk = 0.0		
    do while (abs(kk).le.abs(real_length_pos_y))
      ii = 0.0	
      do while (abs(ii).le.abs(real_length_pos_y))	
        call find_non_screened_Ifront_cells('+','x','+','y','+','z',jj,ii,kk,'y', &
                                               size(pos_x_pos_y_pos_z_dir_y_triangular_array))
      enddo
      call update_triangular_array('+','x','+','y','+','z',kk,'y',size(pos_x_pos_y_pos_z_dir_y_triangular_array))
      kk = kk + 1.0
    enddo

    !direction z triangular array
    ii = 0
    jj = 0
    kk = real_length_pos_z	
    do while (abs(ii).le.abs(real_length_pos_z))
      jj = 0.0	
      do while (abs(jj).le.abs(real_length_pos_z))	
        call find_non_screened_Ifront_cells('+','x','+','y','+','z',kk,jj,ii,'z', &
                                               size(pos_x_pos_y_pos_z_dir_z_triangular_array))
      enddo
      call update_triangular_array('+','x','+','y','+','z',ii,'z',size(pos_x_pos_y_pos_z_dir_z_triangular_array))
      ii = ii + 1.0
    enddo

  end subroutine ray_pass_through_pos_x_pos_y_pos_z_octant
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
  subroutine ray_pass_through_pos_x_pos_y_neg_z_octant(i_source)
		
    implicit none
		
    integer, intent(in) :: i_source

    real :: ii,jj,kk 
    logical :: pass
    integer :: passed_cell
	
    !direction x triangular array
    ii = real_length_pos_x 
    jj = 0.0
    kk = 0.0		
    do while (abs(jj).le.abs(real_length_pos_x))
      kk = 0.0	
      do while (abs(kk).le.abs(real_length_pos_x))	
        call find_non_screened_Ifront_cells('+','x','+','y','-','z',ii,kk,jj,'x', &
                                               size(pos_x_pos_y_neg_z_dir_x_triangular_array))
      enddo
      call update_triangular_array('+','x','+','y','-','z',jj,'x',size(pos_x_pos_y_neg_z_dir_x_triangular_array))
      jj = jj + 1.0
    enddo
	
    !direction y triangular array
    ii = 0.0
    jj = real_length_pos_y
    kk = 0.0		
    do while (abs(kk).le.abs(real_length_pos_y))   
    ii = 0.0	
      do while (abs(ii).le.abs(real_length_pos_y))	
        call find_non_screened_Ifront_cells('+','x','+','y','-','z',jj,ii,kk,'y', &
                                               size(pos_x_pos_y_neg_z_dir_y_triangular_array))
      enddo
      call update_triangular_array('+','x','+','y','-','z',kk,'y',size(pos_x_pos_y_neg_z_dir_y_triangular_array))
      kk = kk - 1.0
    enddo
		
    !direction z triangular array
    ii = 0.0
    jj = 0.0
    kk = real_length_neg_z	
    do while (abs(ii).le.abs(real_length_pos_z))
      jj = 0.0	
      do while (abs(jj).le.abs(real_length_pos_z))	
        call find_non_screened_Ifront_cells('+','x','+','y','-','z',kk,jj,ii,'z', &
                                               size(pos_x_pos_y_neg_z_dir_z_triangular_array))
      enddo
      call update_triangular_array('+','x','+','y','-','z',ii,'z',size(pos_x_pos_y_neg_z_dir_z_triangular_array))
      ii = ii + 1.0
    enddo

  end subroutine ray_pass_through_pos_x_pos_y_neg_z_octant
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ray_pass_through_pos_x_neg_y_pos_z_octant(i_source)
	
    implicit none
		
    integer, intent(in) :: i_source
	
    real :: ii,jj,kk
    logical :: pass
    integer :: passed_cell
			
    !direction x triangular array
    ii = real_length_pos_x
    jj = 0.0
    kk = 0.0		
    do while (abs(jj).le.abs(real_length_pos_x))
      kk = 0.0	
      do while (abs(kk).le.abs(real_length_pos_x))	
        call find_non_screened_Ifront_cells('+','x','-','y','+','z',ii,kk,jj,'x', &
                                               size(pos_x_neg_y_pos_z_dir_x_triangular_array))
      enddo
      call update_triangular_array('+','x','-','y','+','z',jj,'x',size(pos_x_neg_y_pos_z_dir_x_triangular_array))
      jj = jj - 1.0
    enddo

    !direction y triangular array
    ii = 0.0
    jj = real_length_neg_y
    kk = 0.0		
    do while (abs(kk).le.abs(real_length_pos_y))
      ii = 0.0	
      do while (abs(ii).le.abs(real_length_pos_y))	
        call find_non_screened_Ifront_cells('+','x','-','y','+','z',jj,ii,kk,'y', &
                                               size(pos_x_neg_y_pos_z_dir_y_triangular_array))
      enddo
      call update_triangular_array('+','x','-','y','+','z',kk,'y',size(pos_x_neg_y_pos_z_dir_y_triangular_array))
      kk = kk + 1.0
    enddo
	
    !direction z triangular array
    ii = 0.0
    jj = 0.0
    kk = real_length_pos_z	
    do while (abs(ii).le.abs(real_length_pos_z))
      jj = 0.0	
      do while (abs(jj).le.abs(real_length_pos_z))	
        call find_non_screened_Ifront_cells('+','x','-','y','+','z',kk,jj,ii,'z', &
                                               size(pos_x_neg_y_pos_z_dir_z_triangular_array))
      enddo
      call update_triangular_array('+','x','-','y','+','z',ii,'z',size(pos_x_neg_y_pos_z_dir_z_triangular_array))
      ii = ii + 1.0
    enddo

  end subroutine ray_pass_through_pos_x_neg_y_pos_z_octant
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
  subroutine ray_pass_through_pos_x_neg_y_neg_z_octant(i_source)
	
    implicit none
		
    integer, intent(in) :: i_source
	
    real :: ii,jj,kk 
    logical :: pass
    integer :: passed_cell
			
    !direction x triangular array
    ii = real_length_pos_x
    jj = 0.0
    kk = 0.0		
    do while (abs(jj).le.abs(real_length_pos_x))
      kk = 0.0	
      do while (abs(kk).le.abs(real_length_pos_x))	
        call find_non_screened_Ifront_cells('+','x','-','y','-','z',ii,kk,jj,'x', &
                                               size(pos_x_neg_y_neg_z_dir_x_triangular_array))
      enddo
      call update_triangular_array('+','x','-','y','-','z',jj,'x',size(pos_x_neg_y_neg_z_dir_x_triangular_array))
      jj = jj - 1.0
    enddo
	
    !direction y triangular array
    ii = 0.0
    jj = real_length_neg_y
    kk = 0.0		
    do while (abs(kk).le.abs(real_length_pos_y))
      ii = 0.0	
      do while (abs(ii).le.abs(real_length_pos_y))
        call find_non_screened_Ifront_cells('+','x','-','y','-','z',jj,ii,kk,'y', &
                                               size(pos_x_neg_y_neg_z_dir_y_triangular_array))
      enddo
      call update_triangular_array('+','x','-','y','-','z',kk,'y',size(pos_x_neg_y_neg_z_dir_y_triangular_array))
      kk = kk - 1.0
    enddo
	
    !direction z triangular array
    ii = 0.0
    jj = 0.0
    kk = real_length_neg_z	
    do while (abs(ii).le.abs(real_length_pos_z))
      jj= 0.0	
      do while (abs(jj).le.abs(real_length_pos_z))	
        call find_non_screened_Ifront_cells('+','x','-','y','-','z',kk,jj,ii,'z', &
                                               size(pos_x_neg_y_neg_z_dir_z_triangular_array))
      enddo
      call update_triangular_array('+','x','-','y','-','z',ii,'z',size(pos_x_neg_y_neg_z_dir_z_triangular_array))
      ii = ii + 1.0
    enddo

  end subroutine ray_pass_through_pos_x_neg_y_neg_z_octant
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ray_pass_through_neg_x_pos_y_pos_z_octant(i_source)
	
    implicit none
		
    integer, intent(in) :: i_source
	
    real :: ii,jj,kk
    logical :: pass
    integer :: passed_cell
			
    !direction x triangular array
    ii = real_length_neg_x
    jj = 0.0
    kk = 0.0		
    do while (abs(jj).le.abs(real_length_pos_x))
      kk = 0.0	
      do while (abs(kk).le.abs(real_length_pos_x))	
        call find_non_screened_Ifront_cells('-','x','+','y','+','z',ii,kk,jj,'x', &
                                               size(neg_x_pos_y_pos_z_dir_x_triangular_array))
      enddo
      call update_triangular_array('-','x','+','y','+','z',jj,'x',size(neg_x_pos_y_pos_z_dir_x_triangular_array))
      jj = jj + 1.0
    enddo
		
    !direction y triangular array
    ii = 0.0
    jj = real_length_pos_y
    kk = 0.0		
    do while (abs(kk).le.abs(real_length_pos_y))
      ii = 0.0	
      do while (abs(ii).le.abs(real_length_pos_y))	
        call find_non_screened_Ifront_cells('-','x','+','y','+','z',jj,ii,kk,'y', &
                                               size(neg_x_pos_y_pos_z_dir_y_triangular_array))
      enddo
      call update_triangular_array('-','x','+','y','+','z',kk,'y',size(neg_x_pos_y_pos_z_dir_y_triangular_array))
      kk = kk + 1.0
    enddo
		
    !direction z triangular array
    ii = 0.0
    jj = 0.0
    kk = real_length_pos_z	
    do while (abs(ii).le.abs(real_length_pos_z))
      jj = 0.0	
      do while (abs(jj).le.abs(real_length_pos_z))	
        call find_non_screened_Ifront_cells('-','x','+','y','+','z',kk,jj,ii,'z', &
                                               size(neg_x_pos_y_pos_z_dir_z_triangular_array))
      enddo
      call update_triangular_array('-','x','+','y','+','z',ii,'z',size(neg_x_pos_y_pos_z_dir_z_triangular_array))
      ii = ii - 1.0
    enddo

  end subroutine ray_pass_through_neg_x_pos_y_pos_z_octant
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
  subroutine ray_pass_through_neg_x_pos_y_neg_z_octant(i_source)
	
    implicit none
		
    integer, intent(in) :: i_source
	
    real :: ii,jj,kk 
    logical :: pass
    integer :: passed_cell
		
    !direction x triangular array
    ii = real_length_neg_x
    jj = 0.0
    kk = 0.0		
    do while (abs(jj).le.abs(real_length_pos_x))
      kk = 0.0	
      do while (abs(kk).le.abs(real_length_pos_x))	
        call find_non_screened_Ifront_cells('-','x','+','y','-','z',ii,kk,jj,'x', &
                                               size(neg_x_pos_y_neg_z_dir_x_triangular_array))
      enddo
      call update_triangular_array('-','x','+','y','-','z',jj,'x',size(neg_x_pos_y_neg_z_dir_x_triangular_array))
      jj = jj + 1.0
    enddo
		
    !direction y triangular array
    ii = 0.0
    jj = real_length_pos_y
    kk = 0.0		
    do while (abs(kk).le.abs(real_length_pos_y))
      ii = 0.0	
      do while (abs(ii).le.abs(real_length_pos_y))	
        call find_non_screened_Ifront_cells('-','x','+','y','-','z',jj,ii,kk,'y', &
                                               size(neg_x_pos_y_neg_z_dir_y_triangular_array))
      enddo
      call update_triangular_array('-','x','+','y','-','z',kk,'y',size(neg_x_pos_y_neg_z_dir_y_triangular_array))
      kk = kk - 1.0
    enddo
	
    !direction z triangular array
    ii = 0.0
    jj = 0.0
    kk = real_length_neg_z	
    do while (abs(ii).le.abs(real_length_pos_z))
      jj = 0.0	
      do while (abs(jj).le.abs(real_length_pos_z))	
        call find_non_screened_Ifront_cells('-','x','+','y','-','z',kk,jj,ii,'z', &
                                               size(neg_x_pos_y_neg_z_dir_z_triangular_array))
      enddo
      call update_triangular_array('-','x','+','y','-','z',ii,'z',size(neg_x_pos_y_neg_z_dir_z_triangular_array))
      ii = ii - 1.0
    enddo

  end subroutine ray_pass_through_neg_x_pos_y_neg_z_octant
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ray_pass_through_neg_x_neg_y_pos_z_octant(i_source)
	
    implicit none
		
    integer, intent(in) :: i_source
	
    real :: ii,jj,kk 
    logical :: pass
    integer :: passed_cell
		
    !direction x triangular array
    ii = real_length_neg_x
    jj = 0.0
    kk = 0.0		 
    do while (abs(jj).le.abs(real_length_pos_x))
      kk = 0.0	
      do while (abs(kk).le.abs(real_length_pos_x))	
        call find_non_screened_Ifront_cells('-','x','-','y','+','z',ii,kk,jj,'x', &
                                               size(neg_x_neg_y_pos_z_dir_x_triangular_array))
      enddo
      call update_triangular_array('-','x','-','y','+','z',jj,'x',size(neg_x_neg_y_pos_z_dir_x_triangular_array))
      jj = jj - 1.0
    enddo
	
    !direction y triangular array
    ii = 0.0
    jj = real_length_neg_y
    kk = 0.0		
    do while (abs(kk).le.abs(real_length_pos_y))
      ii = 0.0	
      do while (abs(ii).le.abs(real_length_pos_y))	
        call find_non_screened_Ifront_cells('-','x','-','y','+','z',jj,ii,kk,'y', &
                                               size(neg_x_neg_y_pos_z_dir_y_triangular_array))
      enddo
      call update_triangular_array('-','x','-','y','+','z',kk,'y',size(neg_x_neg_y_pos_z_dir_y_triangular_array))
      kk = kk + 1.0
    enddo

    !direction z triangular array
    ii = 0.0
    jj = 0.0
    kk = real_length_pos_z	
    do while (abs(ii).le.abs(real_length_pos_z))
      jj = 0.0	
      do while (abs(jj).le.abs(real_length_pos_z))	
        call find_non_screened_Ifront_cells('-','x','-','y','+','z',kk,jj,ii,'z', &
                                               size(neg_x_neg_y_pos_z_dir_z_triangular_array))
      enddo
      call update_triangular_array('-','x','-','y','+','z',ii,'z',size(neg_x_neg_y_pos_z_dir_z_triangular_array))
      ii = ii - 1.0
    enddo

  end subroutine ray_pass_through_neg_x_neg_y_pos_z_octant
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
  subroutine ray_pass_through_neg_x_neg_y_neg_z_octant(i_source)
	
    implicit none
			
    integer, intent(in) :: i_source

    real :: ii,jj,kk 
    logical :: pass
    integer :: passed_cell
			
    !direction x triangular array
    ii = real_length_neg_x
    jj = 0.0
    kk = 0.0		
    do while (abs(jj).le.abs(real_length_pos_x))
      kk = 0.0	
      do while (abs(kk).le.abs(real_length_pos_x))	
        call find_non_screened_Ifront_cells('-','x','-','y','-','z',ii,kk,jj,'x', &
                                               size(neg_x_neg_y_neg_z_dir_x_triangular_array))
      enddo
      call update_triangular_array('-','x','-','y','-','z',jj,'x',size(neg_x_neg_y_neg_z_dir_x_triangular_array))
      jj = jj - 1.0
    enddo
		
    !direction y triangular array
    ii = 0.0
    jj = real_length_neg_y
    kk = 0.0		
    do while (abs(kk).le.abs(real_length_pos_y))
      ii = 0.0	
      do while (abs(ii).le.abs(real_length_pos_y))	
        call find_non_screened_Ifront_cells('-','x','-','y','-','z',jj,ii,kk,'y', &
                                               size(neg_x_neg_y_neg_z_dir_y_triangular_array))
      enddo
      call update_triangular_array('-','x','-','y','-','z',kk,'y',size(neg_x_neg_y_neg_z_dir_y_triangular_array))
      kk = kk - 1.0
    enddo
	
    !direction z triangular array
    ii = 0.0
    jj = 0.0
    kk = real_length_neg_z	
    do while (abs(ii).le.abs(real_length_pos_z))
      jj= 0.0	
      do while (abs(jj).le.abs(real_length_pos_z))	
        call find_non_screened_Ifront_cells('-','x','-','y','-','z',kk,jj,ii,'z', &
                                               size(neg_x_neg_y_neg_z_dir_z_triangular_array))
      enddo
      call update_triangular_array('-','x','-','y','-','z',ii,'z',size(neg_x_neg_y_neg_z_dir_z_triangular_array))
      ii = ii - 1.0
    enddo

  end subroutine ray_pass_through_neg_x_neg_y_neg_z_octant
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine find_non_screened_Ifront_cells(parity1,axis1,parity2,axis2,parity3,axis3, &
                                               aa,bb,cc,direction,length)
	
    implicit none
	
    character,intent(in) :: parity1
    character,intent(in) :: axis1
    character,intent(in) :: parity2
    character,intent(in) :: axis2
    character,intent(in) :: parity3
    character,intent(in) :: axis3
    real,intent(in) :: aa !primary direction
    real,intent(inout) :: bb !secondary direction
    real,intent(in) :: cc !tertiary direction
    character,intent(in) :: direction
    integer,intent(in) :: length
    logical :: pass
    integer,dimension(:),pointer :: triangular_array
    integer,dimension(:),pointer :: triangular_a_array
    integer,dimension(:),pointer :: triangular_b_array
    integer,dimension(:),pointer :: triangular_c_array
    integer :: count
    logical :: pass_b
    real :: parity_a, parity_b, parity_c
    real :: a1,a2
    real :: b1,b2
    real :: c1,c2
		
    if (length.ne.0) then

      if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
          parity2.eq.'+' .and. axis2.eq.'y' .and. &
          parity3.eq.'+' .and. axis3.eq.'z' .and. &
          direction.eq.'x') then
         triangular_array => pos_x_pos_y_pos_z_dir_x_triangular_array
         triangular_a_array => pos_x_pos_y_pos_z_dir_x_triangular_x_array
         triangular_b_array => pos_x_pos_y_pos_z_dir_x_triangular_z_array
         triangular_c_array => pos_x_pos_y_pos_z_dir_x_triangular_y_array
      endif
		
      if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
          parity2.eq.'+' .and. axis2.eq.'y' .and. &
          parity3.eq.'+' .and. axis3.eq.'z' .and. &
          direction.eq.'y') then
         triangular_array => pos_x_pos_y_pos_z_dir_y_triangular_array
         triangular_a_array => pos_x_pos_y_pos_z_dir_y_triangular_y_array
         triangular_b_array => pos_x_pos_y_pos_z_dir_y_triangular_x_array
         triangular_c_array => pos_x_pos_y_pos_z_dir_y_triangular_z_array
      endif
		
      if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
          parity2.eq.'+' .and. axis2.eq.'y' .and. &
          parity3.eq.'+' .and. axis3.eq.'z' .and. &
          direction.eq.'z') then
         triangular_array => pos_x_pos_y_pos_z_dir_z_triangular_array
         triangular_a_array => pos_x_pos_y_pos_z_dir_z_triangular_z_array
         triangular_b_array => pos_x_pos_y_pos_z_dir_z_triangular_y_array
         triangular_c_array => pos_x_pos_y_pos_z_dir_z_triangular_x_array
      endif
	
      if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
          parity2.eq.'+' .and. axis2.eq.'y' .and. &
          parity3.eq.'-' .and. axis3.eq.'z' .and. &
          direction.eq.'x') then
         triangular_array => pos_x_pos_y_neg_z_dir_x_triangular_array
         triangular_a_array => pos_x_pos_y_neg_z_dir_x_triangular_x_array
         triangular_b_array => pos_x_pos_y_neg_z_dir_x_triangular_z_array
         triangular_c_array => pos_x_pos_y_neg_z_dir_x_triangular_y_array
      endif
		
      if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
          parity2.eq.'+' .and. axis2.eq.'y' .and. &
          parity3.eq.'-' .and. axis3.eq.'z' .and. &
          direction.eq.'y') then
         triangular_array => pos_x_pos_y_neg_z_dir_y_triangular_array
         triangular_a_array => pos_x_pos_y_neg_z_dir_y_triangular_y_array
         triangular_b_array => pos_x_pos_y_neg_z_dir_y_triangular_x_array
         triangular_c_array => pos_x_pos_y_neg_z_dir_y_triangular_z_array
      endif
		
      if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
          parity2.eq.'+' .and. axis2.eq.'y' .and. &
          parity3.eq.'-' .and. axis3.eq.'z' .and. &
          direction.eq.'z') then
         triangular_array => pos_x_pos_y_neg_z_dir_z_triangular_array
         triangular_a_array => pos_x_pos_y_neg_z_dir_z_triangular_z_array
         triangular_b_array => pos_x_pos_y_neg_z_dir_z_triangular_y_array
         triangular_c_array => pos_x_pos_y_neg_z_dir_z_triangular_x_array
      endif

      if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
          parity2.eq.'-' .and. axis2.eq.'y' .and. &
          parity3.eq.'+' .and. axis3.eq.'z' .and. &
          direction.eq.'x') then
         triangular_array => pos_x_neg_y_pos_z_dir_x_triangular_array
         triangular_a_array => pos_x_neg_y_pos_z_dir_x_triangular_x_array
         triangular_b_array => pos_x_neg_y_pos_z_dir_x_triangular_z_array
         triangular_c_array => pos_x_neg_y_pos_z_dir_x_triangular_y_array
      endif
		
      if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
          parity2.eq.'-' .and. axis2.eq.'y' .and. &
          parity3.eq.'+' .and. axis3.eq.'z' .and. &
          direction.eq.'y') then
         triangular_array => pos_x_neg_y_pos_z_dir_y_triangular_array
         triangular_a_array => pos_x_neg_y_pos_z_dir_y_triangular_y_array
         triangular_b_array => pos_x_neg_y_pos_z_dir_y_triangular_x_array
         triangular_c_array => pos_x_neg_y_pos_z_dir_y_triangular_z_array
      endif
		
      if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
          parity2.eq.'-' .and. axis2.eq.'y' .and. &
          parity3.eq.'+' .and. axis3.eq.'z' .and. &
          direction.eq.'z') then
         triangular_array => pos_x_neg_y_pos_z_dir_z_triangular_array
         triangular_a_array => pos_x_neg_y_pos_z_dir_z_triangular_z_array
         triangular_b_array => pos_x_neg_y_pos_z_dir_z_triangular_y_array
         triangular_c_array => pos_x_neg_y_pos_z_dir_z_triangular_x_array
      endif
	
      if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
          parity2.eq.'-' .and. axis2.eq.'y' .and. &
          parity3.eq.'-' .and. axis3.eq.'z' .and. &
          direction.eq.'x') then
         triangular_array => pos_x_neg_y_neg_z_dir_x_triangular_array
         triangular_a_array => pos_x_neg_y_neg_z_dir_x_triangular_x_array
         triangular_b_array => pos_x_neg_y_neg_z_dir_x_triangular_z_array
         triangular_c_array => pos_x_neg_y_neg_z_dir_x_triangular_y_array
      endif
		
      if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
          parity2.eq.'-' .and. axis2.eq.'y' .and. &
          parity3.eq.'-' .and. axis3.eq.'z' .and. &
          direction.eq.'y') then
         triangular_array => pos_x_neg_y_neg_z_dir_y_triangular_array
         triangular_a_array => pos_x_neg_y_neg_z_dir_y_triangular_y_array
         triangular_b_array => pos_x_neg_y_neg_z_dir_y_triangular_x_array
         triangular_c_array => pos_x_neg_y_neg_z_dir_y_triangular_z_array
      endif
		
      if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
          parity2.eq.'-' .and. axis2.eq.'y' .and. &
          parity3.eq.'-' .and. axis3.eq.'z' .and. &
          direction.eq.'z') then
         triangular_array => pos_x_neg_y_neg_z_dir_z_triangular_array
         triangular_a_array => pos_x_neg_y_neg_z_dir_z_triangular_z_array
         triangular_b_array => pos_x_neg_y_neg_z_dir_z_triangular_y_array
         triangular_c_array => pos_x_neg_y_neg_z_dir_z_triangular_x_array
      endif

      if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
          parity2.eq.'+' .and. axis2.eq.'y' .and. &
          parity3.eq.'+' .and. axis3.eq.'z' .and. &
          direction.eq.'x') then
         triangular_array => neg_x_pos_y_pos_z_dir_x_triangular_array
         triangular_a_array => neg_x_pos_y_pos_z_dir_x_triangular_x_array
         triangular_b_array => neg_x_pos_y_pos_z_dir_x_triangular_z_array
         triangular_c_array => neg_x_pos_y_pos_z_dir_x_triangular_y_array
      endif
		
      if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
          parity2.eq.'+' .and. axis2.eq.'y' .and. &
          parity3.eq.'+' .and. axis3.eq.'z' .and. &
          direction.eq.'y') then
         triangular_array => neg_x_pos_y_pos_z_dir_y_triangular_array
         triangular_a_array => neg_x_pos_y_pos_z_dir_y_triangular_y_array
         triangular_b_array => neg_x_pos_y_pos_z_dir_y_triangular_x_array
         triangular_c_array => neg_x_pos_y_pos_z_dir_y_triangular_z_array
      endif
		
      if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
          parity2.eq.'+' .and. axis2.eq.'y' .and. &
          parity3.eq.'+' .and. axis3.eq.'z' .and. &
          direction.eq.'z') then
         triangular_array => neg_x_pos_y_pos_z_dir_z_triangular_array
         triangular_a_array => neg_x_pos_y_pos_z_dir_z_triangular_z_array
         triangular_b_array => neg_x_pos_y_pos_z_dir_z_triangular_y_array
         triangular_c_array => neg_x_pos_y_pos_z_dir_z_triangular_x_array
      endif
	
      if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
          parity2.eq.'+' .and. axis2.eq.'y' .and. &
          parity3.eq.'-' .and. axis3.eq.'z' .and. &
          direction.eq.'x') then
         triangular_array => neg_x_pos_y_neg_z_dir_x_triangular_array
         triangular_a_array => neg_x_pos_y_neg_z_dir_x_triangular_x_array
         triangular_b_array => neg_x_pos_y_neg_z_dir_x_triangular_z_array
         triangular_c_array => neg_x_pos_y_neg_z_dir_x_triangular_y_array
      endif

      if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
          parity2.eq.'+' .and. axis2.eq.'y' .and. &
          parity3.eq.'-' .and. axis3.eq.'z' .and. &
          direction.eq.'y') then
         triangular_array => neg_x_pos_y_neg_z_dir_y_triangular_array
         triangular_a_array => neg_x_pos_y_neg_z_dir_y_triangular_y_array
         triangular_b_array => neg_x_pos_y_neg_z_dir_y_triangular_x_array
         triangular_c_array => neg_x_pos_y_neg_z_dir_y_triangular_z_array
      endif
		
      if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
          parity2.eq.'+' .and. axis2.eq.'y' .and. &
          parity3.eq.'-' .and. axis3.eq.'z' .and. &
          direction.eq.'z') then
         triangular_array => neg_x_pos_y_neg_z_dir_z_triangular_array
         triangular_a_array => neg_x_pos_y_neg_z_dir_z_triangular_z_array
         triangular_b_array => neg_x_pos_y_neg_z_dir_z_triangular_y_array
         triangular_c_array => neg_x_pos_y_neg_z_dir_z_triangular_x_array
      endif

      if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
          parity2.eq.'-' .and. axis2.eq.'y' .and. &
          parity3.eq.'+' .and. axis3.eq.'z' .and. &
          direction.eq.'x') then
         triangular_array => neg_x_neg_y_pos_z_dir_x_triangular_array
         triangular_a_array => neg_x_neg_y_pos_z_dir_x_triangular_x_array
         triangular_b_array => neg_x_neg_y_pos_z_dir_x_triangular_z_array
         triangular_c_array => neg_x_neg_y_pos_z_dir_x_triangular_y_array
      endif
		
      if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
          parity2.eq.'-' .and. axis2.eq.'y' .and. &
          parity3.eq.'+' .and. axis3.eq.'z' .and. &
          direction.eq.'y') then
         triangular_array => neg_x_neg_y_pos_z_dir_y_triangular_array
         triangular_a_array => neg_x_neg_y_pos_z_dir_y_triangular_y_array
         triangular_b_array => neg_x_neg_y_pos_z_dir_y_triangular_x_array
         triangular_c_array => neg_x_neg_y_pos_z_dir_y_triangular_z_array
      endif
		
      if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
          parity2.eq.'-' .and. axis2.eq.'y' .and. &
          parity3.eq.'+' .and. axis3.eq.'z' .and. &
          direction.eq.'z') then
         triangular_array => neg_x_neg_y_pos_z_dir_z_triangular_array
         triangular_a_array => neg_x_neg_y_pos_z_dir_z_triangular_z_array
         triangular_b_array => neg_x_neg_y_pos_z_dir_z_triangular_y_array
         triangular_c_array => neg_x_neg_y_pos_z_dir_z_triangular_x_array
      endif
	
      if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
          parity2.eq.'-' .and. axis2.eq.'y' .and. &
          parity3.eq.'-' .and. axis3.eq.'z' .and. &
          direction.eq.'x') then
         triangular_array => neg_x_neg_y_neg_z_dir_x_triangular_array
         triangular_a_array => neg_x_neg_y_neg_z_dir_x_triangular_x_array
         triangular_b_array => neg_x_neg_y_neg_z_dir_x_triangular_z_array
         triangular_c_array => neg_x_neg_y_neg_z_dir_x_triangular_y_array
      endif
	
      if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
          parity2.eq.'-' .and. axis2.eq.'y' .and. &
          parity3.eq.'-' .and. axis3.eq.'z' .and. &
          direction.eq.'y') then
         triangular_array => neg_x_neg_y_neg_z_dir_y_triangular_array
         triangular_a_array => neg_x_neg_y_neg_z_dir_y_triangular_y_array
         triangular_b_array => neg_x_neg_y_neg_z_dir_y_triangular_x_array
         triangular_c_array => neg_x_neg_y_neg_z_dir_y_triangular_z_array
      endif
		
      if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
          parity2.eq.'-' .and. axis2.eq.'y' .and. &
          parity3.eq.'-' .and. axis3.eq.'z' .and. &
          direction.eq.'z') then
         triangular_array => neg_x_neg_y_neg_z_dir_z_triangular_array
         triangular_a_array => neg_x_neg_y_neg_z_dir_z_triangular_z_array
         triangular_b_array => neg_x_neg_y_neg_z_dir_z_triangular_y_array
         triangular_c_array => neg_x_neg_y_neg_z_dir_z_triangular_x_array
      endif
	
    endif

    if (direction.eq.axis1 .and. parity1.eq.'+') parity_a = 1.0
    if (direction.eq.axis1 .and. parity1.eq.'-') parity_a = -1.0
    if (direction.eq.axis2 .and. parity2.eq.'+') parity_a = 1.0
    if (direction.eq.axis2 .and. parity2.eq.'-') parity_a = -1.0
    if (direction.eq.axis3 .and. parity3.eq.'+') parity_a = 1.0
    if (direction.eq.axis3 .and. parity3.eq.'-') parity_a = -1.0
	
    if (direction.eq.axis1 .and. parity3.eq.'+') parity_b = 1.0
    if (direction.eq.axis1 .and. parity3.eq.'-') parity_b = -1.0
    if (direction.eq.axis2 .and. parity1.eq.'+') parity_b = 1.0
    if (direction.eq.axis2 .and. parity1.eq.'-') parity_b = -1.0
    if (direction.eq.axis3 .and. parity2.eq.'+') parity_b = 1.0
    if (direction.eq.axis3 .and. parity2.eq.'-') parity_b = -1.0

    if (direction.eq.axis1 .and. parity2.eq.'+') parity_c = 1.0
    if (direction.eq.axis1 .and. parity2.eq.'-') parity_c = -1.0
    if (direction.eq.axis2 .and. parity3.eq.'+') parity_c = 1.0
    if (direction.eq.axis2 .and. parity3.eq.'-') parity_c = -1.0
    if (direction.eq.axis3 .and. parity1.eq.'+') parity_c = 1.0
    if (direction.eq.axis3 .and. parity1.eq.'-') parity_c = -1.0
		 
    count = 1
    pass = .false.
    do while (count.le.length)

      ! triangular_array(count) .ge. 1 means this i-front cells will be checked in this test
      ! for 1 means this cell may be seen i-front cells
      ! and 2 means this cell is seen i-front cells
      ! triangular_array(count) = 0 means this i-front cells no need to be checked anymore
      if (triangular_array(count) .ge. 1) then

        !front test
        a1 = real(triangular_a_array(count)) - parity_a*half 
        b1 = a1*(bb/aa)
        c1 = a1*(cc/aa)
        if (abs(b1-triangular_b_array(count)).le.half .and. &
            abs(c1-triangular_c_array(count)).le.half) then
          pass = .true.
          bb = aa*(real(triangular_b_array(count)) + parity_b*half)/ &
               (real(triangular_a_array(count)) - parity_a*half) + parity_b*0.01
            triangular_array(count) = 2
            exit
        endif	

        !right test
        c1 = real(triangular_c_array(count)) - parity_c*half 
        a1 = c1*(aa/cc)
        b1 = c1*(bb/cc)
        if (abs(a1-triangular_a_array(count)).le.half .and. &
            abs(b1-triangular_b_array(count)).le.half) then
          pass = .true.
          bb = cc*(real(triangular_b_array(count)) + parity_b*half)/ &
               (real(triangular_c_array(count)) - parity_c*half) + parity_b*0.01
            triangular_array(count) = 2
            exit
        endif

        !bottom test	
        b1 = real(triangular_b_array(count)) - parity_b*half 
        a1 = b1*(aa/bb)
        c1 = b1*(cc/bb)
        if (abs(a1-triangular_a_array(count)).le.half .and. &
            abs(c1-triangular_c_array(count)).le.half) then
          a1 = real(triangular_a_array(count)) - parity_a*half 
          c1 = a1*(cc/aa)
          if (abs(c1-triangular_c_array(count)).le.half) then
            pass = .true.
            bb = aa*(real(triangular_b_array(count)) + parity_b*half)/ &
                 (real(triangular_a_array(count)) - parity_a*half) + parity_b*0.01
              triangular_array(count) = 2
              exit
          else
            pass = .true.
            bb = cc*(real(triangular_b_array(count)) + parity_b*half)/ &
                 (real(triangular_c_array(count)) - parity_c*half) + parity_b*0.01
              triangular_array(count) = 2
              exit
          endif

        endif	

      endif

      count = count+1

    enddo

    !no Ifront cells are passed through
    if (pass.eqv..false.) bb = bb + parity_b*1.0 
		
    if (length.ne.0) then

      nullify(triangular_array)
      nullify(triangular_a_array)
      nullify(triangular_b_array)
      nullify(triangular_c_array)

    endif
	
  end subroutine find_non_screened_Ifront_cells

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
  subroutine update_triangular_array(parity1,axis1,parity2,axis2,parity3,axis3,cc,direction,length)
	
    implicit none
	
    character,intent(in) :: parity1
    character,intent(in) :: axis1
    character,intent(in) :: parity2
    character,intent(in) :: axis2
    character,intent(in) :: parity3
    character,intent(in) :: axis3
    real,intent(in) :: cc !this plane will be handled
    character,intent(in) :: direction
    integer,intent(in) :: length
    integer,dimension(:),pointer :: triangular_array
    integer,dimension(:),pointer :: triangular_a_array
    integer,dimension(:),pointer :: triangular_b_array
    integer,dimension(:),pointer :: triangular_c_array
    integer,dimension(:),allocatable :: temporary_array
    integer,dimension(:),allocatable :: temporary_a_array
    integer,dimension(:),allocatable :: temporary_b_array
    integer,dimension(:),allocatable :: temporary_c_array
    integer,dimension(:),allocatable :: non_screened_cell_position
    integer :: count,count2, length1,number_of_non_screened_cell, number_of_deleted_cell
    real :: parity_a, parity_b, parity_c
    integer :: size
		
    if (length.ne.0) then

     if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
         parity2.eq.'+' .and. axis2.eq.'y' .and. &
         parity3.eq.'+' .and. axis3.eq.'z' .and. &
         direction.eq.'x') then
        triangular_array => pos_x_pos_y_pos_z_dir_x_triangular_array
        triangular_a_array => pos_x_pos_y_pos_z_dir_x_triangular_x_array
        triangular_b_array => pos_x_pos_y_pos_z_dir_x_triangular_z_array
        triangular_c_array => pos_x_pos_y_pos_z_dir_x_triangular_y_array
     endif
		
     if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
         parity2.eq.'+' .and. axis2.eq.'y' .and. &
         parity3.eq.'+' .and. axis3.eq.'z' .and. &
         direction.eq.'y') then
        triangular_array => pos_x_pos_y_pos_z_dir_y_triangular_array
        triangular_a_array => pos_x_pos_y_pos_z_dir_y_triangular_y_array
        triangular_b_array => pos_x_pos_y_pos_z_dir_y_triangular_x_array
        triangular_c_array => pos_x_pos_y_pos_z_dir_y_triangular_z_array
     endif
		
     if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
         parity2.eq.'+' .and. axis2.eq.'y' .and. &
         parity3.eq.'+' .and. axis3.eq.'z' .and. &
         direction.eq.'z') then
        triangular_array => pos_x_pos_y_pos_z_dir_z_triangular_array
        triangular_a_array => pos_x_pos_y_pos_z_dir_z_triangular_z_array
        triangular_b_array => pos_x_pos_y_pos_z_dir_z_triangular_y_array
        triangular_c_array => pos_x_pos_y_pos_z_dir_z_triangular_x_array
     endif

     if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
         parity2.eq.'+' .and. axis2.eq.'y' .and. &
         parity3.eq.'-' .and. axis3.eq.'z' .and. &
         direction.eq.'x') then
        triangular_array => pos_x_pos_y_neg_z_dir_x_triangular_array
        triangular_a_array => pos_x_pos_y_neg_z_dir_x_triangular_x_array
        triangular_b_array => pos_x_pos_y_neg_z_dir_x_triangular_z_array
        triangular_c_array => pos_x_pos_y_neg_z_dir_x_triangular_y_array
     endif
		
     if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
         parity2.eq.'+' .and. axis2.eq.'y' .and. &
         parity3.eq.'-' .and. axis3.eq.'z' .and. &
         direction.eq.'y') then
        triangular_array => pos_x_pos_y_neg_z_dir_y_triangular_array
        triangular_a_array => pos_x_pos_y_neg_z_dir_y_triangular_y_array
        triangular_b_array => pos_x_pos_y_neg_z_dir_y_triangular_x_array
        triangular_c_array => pos_x_pos_y_neg_z_dir_y_triangular_z_array
     endif
		
     if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
         parity2.eq.'+' .and. axis2.eq.'y' .and. &
         parity3.eq.'-' .and. axis3.eq.'z' .and. &
         direction.eq.'z') then
        triangular_array => pos_x_pos_y_neg_z_dir_z_triangular_array
        triangular_a_array => pos_x_pos_y_neg_z_dir_z_triangular_z_array
        triangular_b_array => pos_x_pos_y_neg_z_dir_z_triangular_y_array
        triangular_c_array => pos_x_pos_y_neg_z_dir_z_triangular_x_array
     endif

     if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
         parity2.eq.'-' .and. axis2.eq.'y' .and. &
         parity3.eq.'+' .and. axis3.eq.'z' .and. &
         direction.eq.'x') then
        triangular_array => pos_x_neg_y_pos_z_dir_x_triangular_array
        triangular_a_array => pos_x_neg_y_pos_z_dir_x_triangular_x_array
        triangular_b_array => pos_x_neg_y_pos_z_dir_x_triangular_z_array
        triangular_c_array => pos_x_neg_y_pos_z_dir_x_triangular_y_array
     endif
		
     if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
         parity2.eq.'-' .and. axis2.eq.'y' .and. &
         parity3.eq.'+' .and. axis3.eq.'z' .and. &
         direction.eq.'y') then
        triangular_array => pos_x_neg_y_pos_z_dir_y_triangular_array
        triangular_a_array => pos_x_neg_y_pos_z_dir_y_triangular_y_array
        triangular_b_array => pos_x_neg_y_pos_z_dir_y_triangular_x_array
        triangular_c_array => pos_x_neg_y_pos_z_dir_y_triangular_z_array
     endif
		
     if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
         parity2.eq.'-' .and. axis2.eq.'y' .and. &
         parity3.eq.'+' .and. axis3.eq.'z' .and. &
         direction.eq.'z') then
        triangular_array => pos_x_neg_y_pos_z_dir_z_triangular_array
        triangular_a_array => pos_x_neg_y_pos_z_dir_z_triangular_z_array
        triangular_b_array => pos_x_neg_y_pos_z_dir_z_triangular_y_array
        triangular_c_array => pos_x_neg_y_pos_z_dir_z_triangular_x_array
     endif

     if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
         parity2.eq.'-' .and. axis2.eq.'y' .and. &
         parity3.eq.'-' .and. axis3.eq.'z' .and. &
         direction.eq.'x') then
        triangular_array => pos_x_neg_y_neg_z_dir_x_triangular_array
        triangular_a_array => pos_x_neg_y_neg_z_dir_x_triangular_x_array
        triangular_b_array => pos_x_neg_y_neg_z_dir_x_triangular_z_array
        triangular_c_array => pos_x_neg_y_neg_z_dir_x_triangular_y_array
     endif
		
     if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
         parity2.eq.'-' .and. axis2.eq.'y' .and. &
         parity3.eq.'-' .and. axis3.eq.'z' .and. &
         direction.eq.'y') then
        triangular_array => pos_x_neg_y_neg_z_dir_y_triangular_array
        triangular_a_array => pos_x_neg_y_neg_z_dir_y_triangular_y_array
        triangular_b_array => pos_x_neg_y_neg_z_dir_y_triangular_x_array
        triangular_c_array => pos_x_neg_y_neg_z_dir_y_triangular_z_array
     endif
		
     if (parity1.eq.'+' .and. axis1.eq.'x' .and. &
         parity2.eq.'-' .and. axis2.eq.'y' .and. &
         parity3.eq.'-' .and. axis3.eq.'z' .and. &
         direction.eq.'z') then
        triangular_array => pos_x_neg_y_neg_z_dir_z_triangular_array
        triangular_a_array => pos_x_neg_y_neg_z_dir_z_triangular_z_array
        triangular_b_array => pos_x_neg_y_neg_z_dir_z_triangular_y_array
        triangular_c_array => pos_x_neg_y_neg_z_dir_z_triangular_x_array
     endif

     if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
         parity2.eq.'+' .and. axis2.eq.'y' .and. &
         parity3.eq.'+' .and. axis3.eq.'z' .and. &
         direction.eq.'x') then
        triangular_array => neg_x_pos_y_pos_z_dir_x_triangular_array
        triangular_a_array => neg_x_pos_y_pos_z_dir_x_triangular_x_array
        triangular_b_array => neg_x_pos_y_pos_z_dir_x_triangular_z_array
        triangular_c_array => neg_x_pos_y_pos_z_dir_x_triangular_y_array
     endif
		
     if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
         parity2.eq.'+' .and. axis2.eq.'y' .and. &
         parity3.eq.'+' .and. axis3.eq.'z' .and. &
         direction.eq.'y') then
        triangular_array => neg_x_pos_y_pos_z_dir_y_triangular_array
        triangular_a_array => neg_x_pos_y_pos_z_dir_y_triangular_y_array
        triangular_b_array => neg_x_pos_y_pos_z_dir_y_triangular_x_array
        triangular_c_array => neg_x_pos_y_pos_z_dir_y_triangular_z_array
     endif
		
     if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
         parity2.eq.'+' .and. axis2.eq.'y' .and. &
         parity3.eq.'+' .and. axis3.eq.'z' .and. &
         direction.eq.'z') then
        triangular_array => neg_x_pos_y_pos_z_dir_z_triangular_array
        triangular_a_array => neg_x_pos_y_pos_z_dir_z_triangular_z_array
        triangular_b_array => neg_x_pos_y_pos_z_dir_z_triangular_y_array
        triangular_c_array => neg_x_pos_y_pos_z_dir_z_triangular_x_array
     endif

     if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
         parity2.eq.'+' .and. axis2.eq.'y' .and. &
         parity3.eq.'-' .and. axis3.eq.'z' .and. &
         direction.eq.'x') then
        triangular_array => neg_x_pos_y_neg_z_dir_x_triangular_array
        triangular_a_array => neg_x_pos_y_neg_z_dir_x_triangular_x_array
        triangular_b_array => neg_x_pos_y_neg_z_dir_x_triangular_z_array
        triangular_c_array => neg_x_pos_y_neg_z_dir_x_triangular_y_array
     endif
		
     if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
         parity2.eq.'+' .and. axis2.eq.'y' .and. &
         parity3.eq.'-' .and. axis3.eq.'z' .and. &
         direction.eq.'y') then
        triangular_array => neg_x_pos_y_neg_z_dir_y_triangular_array
        triangular_a_array => neg_x_pos_y_neg_z_dir_y_triangular_y_array
        triangular_b_array => neg_x_pos_y_neg_z_dir_y_triangular_x_array
        triangular_c_array => neg_x_pos_y_neg_z_dir_y_triangular_z_array
     endif
		
     if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
         parity2.eq.'+' .and. axis2.eq.'y' .and. &
         parity3.eq.'-' .and. axis3.eq.'z' .and. &
         direction.eq.'z') then
        triangular_array => neg_x_pos_y_neg_z_dir_z_triangular_array
        triangular_a_array => neg_x_pos_y_neg_z_dir_z_triangular_z_array
        triangular_b_array => neg_x_pos_y_neg_z_dir_z_triangular_y_array
        triangular_c_array => neg_x_pos_y_neg_z_dir_z_triangular_x_array
     endif

    if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
        parity2.eq.'-' .and. axis2.eq.'y' .and. &
        parity3.eq.'+' .and. axis3.eq.'z' .and. &
        direction.eq.'x') then
       triangular_array => neg_x_neg_y_pos_z_dir_x_triangular_array
       triangular_a_array => neg_x_neg_y_pos_z_dir_x_triangular_x_array
       triangular_b_array => neg_x_neg_y_pos_z_dir_x_triangular_z_array
       triangular_c_array => neg_x_neg_y_pos_z_dir_x_triangular_y_array
    endif
		
    if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
        parity2.eq.'-' .and. axis2.eq.'y' .and. &
        parity3.eq.'+' .and. axis3.eq.'z' .and. &
        direction.eq.'y') then
       triangular_array => neg_x_neg_y_pos_z_dir_y_triangular_array
       triangular_a_array => neg_x_neg_y_pos_z_dir_y_triangular_y_array
       triangular_b_array => neg_x_neg_y_pos_z_dir_y_triangular_x_array
       triangular_c_array => neg_x_neg_y_pos_z_dir_y_triangular_z_array
    endif
		
    if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
        parity2.eq.'-' .and. axis2.eq.'y' .and. &
        parity3.eq.'+' .and. axis3.eq.'z' .and. &
        direction.eq.'z') then
       triangular_array => neg_x_neg_y_pos_z_dir_z_triangular_array
       triangular_a_array => neg_x_neg_y_pos_z_dir_z_triangular_z_array
       triangular_b_array => neg_x_neg_y_pos_z_dir_z_triangular_y_array
       triangular_c_array => neg_x_neg_y_pos_z_dir_z_triangular_x_array
    endif

    if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
        parity2.eq.'-' .and. axis2.eq.'y' .and. &
        parity3.eq.'-' .and. axis3.eq.'z' .and. &
        direction.eq.'x') then
       triangular_array => neg_x_neg_y_neg_z_dir_x_triangular_array
       triangular_a_array => neg_x_neg_y_neg_z_dir_x_triangular_x_array
       triangular_b_array => neg_x_neg_y_neg_z_dir_x_triangular_z_array
       triangular_c_array => neg_x_neg_y_neg_z_dir_x_triangular_y_array
    endif
		
    if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
        parity2.eq.'-' .and. axis2.eq.'y' .and. &
        parity3.eq.'-' .and. axis3.eq.'z' .and. &
        direction.eq.'y') then
       triangular_array => neg_x_neg_y_neg_z_dir_y_triangular_array
       triangular_a_array => neg_x_neg_y_neg_z_dir_y_triangular_y_array
       triangular_b_array => neg_x_neg_y_neg_z_dir_y_triangular_x_array
       triangular_c_array => neg_x_neg_y_neg_z_dir_y_triangular_z_array
    endif
		
    if (parity1.eq.'-' .and. axis1.eq.'x' .and. &
        parity2.eq.'-' .and. axis2.eq.'y' .and. &
        parity3.eq.'-' .and. axis3.eq.'z' .and. &
        direction.eq.'z') then
       triangular_array => neg_x_neg_y_neg_z_dir_z_triangular_array
       triangular_a_array => neg_x_neg_y_neg_z_dir_z_triangular_z_array
       triangular_b_array => neg_x_neg_y_neg_z_dir_z_triangular_y_array
       triangular_c_array => neg_x_neg_y_neg_z_dir_z_triangular_x_array
    endif
		
    if (direction.eq.axis1 .and. parity1.eq.'+') parity_a = 1.0
    if (direction.eq.axis1 .and. parity1.eq.'-') parity_a = -1.0
    if (direction.eq.axis2 .and. parity2.eq.'+') parity_a = 1.0
    if (direction.eq.axis2 .and. parity2.eq.'-') parity_a = -1.0
    if (direction.eq.axis3 .and. parity3.eq.'+') parity_a = 1.0
    if (direction.eq.axis3 .and. parity3.eq.'-') parity_a = -1.0
		
    if (direction.eq.axis1 .and. parity3.eq.'+') parity_b = 1.0
    if (direction.eq.axis1 .and. parity3.eq.'-') parity_b = -1.0
    if (direction.eq.axis2 .and. parity1.eq.'+') parity_b = 1.0
    if (direction.eq.axis2 .and. parity1.eq.'-') parity_b = -1.0
    if (direction.eq.axis3 .and. parity2.eq.'+') parity_b = 1.0
    if (direction.eq.axis3 .and. parity2.eq.'-') parity_b = -1.0 
		
    if (direction.eq.axis1 .and. parity2.eq.'+') parity_c = 1.0
    if (direction.eq.axis1 .and. parity2.eq.'-') parity_c = -1.0
    if (direction.eq.axis2 .and. parity3.eq.'+') parity_c = 1.0
    if (direction.eq.axis2 .and. parity3.eq.'-') parity_c = -1.0
    if (direction.eq.axis3 .and. parity1.eq.'+') parity_c = 1.0
    if (direction.eq.axis3 .and. parity1.eq.'-') parity_c = -1.0
		
    allocate(non_screened_cell_position(1:size(triangular_array)))

    non_screened_cell_position = 0		
    number_of_non_screened_cell = 0
	   	
    ! this loop find outs where all the seen I-front cells (number 2 cells) are for slices level .le. abs(cc)
    do count = 1,size(triangular_array)
      if (abs(triangular_c_array(count)).le.abs(cc) .and. triangular_array(count).eq.2) then
        number_of_non_screened_cell = number_of_non_screened_cell + 1
	non_screened_cell_position(number_of_non_screened_cell) = count
      endif
    enddo

    if (number_of_non_screened_cell.ne.0) then
      do count = 1,number_of_non_screened_cell-1
        if ((abs(triangular_b_array(non_screened_cell_position(count+1))).gt. &
          abs(triangular_b_array(non_screened_cell_position(count)))) .and. &
          (triangular_c_array(non_screened_cell_position(count+1)).eq. &
          triangular_c_array(non_screened_cell_position(count)))) then
          do count2 = non_screened_cell_position(count)+1,non_screened_cell_position(count+1)-1
              if (triangular_array(count2).ne.0 .and. &
                 parity_a*triangular_a_array(count2).ge. &
                 parity_a*triangular_a_array(non_screened_cell_position(count+1))) then
                triangular_array(count2) = 0
              endif      
          enddo
        endif
      enddo
    endif

endif
	
  end subroutine update_triangular_array
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_non_screened_Ifront_array
	
    implicit none
	
    integer :: count


!$OMP SECTIONS PRIVATE(count)

!$OMP SECTION
	
    do count = 1,size(pos_x_pos_y_pos_z_dir_x_triangular_array)
      if (pos_x_pos_y_pos_z_dir_x_triangular_array(count).eq.2) then
        non_screened_Ifront_array(pos_x_pos_y_pos_z_dir_x_triangular_x_array(count), &
                                  pos_x_pos_y_pos_z_dir_x_triangular_y_array(count), &
                                  pos_x_pos_y_pos_z_dir_x_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(pos_x_pos_y_pos_z_dir_y_triangular_array)
      if (pos_x_pos_y_pos_z_dir_y_triangular_array(count).eq.2) then
        non_screened_Ifront_array(pos_x_pos_y_pos_z_dir_y_triangular_x_array(count), &
                                  pos_x_pos_y_pos_z_dir_y_triangular_y_array(count), &
                                  pos_x_pos_y_pos_z_dir_y_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(pos_x_pos_y_pos_z_dir_z_triangular_array)
      if (pos_x_pos_y_pos_z_dir_z_triangular_array(count).eq.2) then
        non_screened_Ifront_array(pos_x_pos_y_pos_z_dir_z_triangular_x_array(count), &
                                  pos_x_pos_y_pos_z_dir_z_triangular_y_array(count), &
                                  pos_x_pos_y_pos_z_dir_z_triangular_z_array(count)) = 1
      endif
    enddo

!$OMP SECTION	

    do count = 1,size(pos_x_pos_y_neg_z_dir_x_triangular_array)
      if (pos_x_pos_y_neg_z_dir_x_triangular_array(count).eq.2) then
        non_screened_Ifront_array(pos_x_pos_y_neg_z_dir_x_triangular_x_array(count), &
                                  pos_x_pos_y_neg_z_dir_x_triangular_y_array(count), &
                                  pos_x_pos_y_neg_z_dir_x_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(pos_x_pos_y_neg_z_dir_y_triangular_array)
      if (pos_x_pos_y_neg_z_dir_y_triangular_array(count).eq.2) then
        non_screened_Ifront_array(pos_x_pos_y_neg_z_dir_y_triangular_x_array(count), &
                                  pos_x_pos_y_neg_z_dir_y_triangular_y_array(count), &
                                  pos_x_pos_y_neg_z_dir_y_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(pos_x_pos_y_neg_z_dir_z_triangular_array)
      if (pos_x_pos_y_neg_z_dir_z_triangular_array(count).eq.2) then
        non_screened_Ifront_array(pos_x_pos_y_neg_z_dir_z_triangular_x_array(count), &
                                  pos_x_pos_y_neg_z_dir_z_triangular_y_array(count), &
                                  pos_x_pos_y_neg_z_dir_z_triangular_z_array(count)) = 1
      endif
    enddo

!$OMP SECTION	

    do count = 1,size(pos_x_neg_y_pos_z_dir_x_triangular_array)
      if (pos_x_neg_y_pos_z_dir_x_triangular_array(count).eq.2) then
        non_screened_Ifront_array(pos_x_neg_y_pos_z_dir_x_triangular_x_array(count), &
                                  pos_x_neg_y_pos_z_dir_x_triangular_y_array(count), &
                                  pos_x_neg_y_pos_z_dir_x_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(pos_x_neg_y_pos_z_dir_y_triangular_array)
      if (pos_x_neg_y_pos_z_dir_y_triangular_array(count).eq.2) then
        non_screened_Ifront_array(pos_x_neg_y_pos_z_dir_y_triangular_x_array(count), &
                                  pos_x_neg_y_pos_z_dir_y_triangular_y_array(count), &
                                  pos_x_neg_y_pos_z_dir_y_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(pos_x_neg_y_pos_z_dir_z_triangular_array)
      if (pos_x_neg_y_pos_z_dir_z_triangular_array(count).eq.2) then
        non_screened_Ifront_array(pos_x_neg_y_pos_z_dir_z_triangular_x_array(count), &
                                  pos_x_neg_y_pos_z_dir_z_triangular_y_array(count), &
                                  pos_x_neg_y_pos_z_dir_z_triangular_z_array(count)) = 1
      endif
    enddo

!$OMP SECTION
	
    do count = 1,size(pos_x_neg_y_neg_z_dir_x_triangular_array)
      if (pos_x_neg_y_neg_z_dir_x_triangular_array(count).eq.2) then
        non_screened_Ifront_array(pos_x_neg_y_neg_z_dir_x_triangular_x_array(count), &
                                  pos_x_neg_y_neg_z_dir_x_triangular_y_array(count), &
                                  pos_x_neg_y_neg_z_dir_x_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(pos_x_neg_y_neg_z_dir_y_triangular_array)
      if (pos_x_neg_y_neg_z_dir_y_triangular_array(count).eq.2) then
        non_screened_Ifront_array(pos_x_neg_y_neg_z_dir_y_triangular_x_array(count), &
                                  pos_x_neg_y_neg_z_dir_y_triangular_y_array(count), &
                                  pos_x_neg_y_neg_z_dir_y_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(pos_x_neg_y_neg_z_dir_z_triangular_array)
      if (pos_x_neg_y_neg_z_dir_z_triangular_array(count).eq.2) then
        non_screened_Ifront_array(pos_x_neg_y_neg_z_dir_z_triangular_x_array(count), &
                                  pos_x_neg_y_neg_z_dir_z_triangular_y_array(count), &
                                  pos_x_neg_y_neg_z_dir_z_triangular_z_array(count)) = 1
      endif
    enddo

!$OMP SECTION
	
    do count = 1,size(neg_x_pos_y_pos_z_dir_x_triangular_array)
      if (neg_x_pos_y_pos_z_dir_x_triangular_array(count).eq.2) then
        non_screened_Ifront_array(neg_x_pos_y_pos_z_dir_x_triangular_x_array(count), &
                                  neg_x_pos_y_pos_z_dir_x_triangular_y_array(count), &
                                  neg_x_pos_y_pos_z_dir_x_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(neg_x_pos_y_pos_z_dir_y_triangular_array)
      if (neg_x_pos_y_pos_z_dir_y_triangular_array(count).eq.2) then
        non_screened_Ifront_array(neg_x_pos_y_pos_z_dir_y_triangular_x_array(count), &
                                  neg_x_pos_y_pos_z_dir_y_triangular_y_array(count), &
                                  neg_x_pos_y_pos_z_dir_y_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(neg_x_pos_y_pos_z_dir_z_triangular_array)
      if (neg_x_pos_y_pos_z_dir_z_triangular_array(count).eq.2) then
        non_screened_Ifront_array(neg_x_pos_y_pos_z_dir_z_triangular_x_array(count), &
                                  neg_x_pos_y_pos_z_dir_z_triangular_y_array(count), &
                                  neg_x_pos_y_pos_z_dir_z_triangular_z_array(count)) = 1
      endif
    enddo

!$OMP SECTION
	
    do count = 1,size(neg_x_pos_y_neg_z_dir_x_triangular_array)
      if (neg_x_pos_y_neg_z_dir_x_triangular_array(count).eq.2) then
        non_screened_Ifront_array(neg_x_pos_y_neg_z_dir_x_triangular_x_array(count), &
                                  neg_x_pos_y_neg_z_dir_x_triangular_y_array(count), &
                                  neg_x_pos_y_neg_z_dir_x_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(neg_x_pos_y_neg_z_dir_y_triangular_array)
      if (neg_x_pos_y_neg_z_dir_y_triangular_array(count).eq.2) then
        non_screened_Ifront_array(neg_x_pos_y_neg_z_dir_y_triangular_x_array(count), &
                                  neg_x_pos_y_neg_z_dir_y_triangular_y_array(count), &
                                  neg_x_pos_y_neg_z_dir_y_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(neg_x_pos_y_neg_z_dir_z_triangular_array)
      if (neg_x_pos_y_neg_z_dir_z_triangular_array(count).eq.2) then
        non_screened_Ifront_array(neg_x_pos_y_neg_z_dir_z_triangular_x_array(count), &
                                  neg_x_pos_y_neg_z_dir_z_triangular_y_array(count), &
                                  neg_x_pos_y_neg_z_dir_z_triangular_z_array(count)) = 1
      endif
    enddo

!$OMP SECTION
	
    do count = 1,size(neg_x_neg_y_pos_z_dir_x_triangular_array)
      if (neg_x_neg_y_pos_z_dir_x_triangular_array(count).eq.2) then
        non_screened_Ifront_array(neg_x_neg_y_pos_z_dir_x_triangular_x_array(count), &
                                  neg_x_neg_y_pos_z_dir_x_triangular_y_array(count), &
                                  neg_x_neg_y_pos_z_dir_x_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(neg_x_neg_y_pos_z_dir_y_triangular_array)
      if (neg_x_neg_y_pos_z_dir_y_triangular_array(count).eq.2) then
        non_screened_Ifront_array(neg_x_neg_y_pos_z_dir_y_triangular_x_array(count), &
                                    neg_x_neg_y_pos_z_dir_y_triangular_y_array(count), &
                                  neg_x_neg_y_pos_z_dir_y_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(neg_x_neg_y_pos_z_dir_z_triangular_array)
      if (neg_x_neg_y_pos_z_dir_z_triangular_array(count).eq.2) then
        non_screened_Ifront_array(neg_x_neg_y_pos_z_dir_z_triangular_x_array(count), &
                                  neg_x_neg_y_pos_z_dir_z_triangular_y_array(count), &
                                  neg_x_neg_y_pos_z_dir_z_triangular_z_array(count)) = 1
      endif
    enddo

!$OMP SECTION	

    do count = 1,size(neg_x_neg_y_neg_z_dir_x_triangular_array)
      if (neg_x_neg_y_neg_z_dir_x_triangular_array(count).eq.2) then
        non_screened_Ifront_array(neg_x_neg_y_neg_z_dir_x_triangular_x_array(count), &
                                  neg_x_neg_y_neg_z_dir_x_triangular_y_array(count), &
                                  neg_x_neg_y_neg_z_dir_x_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(neg_x_neg_y_neg_z_dir_y_triangular_array)
      if (neg_x_neg_y_neg_z_dir_y_triangular_array(count).eq.2) then
        non_screened_Ifront_array(neg_x_neg_y_neg_z_dir_y_triangular_x_array(count), &
                                  neg_x_neg_y_neg_z_dir_y_triangular_y_array(count), &
                                  neg_x_neg_y_neg_z_dir_y_triangular_z_array(count)) = 1
      endif
    enddo

    do count = 1,size(neg_x_neg_y_neg_z_dir_z_triangular_array)
      if (neg_x_neg_y_neg_z_dir_z_triangular_array(count).eq.2) then
        non_screened_Ifront_array(neg_x_neg_y_neg_z_dir_z_triangular_x_array(count), &
                                  neg_x_neg_y_neg_z_dir_z_triangular_y_array(count), &
                                  neg_x_neg_y_neg_z_dir_z_triangular_z_array(count)) = 1
      endif
    enddo

!$OMP END SECTIONS
	
  end subroutine update_non_screened_Ifront_array
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module AT_ray
