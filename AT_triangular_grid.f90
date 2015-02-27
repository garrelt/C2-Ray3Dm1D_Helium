module AT_triangular_grid

  use my_mpi
  use parameter, only: length_pos_x, length_neg_x, &
                       length_pos_y, length_neg_y, &
                       length_pos_z, length_neg_z
  use array, only: pos_x_pos_y_pos_z_octant, pos_x_pos_y_neg_z_octant, &
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
                   neg_x_neg_y_neg_z_dir_z_triangular_z_array, &
                   source_Ifront_array

  implicit none

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AT_triangular_grid_initialization()

    implicit none

    integer :: pos_x_pos_y_pos_z_dir_x_length
    integer :: pos_x_pos_y_pos_z_dir_y_length
    integer :: pos_x_pos_y_pos_z_dir_z_length
    integer :: pos_x_pos_y_neg_z_dir_x_length
    integer :: pos_x_pos_y_neg_z_dir_y_length
    integer :: pos_x_pos_y_neg_z_dir_z_length	
    integer :: pos_x_neg_y_pos_z_dir_x_length
    integer :: pos_x_neg_y_pos_z_dir_y_length
    integer :: pos_x_neg_y_pos_z_dir_z_length
    integer :: pos_x_neg_y_neg_z_dir_x_length
    integer :: pos_x_neg_y_neg_z_dir_y_length
    integer :: pos_x_neg_y_neg_z_dir_z_length
    integer :: neg_x_pos_y_pos_z_dir_x_length
    integer :: neg_x_pos_y_pos_z_dir_y_length
    integer :: neg_x_pos_y_pos_z_dir_z_length
    integer :: neg_x_pos_y_neg_z_dir_x_length
    integer :: neg_x_pos_y_neg_z_dir_y_length
    integer :: neg_x_pos_y_neg_z_dir_z_length
    integer :: neg_x_neg_y_pos_z_dir_x_length
    integer :: neg_x_neg_y_pos_z_dir_y_length
    integer :: neg_x_neg_y_pos_z_dir_z_length
    integer :: neg_x_neg_y_neg_z_dir_x_length
    integer :: neg_x_neg_y_neg_z_dir_y_length
    integer :: neg_x_neg_y_neg_z_dir_z_length	
    integer :: i, j ,k
    integer :: count
		
    pos_x_pos_y_pos_z_dir_x_length = 0
    pos_x_pos_y_pos_z_dir_y_length = 0
    pos_x_pos_y_pos_z_dir_z_length = 0
    pos_x_pos_y_neg_z_dir_x_length = 0
    pos_x_pos_y_neg_z_dir_y_length = 0
    pos_x_pos_y_neg_z_dir_z_length = 0		
    pos_x_neg_y_pos_z_dir_x_length = 0
    pos_x_neg_y_pos_z_dir_y_length = 0
    pos_x_neg_y_pos_z_dir_z_length = 0		
    pos_x_neg_y_neg_z_dir_x_length = 0
    pos_x_neg_y_neg_z_dir_y_length = 0
    pos_x_neg_y_neg_z_dir_z_length = 0
    neg_x_pos_y_pos_z_dir_x_length = 0
    neg_x_pos_y_pos_z_dir_y_length = 0
    neg_x_pos_y_pos_z_dir_z_length = 0		
    neg_x_pos_y_neg_z_dir_x_length = 0
    neg_x_pos_y_neg_z_dir_y_length = 0
    neg_x_pos_y_neg_z_dir_z_length = 0		
    neg_x_neg_y_pos_z_dir_x_length = 0
    neg_x_neg_y_pos_z_dir_y_length = 0
    neg_x_neg_y_pos_z_dir_z_length = 0		
    neg_x_neg_y_neg_z_dir_x_length = 0
    neg_x_neg_y_neg_z_dir_y_length = 0
    neg_x_neg_y_neg_z_dir_z_length = 0

!$OMP PARALLEL PRIVATE(i,j,k,count)

!$OMP SECTIONS 

!$OMP SECTION	

    do j = 0, length_pos_y
      do k = 0, length_pos_z
        do i = k, length_pos_x
          if (i.ge.j .or. k.ge.j) then
            if (pos_x_pos_y_pos_z_octant(i,j,k).eq.1) then
              pos_x_pos_y_pos_z_dir_x_length = pos_x_pos_y_pos_z_dir_x_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(pos_x_pos_y_pos_z_dir_x_triangular_array(1:pos_x_pos_y_pos_z_dir_x_length))
    allocate(pos_x_pos_y_pos_z_dir_x_triangular_x_array(1:pos_x_pos_y_pos_z_dir_x_length))
    allocate(pos_x_pos_y_pos_z_dir_x_triangular_y_array(1:pos_x_pos_y_pos_z_dir_x_length))
    allocate(pos_x_pos_y_pos_z_dir_x_triangular_z_array(1:pos_x_pos_y_pos_z_dir_x_length))		
    count = 1
    do j = 0, length_pos_y
      do k = 0, length_pos_z
        do i = k, length_pos_x
          if (i.ge.j .or. k.ge.j) then
            if (pos_x_pos_y_pos_z_octant(i,j,k).eq.1) then
              pos_x_pos_y_pos_z_dir_x_triangular_array(count) = 1
              pos_x_pos_y_pos_z_dir_x_triangular_x_array(count) = i
              pos_x_pos_y_pos_z_dir_x_triangular_y_array(count) = j
              pos_x_pos_y_pos_z_dir_x_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo

    do k = 0, length_pos_z
      do i = 0, length_pos_x
        do j = i, length_pos_y
          if (j.ge.k .or. i.ge.k) then
            if (pos_x_pos_y_pos_z_octant(i,j,k).eq.1) then
              pos_x_pos_y_pos_z_dir_y_length = pos_x_pos_y_pos_z_dir_y_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(pos_x_pos_y_pos_z_dir_y_triangular_array(1:pos_x_pos_y_pos_z_dir_y_length))
    allocate(pos_x_pos_y_pos_z_dir_y_triangular_x_array(1:pos_x_pos_y_pos_z_dir_y_length))
    allocate(pos_x_pos_y_pos_z_dir_y_triangular_y_array(1:pos_x_pos_y_pos_z_dir_y_length))
    allocate(pos_x_pos_y_pos_z_dir_y_triangular_z_array(1:pos_x_pos_y_pos_z_dir_y_length))	
    count = 1
    do k = 0, length_pos_z
      do i = 0, length_pos_x
        do j = i, length_pos_y
          if (j.ge.k .or. i.ge.k) then
            if (pos_x_pos_y_pos_z_octant(i,j,k).eq.1) then
              pos_x_pos_y_pos_z_dir_y_triangular_array(count) = 1
              pos_x_pos_y_pos_z_dir_y_triangular_x_array(count) = i
              pos_x_pos_y_pos_z_dir_y_triangular_y_array(count) = j
              pos_x_pos_y_pos_z_dir_y_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
	
    do i = 0, length_pos_x
      do j = 0, length_pos_y
        do k = j, length_pos_z
          if (k.ge.i .or. j.ge.i) then
            if (pos_x_pos_y_pos_z_octant(i,j,k).eq.1) then
              pos_x_pos_y_pos_z_dir_z_length = pos_x_pos_y_pos_z_dir_z_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(pos_x_pos_y_pos_z_dir_z_triangular_array(1:pos_x_pos_y_pos_z_dir_z_length))
    allocate(pos_x_pos_y_pos_z_dir_z_triangular_x_array(1:pos_x_pos_y_pos_z_dir_z_length))
    allocate(pos_x_pos_y_pos_z_dir_z_triangular_y_array(1:pos_x_pos_y_pos_z_dir_z_length))
    allocate(pos_x_pos_y_pos_z_dir_z_triangular_z_array(1:pos_x_pos_y_pos_z_dir_z_length))	
    count = 1
    do i = 0, length_pos_x
      do j = 0, length_pos_y
        do k = j, length_pos_z
          if (k.ge.i .or. j.ge.i) then
            if (pos_x_pos_y_pos_z_octant(i,j,k).eq.1) then
              pos_x_pos_y_pos_z_dir_z_triangular_array(count) = 1
              pos_x_pos_y_pos_z_dir_z_triangular_x_array(count) = i
              pos_x_pos_y_pos_z_dir_z_triangular_y_array(count) = j
              pos_x_pos_y_pos_z_dir_z_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
!$OMP SECTION		
    do j = 0, length_pos_y
      do k = 0, length_neg_z, -1
        do i = -k, length_pos_x
          if (i.ge.j .or. -k.ge.j) then
            if (pos_x_pos_y_neg_z_octant(i,j,k).eq.1) then
              pos_x_pos_y_neg_z_dir_x_length = pos_x_pos_y_neg_z_dir_x_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(pos_x_pos_y_neg_z_dir_x_triangular_array(1:pos_x_pos_y_neg_z_dir_x_length))
    allocate(pos_x_pos_y_neg_z_dir_x_triangular_x_array(1:pos_x_pos_y_neg_z_dir_x_length))
    allocate(pos_x_pos_y_neg_z_dir_x_triangular_y_array(1:pos_x_pos_y_neg_z_dir_x_length))
    allocate(pos_x_pos_y_neg_z_dir_x_triangular_z_array(1:pos_x_pos_y_neg_z_dir_x_length))		
    count = 1
    do j = 0, length_pos_y
      do k = 0, length_neg_z, -1
        do i = -k, length_pos_x
          if (i.ge.j .or. -k.ge.j) then
            if (pos_x_pos_y_neg_z_octant(i,j,k).eq.1) then
              pos_x_pos_y_neg_z_dir_x_triangular_array(count) = 1
              pos_x_pos_y_neg_z_dir_x_triangular_x_array(count) = i
              pos_x_pos_y_neg_z_dir_x_triangular_y_array(count) = j
              pos_x_pos_y_neg_z_dir_x_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
	
    do k = 0, length_neg_z, -1
      do i = 0, length_pos_x
        do j = i, length_pos_y
          if (j.ge.-k .or. i.ge.-k) then
            if (pos_x_pos_y_neg_z_octant(i,j,k).eq.1) then
              pos_x_pos_y_neg_z_dir_y_length = pos_x_pos_y_neg_z_dir_y_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(pos_x_pos_y_neg_z_dir_y_triangular_array(1:pos_x_pos_y_neg_z_dir_y_length))
    allocate(pos_x_pos_y_neg_z_dir_y_triangular_x_array(1:pos_x_pos_y_neg_z_dir_y_length))
    allocate(pos_x_pos_y_neg_z_dir_y_triangular_y_array(1:pos_x_pos_y_neg_z_dir_y_length))
    allocate(pos_x_pos_y_neg_z_dir_y_triangular_z_array(1:pos_x_pos_y_neg_z_dir_y_length))	
    count = 1
    do k = 0, length_neg_z, -1
      do i = 0, length_pos_x
        do j = i, length_pos_y
          if (j.ge.-k .or. i.ge.-k) then
            if (pos_x_pos_y_neg_z_octant(i,j,k).eq.1) then
              pos_x_pos_y_neg_z_dir_y_triangular_array(count) = 1
              pos_x_pos_y_neg_z_dir_y_triangular_x_array(count) = i
              pos_x_pos_y_neg_z_dir_y_triangular_y_array(count) = j
              pos_x_pos_y_neg_z_dir_y_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
	
    do i = 0, length_pos_x
      do j = 0, length_pos_y
        do k = -j, length_neg_z, -1
          if (-k.ge.i .or. j.ge.i) then
            if (pos_x_pos_y_neg_z_octant(i,j,k).eq.1) then
              pos_x_pos_y_neg_z_dir_z_length = pos_x_pos_y_neg_z_dir_z_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(pos_x_pos_y_neg_z_dir_z_triangular_array(1:pos_x_pos_y_neg_z_dir_z_length))
    allocate(pos_x_pos_y_neg_z_dir_z_triangular_x_array(1:pos_x_pos_y_neg_z_dir_z_length))
    allocate(pos_x_pos_y_neg_z_dir_z_triangular_y_array(1:pos_x_pos_y_neg_z_dir_z_length))
    allocate(pos_x_pos_y_neg_z_dir_z_triangular_z_array(1:pos_x_pos_y_neg_z_dir_z_length))	
    count = 1
    do i = 0, length_pos_x
      do j = 0, length_pos_y
        do k = -j, length_neg_z, -1
          if (-k.ge.i .or. j.ge.i) then
            if (pos_x_pos_y_neg_z_octant(i,j,k).eq.1) then
              pos_x_pos_y_neg_z_dir_z_triangular_array(count) = 1
              pos_x_pos_y_neg_z_dir_z_triangular_x_array(count) = i
              pos_x_pos_y_neg_z_dir_z_triangular_y_array(count) = j
              pos_x_pos_y_neg_z_dir_z_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
!$OMP SECTION
    do j = 0, length_neg_y, -1
      do k = 0, length_pos_z
        do i = k, length_pos_x
          if (i.ge.-j .or. k.ge.-j) then
            if (pos_x_neg_y_pos_z_octant(i,j,k).eq.1) then
              pos_x_neg_y_pos_z_dir_x_length = pos_x_neg_y_pos_z_dir_x_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(pos_x_neg_y_pos_z_dir_x_triangular_array(1:pos_x_neg_y_pos_z_dir_x_length))
    allocate(pos_x_neg_y_pos_z_dir_x_triangular_x_array(1:pos_x_neg_y_pos_z_dir_x_length))
    allocate(pos_x_neg_y_pos_z_dir_x_triangular_y_array(1:pos_x_neg_y_pos_z_dir_x_length))
    allocate(pos_x_neg_y_pos_z_dir_x_triangular_z_array(1:pos_x_neg_y_pos_z_dir_x_length))		
    count = 1
    do j = 0, length_neg_y, -1
      do k = 0, length_pos_z
        do i = k, length_pos_x
          if (i.ge.-j .or. k.ge.-j) then
            if (pos_x_neg_y_pos_z_octant(i,j,k).eq.1) then
              pos_x_neg_y_pos_z_dir_x_triangular_array(count) = 1
              pos_x_neg_y_pos_z_dir_x_triangular_x_array(count) = i
              pos_x_neg_y_pos_z_dir_x_triangular_y_array(count) = j
              pos_x_neg_y_pos_z_dir_x_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
	
    do k = 0, length_pos_z
      do i = 0, length_pos_x
        do j = -i, length_neg_y, -1
          if (-j.ge.k .or. i.ge.k) then
            if (pos_x_neg_y_pos_z_octant(i,j,k).eq.1) then
              pos_x_neg_y_pos_z_dir_y_length = pos_x_neg_y_pos_z_dir_y_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(pos_x_neg_y_pos_z_dir_y_triangular_array(1:pos_x_neg_y_pos_z_dir_y_length))
    allocate(pos_x_neg_y_pos_z_dir_y_triangular_x_array(1:pos_x_neg_y_pos_z_dir_y_length))
    allocate(pos_x_neg_y_pos_z_dir_y_triangular_y_array(1:pos_x_neg_y_pos_z_dir_y_length))
    allocate(pos_x_neg_y_pos_z_dir_y_triangular_z_array(1:pos_x_neg_y_pos_z_dir_y_length))	
    count = 1
    do k = 0, length_pos_z
      do i = 0, length_pos_x
        do j = -i, length_neg_y, -1
          if (-j.ge.k .or. i.ge.k) then
            if (pos_x_neg_y_pos_z_octant(i,j,k).eq.1) then
              pos_x_neg_y_pos_z_dir_y_triangular_array(count) = 1
              pos_x_neg_y_pos_z_dir_y_triangular_x_array(count) = i
              pos_x_neg_y_pos_z_dir_y_triangular_y_array(count) = j
              pos_x_neg_y_pos_z_dir_y_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
	
    do i = 0, length_pos_x
      do j = 0, length_neg_y, -1
        do k = -j, length_pos_z
          if (k.ge.i .or. -j.ge.i) then
            if (pos_x_neg_y_pos_z_octant(i,j,k).eq.1) then
              pos_x_neg_y_pos_z_dir_z_length = pos_x_neg_y_pos_z_dir_z_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(pos_x_neg_y_pos_z_dir_z_triangular_array(1:pos_x_neg_y_pos_z_dir_z_length))
    allocate(pos_x_neg_y_pos_z_dir_z_triangular_x_array(1:pos_x_neg_y_pos_z_dir_z_length))
    allocate(pos_x_neg_y_pos_z_dir_z_triangular_y_array(1:pos_x_neg_y_pos_z_dir_z_length))
    allocate(pos_x_neg_y_pos_z_dir_z_triangular_z_array(1:pos_x_neg_y_pos_z_dir_z_length))	
    count = 1
    do i = 0, length_pos_x
      do j = 0, length_neg_y, -1
        do k = -j, length_pos_z
          if (k.ge.i .or. -j.ge.i) then
            if (pos_x_neg_y_pos_z_octant(i,j,k).eq.1) then
              pos_x_neg_y_pos_z_dir_z_triangular_array(count) = 1
              pos_x_neg_y_pos_z_dir_z_triangular_x_array(count) = i
              pos_x_neg_y_pos_z_dir_z_triangular_y_array(count) = j
              pos_x_neg_y_pos_z_dir_z_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
!$OMP SECTION		
    do j = 0, length_neg_y, -1
      do k = 0, length_neg_z, -1
        do i = -k, length_pos_x
          if (i.ge.-j .or. -k.ge.-j) then
            if (pos_x_neg_y_neg_z_octant(i,j,k).eq.1) then
              pos_x_neg_y_neg_z_dir_x_length = pos_x_neg_y_neg_z_dir_x_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(pos_x_neg_y_neg_z_dir_x_triangular_array(1:pos_x_neg_y_neg_z_dir_x_length))
    allocate(pos_x_neg_y_neg_z_dir_x_triangular_x_array(1:pos_x_neg_y_neg_z_dir_x_length))
    allocate(pos_x_neg_y_neg_z_dir_x_triangular_y_array(1:pos_x_neg_y_neg_z_dir_x_length))
    allocate(pos_x_neg_y_neg_z_dir_x_triangular_z_array(1:pos_x_neg_y_neg_z_dir_x_length))		
    count = 1
    do j = 0, length_neg_y, -1
      do k = 0, length_neg_z, -1
        do i = -k, length_pos_x
          if (i.ge.-j .or. -k.ge.-j) then
            if (pos_x_neg_y_neg_z_octant(i,j,k).eq.1) then
              pos_x_neg_y_neg_z_dir_x_triangular_array(count) = 1
              pos_x_neg_y_neg_z_dir_x_triangular_x_array(count) = i
              pos_x_neg_y_neg_z_dir_x_triangular_y_array(count) = j
              pos_x_neg_y_neg_z_dir_x_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
		
    do k = 0, length_neg_z, -1
      do i = 0, length_pos_x
        do j = -i, length_neg_y, -1
          if (-j.ge.-k .or. i.ge.-k) then
            if (pos_x_neg_y_neg_z_octant(i,j,k).eq.1) then
              pos_x_neg_y_neg_z_dir_y_length = pos_x_neg_y_neg_z_dir_y_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(pos_x_neg_y_neg_z_dir_y_triangular_array(1:pos_x_neg_y_neg_z_dir_y_length))
    allocate(pos_x_neg_y_neg_z_dir_y_triangular_x_array(1:pos_x_neg_y_neg_z_dir_y_length))
    allocate(pos_x_neg_y_neg_z_dir_y_triangular_y_array(1:pos_x_neg_y_neg_z_dir_y_length))
    allocate(pos_x_neg_y_neg_z_dir_y_triangular_z_array(1:pos_x_neg_y_neg_z_dir_y_length))	
    count = 1
    do k = 0, length_neg_z, -1
      do i = 0, length_pos_x
        do j = -i, length_neg_y, -1
          if (-j.ge.-k .or. i.ge.-k) then
            if (pos_x_neg_y_neg_z_octant(i,j,k).eq.1) then
              pos_x_neg_y_neg_z_dir_y_triangular_array(count) = 1
              pos_x_neg_y_neg_z_dir_y_triangular_x_array(count) = i
              pos_x_neg_y_neg_z_dir_y_triangular_y_array(count) = j
              pos_x_neg_y_neg_z_dir_y_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
	
    do i = 0, length_pos_x
      do j = 0, length_neg_y, -1
        do k = j, length_neg_z, -1
          if (-k.ge.i .or. -j.ge.i) then
            if (pos_x_neg_y_neg_z_octant(i,j,k).eq.1) then
                pos_x_neg_y_neg_z_dir_z_length = pos_x_neg_y_neg_z_dir_z_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(pos_x_neg_y_neg_z_dir_z_triangular_array(1:pos_x_neg_y_neg_z_dir_z_length))
    allocate(pos_x_neg_y_neg_z_dir_z_triangular_x_array(1:pos_x_neg_y_neg_z_dir_z_length))
    allocate(pos_x_neg_y_neg_z_dir_z_triangular_y_array(1:pos_x_neg_y_neg_z_dir_z_length))
    allocate(pos_x_neg_y_neg_z_dir_z_triangular_z_array(1:pos_x_neg_y_neg_z_dir_z_length))	
    count = 1
    do i = 0, length_pos_x
      do j = 0, length_neg_y, -1
        do k = j, length_neg_z, -1
          if (-k.ge.i .or. -j.ge.i) then
            if (pos_x_neg_y_neg_z_octant(i,j,k).eq.1) then
              pos_x_neg_y_neg_z_dir_z_triangular_array(count) = 1
              pos_x_neg_y_neg_z_dir_z_triangular_x_array(count) = i
              pos_x_neg_y_neg_z_dir_z_triangular_y_array(count) = j
              pos_x_neg_y_neg_z_dir_z_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
!$OMP SECTION
    do j = 0, length_pos_y
      do k = 0, length_pos_z
        do i = -k, length_neg_x, -1
          if (-i.ge.j .or. k.ge.j) then
            if (neg_x_pos_y_pos_z_octant(i,j,k).eq.1) then
              neg_x_pos_y_pos_z_dir_x_length = neg_x_pos_y_pos_z_dir_x_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(neg_x_pos_y_pos_z_dir_x_triangular_array(1:neg_x_pos_y_pos_z_dir_x_length))
    allocate(neg_x_pos_y_pos_z_dir_x_triangular_x_array(1:neg_x_pos_y_pos_z_dir_x_length))
    allocate(neg_x_pos_y_pos_z_dir_x_triangular_y_array(1:neg_x_pos_y_pos_z_dir_x_length))
    allocate(neg_x_pos_y_pos_z_dir_x_triangular_z_array(1:neg_x_pos_y_pos_z_dir_x_length))		
    count = 1
    do j = 0, length_pos_y
      do k = 0, length_pos_z
        do i = -k, length_neg_x, -1
          if (-i.ge.j .or. k.ge.j) then
            if (neg_x_pos_y_pos_z_octant(i,j,k).eq.1) then
              neg_x_pos_y_pos_z_dir_x_triangular_array(count) = 1
              neg_x_pos_y_pos_z_dir_x_triangular_x_array(count) = i
              neg_x_pos_y_pos_z_dir_x_triangular_y_array(count) = j
              neg_x_pos_y_pos_z_dir_x_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
	
    do k = 0, length_pos_z
      do i = 0, length_neg_x, -1
        do j = -i, length_pos_y
          if (j.ge.k .or. -i.ge.k) then
            if (neg_x_pos_y_pos_z_octant(i,j,k).eq.1) then
              neg_x_pos_y_pos_z_dir_y_length = neg_x_pos_y_pos_z_dir_y_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(neg_x_pos_y_pos_z_dir_y_triangular_array(1:neg_x_pos_y_pos_z_dir_y_length))
    allocate(neg_x_pos_y_pos_z_dir_y_triangular_x_array(1:neg_x_pos_y_pos_z_dir_y_length))
    allocate(neg_x_pos_y_pos_z_dir_y_triangular_y_array(1:neg_x_pos_y_pos_z_dir_y_length))
    allocate(neg_x_pos_y_pos_z_dir_y_triangular_z_array(1:neg_x_pos_y_pos_z_dir_y_length))	
    count = 1
    do k = 0, length_pos_z
      do i = 0, length_neg_x, -1
        do j = -i, length_pos_y
          if (j.ge.k .or. -i.ge.k) then
            if (neg_x_pos_y_pos_z_octant(i,j,k).eq.1) then
              neg_x_pos_y_pos_z_dir_y_triangular_array(count) = 1
              neg_x_pos_y_pos_z_dir_y_triangular_x_array(count) = i
              neg_x_pos_y_pos_z_dir_y_triangular_y_array(count) = j
              neg_x_pos_y_pos_z_dir_y_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
	
    do i = 0, length_neg_x, -1
      do j = 0, length_pos_y
        do k = j, length_pos_z
          if (k.ge.-i .or. j.ge.-i) then
            if (neg_x_pos_y_pos_z_octant(i,j,k).eq.1) then
              neg_x_pos_y_pos_z_dir_z_length = neg_x_pos_y_pos_z_dir_z_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(neg_x_pos_y_pos_z_dir_z_triangular_array(1:neg_x_pos_y_pos_z_dir_z_length))
    allocate(neg_x_pos_y_pos_z_dir_z_triangular_x_array(1:neg_x_pos_y_pos_z_dir_z_length))
    allocate(neg_x_pos_y_pos_z_dir_z_triangular_y_array(1:neg_x_pos_y_pos_z_dir_z_length))
    allocate(neg_x_pos_y_pos_z_dir_z_triangular_z_array(1:neg_x_pos_y_pos_z_dir_z_length))	
    count = 1
    do i = 0, length_neg_x, -1
      do j = 0, length_pos_y
        do k = j, length_pos_z
          if (k.ge.-i .or. j.ge.-i) then
            if (neg_x_pos_y_pos_z_octant(i,j,k).eq.1) then
              neg_x_pos_y_pos_z_dir_z_triangular_array(count) = 1
              neg_x_pos_y_pos_z_dir_z_triangular_x_array(count) = i
              neg_x_pos_y_pos_z_dir_z_triangular_y_array(count) = j
              neg_x_pos_y_pos_z_dir_z_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
!$OMP SECTION		
    do j = 0, length_pos_y
      do k = 0, length_neg_z, -1
        do i = k, length_neg_x, -1
          if (-i.ge.j .or. -k.ge.j) then
            if (neg_x_pos_y_neg_z_octant(i,j,k).eq.1) then
              neg_x_pos_y_neg_z_dir_x_length = neg_x_pos_y_neg_z_dir_x_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(neg_x_pos_y_neg_z_dir_x_triangular_array(1:neg_x_pos_y_neg_z_dir_x_length))
    allocate(neg_x_pos_y_neg_z_dir_x_triangular_x_array(1:neg_x_pos_y_neg_z_dir_x_length))
    allocate(neg_x_pos_y_neg_z_dir_x_triangular_y_array(1:neg_x_pos_y_neg_z_dir_x_length))
    allocate(neg_x_pos_y_neg_z_dir_x_triangular_z_array(1:neg_x_pos_y_neg_z_dir_x_length))		
    count = 1
    do j = 0, length_pos_y
      do k = 0, length_neg_z, -1
        do i = k, length_neg_x, -1
          if (-i.ge.j .or. -k.ge.j) then
            if (neg_x_pos_y_neg_z_octant(i,j,k).eq.1) then
              neg_x_pos_y_neg_z_dir_x_triangular_array(count) = 1
              neg_x_pos_y_neg_z_dir_x_triangular_x_array(count) = i
              neg_x_pos_y_neg_z_dir_x_triangular_y_array(count) = j
              neg_x_pos_y_neg_z_dir_x_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
		
    do k = 0, length_neg_z, -1
      do i = 0, length_neg_x, -1
        do j = -i, length_pos_y
          if (j.ge.-k .or. -i.ge.-k) then
            if (neg_x_pos_y_neg_z_octant(i,j,k).eq.1) then
              neg_x_pos_y_neg_z_dir_y_length = neg_x_pos_y_neg_z_dir_y_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(neg_x_pos_y_neg_z_dir_y_triangular_array(1:neg_x_pos_y_neg_z_dir_y_length))
    allocate(neg_x_pos_y_neg_z_dir_y_triangular_x_array(1:neg_x_pos_y_neg_z_dir_y_length))
    allocate(neg_x_pos_y_neg_z_dir_y_triangular_y_array(1:neg_x_pos_y_neg_z_dir_y_length))
    allocate(neg_x_pos_y_neg_z_dir_y_triangular_z_array(1:neg_x_pos_y_neg_z_dir_y_length))	
    count = 1
    do k = 0, length_neg_z, -1
      do i = 0, length_neg_x, -1
        do j = -i, length_pos_y
          if (j.ge.-k .or. -i.ge.-k) then
            if (neg_x_pos_y_neg_z_octant(i,j,k).eq.1) then
              neg_x_pos_y_neg_z_dir_y_triangular_array(count) = 1
              neg_x_pos_y_neg_z_dir_y_triangular_x_array(count) = i
              neg_x_pos_y_neg_z_dir_y_triangular_y_array(count) = j
              neg_x_pos_y_neg_z_dir_y_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
	
    do i = 0, length_neg_x, -1
      do j = 0, length_pos_y
        do k = -j, length_neg_z, -1
          if (-k.ge.-i .or. j.ge.-i) then
            if (neg_x_pos_y_neg_z_octant(i,j,k).eq.1) then
              neg_x_pos_y_neg_z_dir_z_length = neg_x_pos_y_neg_z_dir_z_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(neg_x_pos_y_neg_z_dir_z_triangular_array(1:neg_x_pos_y_neg_z_dir_z_length))
    allocate(neg_x_pos_y_neg_z_dir_z_triangular_x_array(1:neg_x_pos_y_neg_z_dir_z_length))
    allocate(neg_x_pos_y_neg_z_dir_z_triangular_y_array(1:neg_x_pos_y_neg_z_dir_z_length))
    allocate(neg_x_pos_y_neg_z_dir_z_triangular_z_array(1:neg_x_pos_y_neg_z_dir_z_length))	
    count = 1
    do i = 0, length_neg_x, -1
      do j = 0, length_pos_y
        do k = -j, length_neg_z, -1
          if (-k.ge.-i .or. j.ge.-i) then
            if (neg_x_pos_y_neg_z_octant(i,j,k).eq.1) then
              neg_x_pos_y_neg_z_dir_z_triangular_array(count) = 1
              neg_x_pos_y_neg_z_dir_z_triangular_x_array(count) = i
              neg_x_pos_y_neg_z_dir_z_triangular_y_array(count) = j
              neg_x_pos_y_neg_z_dir_z_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
!$OMP SECTION		
    do j = 0, length_neg_y, -1
      do k = 0, length_pos_z
        do i = -k, length_neg_x, -1
          if (-i.ge.-j .or. k.ge.-j) then
            if (neg_x_neg_y_pos_z_octant(i,j,k).eq.1) then
              neg_x_neg_y_pos_z_dir_x_length = neg_x_neg_y_pos_z_dir_x_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(neg_x_neg_y_pos_z_dir_x_triangular_array(1:neg_x_neg_y_pos_z_dir_x_length))
    allocate(neg_x_neg_y_pos_z_dir_x_triangular_x_array(1:neg_x_neg_y_pos_z_dir_x_length))
    allocate(neg_x_neg_y_pos_z_dir_x_triangular_y_array(1:neg_x_neg_y_pos_z_dir_x_length))
    allocate(neg_x_neg_y_pos_z_dir_x_triangular_z_array(1:neg_x_neg_y_pos_z_dir_x_length))		
    count = 1
    do j = 0, length_neg_y, -1
      do k = 0, length_pos_z
        do i = -k, length_neg_x, -1
          if (-i.ge.-j .or. k.ge.-j) then
            if (neg_x_neg_y_pos_z_octant(i,j,k).eq.1) then
              neg_x_neg_y_pos_z_dir_x_triangular_array(count) = 1
              neg_x_neg_y_pos_z_dir_x_triangular_x_array(count) = i
              neg_x_neg_y_pos_z_dir_x_triangular_y_array(count) = j
              neg_x_neg_y_pos_z_dir_x_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
	
    do k = 0, length_pos_z
      do i = 0, length_neg_x, -1
        do j = i, length_neg_y, -1
          if (-j.ge.k .or. -i.ge.k) then
            if (neg_x_neg_y_pos_z_octant(i,j,k).eq.1) then
              neg_x_neg_y_pos_z_dir_y_length = neg_x_neg_y_pos_z_dir_y_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(neg_x_neg_y_pos_z_dir_y_triangular_array(1:neg_x_neg_y_pos_z_dir_y_length))
    allocate(neg_x_neg_y_pos_z_dir_y_triangular_x_array(1:neg_x_neg_y_pos_z_dir_y_length))
    allocate(neg_x_neg_y_pos_z_dir_y_triangular_y_array(1:neg_x_neg_y_pos_z_dir_y_length))
    allocate(neg_x_neg_y_pos_z_dir_y_triangular_z_array(1:neg_x_neg_y_pos_z_dir_y_length))	
    count = 1
    do k = 0, length_pos_z
      do i = 0, length_neg_x, -1
        do j = i, length_neg_y, -1
          if (-j.ge.k .or. -i.ge.k) then
            if (neg_x_neg_y_pos_z_octant(i,j,k).eq.1) then
              neg_x_neg_y_pos_z_dir_y_triangular_array(count) = 1
              neg_x_neg_y_pos_z_dir_y_triangular_x_array(count) = i
              neg_x_neg_y_pos_z_dir_y_triangular_y_array(count) = j
              neg_x_neg_y_pos_z_dir_y_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
	
    do i = 0, length_neg_x, -1
      do j = 0, length_neg_y, -1
        do k = -j, length_pos_z
          if (k.ge.-i .or. -j.ge.-i) then
            if (neg_x_neg_y_pos_z_octant(i,j,k).eq.1) then
              neg_x_neg_y_pos_z_dir_z_length = neg_x_neg_y_pos_z_dir_z_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(neg_x_neg_y_pos_z_dir_z_triangular_array(1:neg_x_neg_y_pos_z_dir_z_length))
    allocate(neg_x_neg_y_pos_z_dir_z_triangular_x_array(1:neg_x_neg_y_pos_z_dir_z_length))
    allocate(neg_x_neg_y_pos_z_dir_z_triangular_y_array(1:neg_x_neg_y_pos_z_dir_z_length))
    allocate(neg_x_neg_y_pos_z_dir_z_triangular_z_array(1:neg_x_neg_y_pos_z_dir_z_length))	
    count = 1
    do i = 0, length_neg_x, -1
      do j = 0, length_neg_y, -1
        do k = -j, length_pos_z
          if (k.ge.-i .or. -j.ge.-i) then
            if (neg_x_neg_y_pos_z_octant(i,j,k).eq.1) then
              neg_x_neg_y_pos_z_dir_z_triangular_array(count) = 1
              neg_x_neg_y_pos_z_dir_z_triangular_x_array(count) = i
              neg_x_neg_y_pos_z_dir_z_triangular_y_array(count) = j
              neg_x_neg_y_pos_z_dir_z_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
!$OMP SECTION		
    do j = 0, length_neg_y, -1
      do k = 0, length_neg_z, -1
        do i = k, length_neg_x, -1
          if (-i.ge.-j .or. -k.ge.-j) then
            if (neg_x_neg_y_neg_z_octant(i,j,k).eq.1) then
              neg_x_neg_y_neg_z_dir_x_length = neg_x_neg_y_neg_z_dir_x_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(neg_x_neg_y_neg_z_dir_x_triangular_array(1:neg_x_neg_y_neg_z_dir_x_length))
    allocate(neg_x_neg_y_neg_z_dir_x_triangular_x_array(1:neg_x_neg_y_neg_z_dir_x_length))
    allocate(neg_x_neg_y_neg_z_dir_x_triangular_y_array(1:neg_x_neg_y_neg_z_dir_x_length))
    allocate(neg_x_neg_y_neg_z_dir_x_triangular_z_array(1:neg_x_neg_y_neg_z_dir_x_length))		
    count = 1
    do j = 0, length_neg_y, -1
      do k = 0, length_neg_z, -1
        do i =k , length_neg_x, -1
          if (-i.ge.-j .or. -k.ge.-j) then
            if (neg_x_neg_y_neg_z_octant(i,j,k).eq.1) then
              neg_x_neg_y_neg_z_dir_x_triangular_array(count) = 1
              neg_x_neg_y_neg_z_dir_x_triangular_x_array(count) = i
              neg_x_neg_y_neg_z_dir_x_triangular_y_array(count) = j
              neg_x_neg_y_neg_z_dir_x_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
		
    do k = 0, length_neg_z, -1
      do i = 0, length_neg_x, -1
        do j = i, length_neg_y, -1
          if (-j.ge.-k .or. -i.ge.-k) then
            if (neg_x_neg_y_neg_z_octant(i,j,k).eq.1) then
              neg_x_neg_y_neg_z_dir_y_length = neg_x_neg_y_neg_z_dir_y_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(neg_x_neg_y_neg_z_dir_y_triangular_array(1:neg_x_neg_y_neg_z_dir_y_length))
    allocate(neg_x_neg_y_neg_z_dir_y_triangular_x_array(1:neg_x_neg_y_neg_z_dir_y_length))
    allocate(neg_x_neg_y_neg_z_dir_y_triangular_y_array(1:neg_x_neg_y_neg_z_dir_y_length))
    allocate(neg_x_neg_y_neg_z_dir_y_triangular_z_array(1:neg_x_neg_y_neg_z_dir_y_length))	
    count = 1
    do k = 0, length_neg_z, -1
      do i = 0, length_neg_x, -1
        do j = i, length_neg_y, -1
          if (-j.ge.-k .or. -i.ge.-k) then
            if (neg_x_neg_y_neg_z_octant(i,j,k).eq.1) then
              neg_x_neg_y_neg_z_dir_y_triangular_array(count) = 1
              neg_x_neg_y_neg_z_dir_y_triangular_x_array(count) = i
              neg_x_neg_y_neg_z_dir_y_triangular_y_array(count) = j
              neg_x_neg_y_neg_z_dir_y_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo
	
    do i = 0, length_neg_x, -1
      do j = 0, length_neg_y, -1
        do k = j, length_neg_z, -1
          if (-k.ge.-i .or. -j.ge.-i) then
            if (neg_x_neg_y_neg_z_octant(i,j,k).eq.1) then
              neg_x_neg_y_neg_z_dir_z_length = neg_x_neg_y_neg_z_dir_z_length + 1
            endif
          endif
        enddo	
      enddo
    enddo	
    allocate(neg_x_neg_y_neg_z_dir_z_triangular_array(1:neg_x_neg_y_neg_z_dir_z_length))
    allocate(neg_x_neg_y_neg_z_dir_z_triangular_x_array(1:neg_x_neg_y_neg_z_dir_z_length))
    allocate(neg_x_neg_y_neg_z_dir_z_triangular_y_array(1:neg_x_neg_y_neg_z_dir_z_length))
    allocate(neg_x_neg_y_neg_z_dir_z_triangular_z_array(1:neg_x_neg_y_neg_z_dir_z_length))	
    count = 1
    do i = 0, length_neg_x, -1
      do j = 0, length_neg_y, -1
        do k = j, length_neg_z, -1
          if (-k.ge.-i .or. -j.ge.-i) then
            if (neg_x_neg_y_neg_z_octant(i,j,k).eq.1) then
              neg_x_neg_y_neg_z_dir_z_triangular_array(count) = 1
              neg_x_neg_y_neg_z_dir_z_triangular_x_array(count) = i
              neg_x_neg_y_neg_z_dir_z_triangular_y_array(count) = j
              neg_x_neg_y_neg_z_dir_z_triangular_z_array(count) = k
              count = count + 1
            endif
          endif
        enddo	
      enddo
    enddo

!$OMP END SECTIONS

!$OMP END PARALLEL 

  end subroutine AT_triangular_grid_initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module AT_triangular_grid
