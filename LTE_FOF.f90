module LTE_FOF

  use precision, only: dp 
  use parameter, only: mesh,mesh3
  use array, only: global_LTE_array,equivalent_class_index_array,xHI_array, source_class_index
  use sourceprops, only: number_of_source, srcpos

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LTE_equivalent_class_finder()

    implicit none

    integer :: i,j,k,m1,m2
    integer :: m
    integer :: ip1,jp1,kp1
    !integer, dimension(1:mesh3) :: class
    !integer, dimension(1:mesh,1:mesh,1:mesh) :: ijk_to_m
    integer, dimension(:), allocatable :: class   
    integer, dimension(:,:,:), allocatable :: ijk_to_m

allocate(class(1:mesh3))
allocate(ijk_to_m(1:mesh,1:mesh,1:mesh))

    m = 1
    do k = 1,mesh
       do j = 1,mesh
          do i = 1,mesh

          ! Link coordinate (i,j,k) to a number m
          ijk_to_m(i,j,k) = m
          ! Initialize the equivalent class index of each cell 
          class(m) = m
          ! Next number
          m = m+1

        enddo
      enddo
    enddo

    ! In search of the connected component cell-by-cell
    do k = 1,mesh
       do j = 1,mesh
          do i = 1,mesh

          ! If the cell is a LTE cell
          if (global_LTE_array(i,j,k).eq.1) then

            ! Find out the equivalent class index of cell (i,j,k) 
            m1 = ijk_to_m(i,j,k)          
            do while (m1 .ne. class(m1))
              m1 = class(m1)
            enddo

            ! In search of the face-to-face neighbouring cells at positive directions
            ip1=i+1
            jp1=j+1
            kp1=k+1

            ! For face-to-face neighbouring cells across the boundary
            if (ip1.gt.mesh) ip1=1
            if (jp1.gt.mesh) jp1=1
            if (kp1.gt.mesh) kp1=1
            
            if (global_LTE_array(i,j,kp1).eq.1) then

              ! Find out the equivalent class index of neighbouring cells (+k direction)
              m2 = ijk_to_m(i,j,kp1)
              do while (m2 .ne. class(m2))
                m2 = class(m2)
              enddo

              ! Set the equivalent class index of neighbouring cells to its own equivalent class index 
              class(m2) = m1

            endif

            if (global_LTE_array(i,jp1,k).eq.1) then

              ! Find out the equivalent class index of neighbouring cells (+j direction)
              m2 = ijk_to_m(i,jp1,k)
              do while (m2 .ne. class(m2))
                m2 = class(m2)
              enddo

              ! Set the equivalent class index of neighbouring cells to its own equivalent class index 
              class(m2) = m1

            endif

            if (global_LTE_array(ip1,j,k).eq.1) then

              ! Find out the equivalent class index of neighbouring cells (+i direction)
              m2 = ijk_to_m(ip1,j,k)
              do while (m2 .ne. class(m2))
                m2 = class(m2)
              enddo

              ! Set the equivalent class index of neighbouring cells to its own equivalent class index 
              class(m2) = m1

            endif

          endif

        enddo
      enddo	
    enddo
     
    ! Loop all the cells and equilize the equivalent class indices of cells in the same equivalent class
    do k = 1,mesh
       do j = 1,mesh
          do i = 1,mesh

          if (global_LTE_array(i,j,k).eq.1) then

            m = ijk_to_m(i,j,k)
            do while (class(m) .ne. class(class(m)))
              class(m) = class(class(m))
            enddo

          endif

        enddo
      enddo
    enddo

    !Copy the result to equivalent class index array
    do k = 1,mesh 
       do j = 1,mesh
          do i = 1,mesh

          if (global_LTE_array(i,j,k).eq.1) then
            equivalent_class_index_array(i,j,k) = class(ijk_to_m(i,j,k))
          endif

        enddo
      enddo
    enddo

deallocate(class)
deallocate(ijk_to_m)

  end subroutine LTE_equivalent_class_finder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LTE_source_class_index_finder()

    implicit none

    integer :: i_source

    do i_source = 1,number_of_source

      source_class_index(i_source) = &
            equivalent_class_index_array(srcpos(1,i_source),srcpos(2,i_source),srcpos(3,i_source))

    enddo

  end subroutine LTE_source_class_index_finder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module LTE_FOF
