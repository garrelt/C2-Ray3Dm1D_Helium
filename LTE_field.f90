module LTE_field
 
  use precision, only: dp 
  use my_mpi
  use file_admin, only: logf
  use array, only: global_LTE_array, global_evolved_LTE_array, &
                   xHI_array,xHeI_array,xHeII_array,xHeIII_array, number_density_array
  use sourceprops, only: srcpos
  use parameter, only: mesh, LTE_size_of_m, &
                       length_pos_x, length_neg_x, &
                       length_pos_y, length_neg_y, &
                       length_pos_z, length_neg_z,LTE_non_evolved_cells_exist
  use sourceprops, only: number_of_source

  implicit none 
	
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LTE_global_LTE_field_generator(time)

    implicit none

    real(kind=dp), intent(in) :: time
    integer :: i,j,k

!if (global_LTE_array(5,5,5).lt.0) write(*,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
!write(*,*) 'LTE field xHI_array(5,5,5)',xHI_array(5,5,5)
!write(*,*) 'LTE field xHeII_array(5,5,5)',xHeII_array(5,5,5)
    do i = 1,mesh
      do j = 1,mesh
        do k = 1,mesh

          ! for mild source, e.g. 5e4 BB
          !if (xHI_array(i,j,k).le.1.0e-4) global_LTE_array(i,j,k) = 1
          if (global_LTE_array(i,j,k).lt.0) then

          ! for strong source, e.g. power law or 1e5 BB
            if (xHI_array(i,j,k).le.1.0e-4 .and. xHeII_array(i,j,k).le.3.0e-4) then

              global_LTE_array(i,j,k) = time
             !if (i.eq.4.and.j.eq.4.and.k.eq.4) write(*,*) 'L time is', global_LTE_array(i,j,k)/3.1536e13,' Myr'
            endif
          endif

        enddo
      enddo
    enddo

  end subroutine LTE_global_LTE_field_generator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LTE_non_evolved_cells_finder()

    implicit none

    if (LTE_size_of_m.ge.1) then
       LTE_non_evolved_cells_exist = .true.
    else
       LTE_non_evolved_cells_exist = .false.
    endif

    if (rank .eq. 0) then
      write(logf,*) 'LTE non evolved cells exist? ',LTE_non_evolved_cells_exist
      flush(logf) 
    endif

  end subroutine LTE_non_evolved_cells_finder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LTE_global_evolved_LTE_field_generator()

    implicit none

    integer :: i,j,k

    do i = 1,mesh
      do j = 1,mesh
        do k = 1,mesh

          global_evolved_LTE_array(i,j,k) = global_LTE_array(i,j,k)

        enddo
      enddo
    enddo

  end subroutine LTE_global_evolved_LTE_field_generator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module LTE_field
