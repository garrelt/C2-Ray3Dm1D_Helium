module AT_transformation
  
  use parameter, only: length_pos_x, length_neg_x, &
                       length_pos_y, length_neg_y, &
                       length_pos_z, length_neg_z, &
                       do_periodic
  use my_mpi
  use parameter, only: mesh
  use array, only: global_Ifront_array, source_Ifront_array
  use sourceprops, only: srcpos	

  implicit none
	
contains
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AT_global_to_source_grid_transformation(ns)

    implicit none

    integer, intent(in) :: ns
    integer :: i_source, j_source, k_source
    integer :: i_global, j_global, k_global
 ! I find a problem here,
 ! if mesh < OMP_get_num_threads(), it hangs
!$OMP PARALLEL PRIVATE(i_global,j_global,k_global)
!$OMP DO SCHEDULE(STATIC,mesh/OMP_get_num_threads())

    do i_source = length_neg_x, length_pos_x
      do j_source = length_neg_y, length_pos_y
        do k_source = length_neg_z, length_pos_z

          if (do_periodic) then
            i_global = modulo(i_source+srcpos(1,ns)-1, mesh) + 1
            j_global = modulo(j_source+srcpos(2,ns)-1, mesh) + 1
            k_global = modulo(k_source+srcpos(3,ns)-1, mesh) + 1
          else
            i_global = i_source + srcpos(1,ns)
            j_global = j_source + srcpos(2,ns)
            k_global = k_source + srcpos(3,ns)
          endif

          source_Ifront_array(i_source,j_source,k_source) = global_Ifront_array(i_global,j_global,k_global)

        enddo
      enddo
    enddo

!$OMP END DO
!$OMP END PARALLEL
 
  end subroutine AT_global_to_source_grid_transformation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module AT_transformation
