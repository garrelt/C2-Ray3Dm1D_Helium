module global_column_density

  use precision, only: dp 
  use my_mpi
  use array, only: coldens_in_HI, &
                   coldens_in_HeI, &  
                   coldens_in_HeII, &
                   coldens_out_HI, &
                   coldens_out_HeI, &
                   coldens_out_HeII, & 
                   ring_vol
  use column_density_0D, only: column_density_0D_point
  use column_density_1D, only: column_density_1D_axis
  use column_density_2D, only: column_density_2D_plane
  use column_density_3D, only: column_density_3D_octant

contains

  subroutine global_get_column_density(source_index,case_of_evolution,ADP_help_do_LTE_ray_tracing)

    implicit none

    integer,intent(in) :: source_index
    character,intent(in) :: case_of_evolution    
    logical,intent(in) :: ADP_help_do_LTE_ray_tracing
    integer :: i_axis,i_plane,i_octant

    ! Initialization of some global arrays
    !coldens_in_HI(:,:,:) = 0.0
    !coldens_in_HeI(:,:,:) = 0.0  
    !coldens_in_HeII(:,:,:) = 0.0  
    !coldens_out_HI(:,:,:) = 0.0
    !coldens_out_HeI(:,:,:) = 0.0  
    !coldens_out_HeII(:,:,:) = 0.0  
    !ring_vol(:,:,:) = 0.0

    ! Do the source point
    call column_density_0D_point(source_index,case_of_evolution,ADP_help_do_LTE_ray_tracing)

!$OMP PARALLEL

!$OMP DO SCHEDULE (DYNAMIC,1)
    ! Do along the 6 axis
    do i_axis = 1,6
      call column_density_1D_axis(source_index,i_axis,case_of_evolution,ADP_help_do_LTE_ray_tracing)
    enddo
!$OMP END DO

!$OMP DO SCHEDULE (DYNAMIC,1)
    ! Do across the 12 planes
    do i_plane = 1,12
      call column_density_2D_plane(source_index,i_plane,case_of_evolution,ADP_help_do_LTE_ray_tracing)
    end do
!$OMP END DO

!$OMP DO SCHEDULE (DYNAMIC,1)
    ! Do all the 8 octants
    do i_octant = 1,8
      call column_density_3D_octant(source_index,i_octant,case_of_evolution,ADP_help_do_LTE_ray_tracing)
    end do
!$OMP END DO

!$OMP END PARALLEL
!if (ADP_help_do_LTE_ray_tracing.eqv..true.) then
!write(*,*) 'coldens_in_HI is ',coldens_in_HI(5,5,5)
!endif

  end subroutine global_get_column_density

end module global_column_density
