module global_photo_rate
 
  use precision, only: dp 
  use my_mpi
  use parameter, only: mesh, do_periodic, do_LTE
  use array, only: global_LTE_array, global_evolved_LTE_array, source_class_index, equivalent_class_index_array
  use array, only: intermediate_average_xHI,xHI_array, &
                   intermediate_average_xHII,xHII_array, &
                   intermediate_average_xHeI,xHeI_array, &
                   intermediate_average_xHeII,xHeII_array
  use sourceprops, only: srcpos
  use c2ray_parameters, only: epsilon_value
  use local_photo_rate, only: photo_rate_computation

contains

  subroutine global_get_photo_rate(source_index,case_of_evolution,ADP_help_do_LTE_ray_tracing)

    implicit none

    integer,intent(in) :: source_index
    character,intent(in) :: case_of_evolution
    logical,intent(in) :: ADP_help_do_LTE_ray_tracing
    integer,dimension(1:3) :: global_pos
    integer,dimension(1:3) :: neg_dir_from_source_pos
    integer,dimension(1:3) :: pos_dir_from_source_pos
    real(kind=dp) :: xHI,xHII,xHeI,xHeII
    logical :: do_photo_rate_computation
    logical :: boundary_i, boundary_j, boundary_k, boundary
    integer :: i, j, k

    ! The left/right, up/down, front/back position limit assuming infinite box
    if (do_periodic) then
      neg_dir_from_source_pos(:) = srcpos(:,source_index)-mesh/2
      pos_dir_from_source_pos(:) = srcpos(:,source_index)+mesh/2-1+mod(mesh,2)
    else
      neg_dir_from_source_pos(:) = 1
      pos_dir_from_source_pos(:) = mesh
    endif

!$OMP PARALLEL PRIVATE(global_pos,i,j,xHI,xHII,xHeI,xHeII,boundary_i,boundary_j,boundary_k,boundary) 
!$OMP DO SCHEDULE (DYNAMIC,1)
    do k = neg_dir_from_source_pos(3),pos_dir_from_source_pos(3)

      boundary_k = .false.
      if (k.eq.neg_dir_from_source_pos(3)) boundary_k = .true.
      if (k.eq.pos_dir_from_source_pos(3)) boundary_k = .true.

      do j = neg_dir_from_source_pos(2),pos_dir_from_source_pos(2)

        boundary_j = .false.
        if (j.eq.neg_dir_from_source_pos(2)) boundary_j = .true.
        if (j.eq.pos_dir_from_source_pos(2)) boundary_j = .true.

        do i = neg_dir_from_source_pos(1),pos_dir_from_source_pos(1)

          boundary_i = .false.
          if (i.eq.neg_dir_from_source_pos(1)) boundary_i = .true.
          if (i.eq.pos_dir_from_source_pos(1)) boundary_i = .true.

         ! Initialize the logical parameter
          boundary = .false.

          if (boundary_i .or. boundary_j .or. boundary_k) then
            boundary = .true.
          endif

          global_pos(1) = modulo(i-1,mesh)+1
          global_pos(2) = modulo(j-1,mesh)+1
          global_pos(3) = modulo(k-1,mesh)+1

          select case (case_of_evolution)

case ('L')

            if (equivalent_class_index_array(global_pos(1),global_pos(2),global_pos(3)).eq.source_class_index(source_index) .and. &
                global_evolved_LTE_array(global_pos(1),global_pos(2),global_pos(3)).eq.0) then

              xHI = xHI_array(global_pos(1),global_pos(2),global_pos(3))
              xHII = xHII_array(global_pos(1),global_pos(2),global_pos(3))
              xHeI = xHeI_array(global_pos(1),global_pos(2),global_pos(3))
              xHeII = xHeII_array(global_pos(1),global_pos(2),global_pos(3))

              call photo_rate_computation('L',global_pos,source_index,xHI,xHII,xHeI,xHeII,boundary)

            endif

case ('A')

if (do_LTE .eqv. .true.) then

 if (ADP_help_do_LTE_ray_tracing .eqv. .false.) then            

       if (global_LTE_array(global_pos(1),global_pos(2),global_pos(3)).lt.0) then

              xHI=max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHII=max(intermediate_average_xHII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI=max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII=max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

              call photo_rate_computation('A',global_pos,source_index,xHI,xHII,xHeI,xHeII,boundary)

        endif

 else if (ADP_help_do_LTE_ray_tracing .eqv. .true.) then

        if (global_LTE_array(global_pos(1),global_pos(2),global_pos(3)) .ge. 0.0) then

              xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHII = max(xHII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

        else

              xHI=max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHII=max(intermediate_average_xHII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI=max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII=max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

        endif

            call photo_rate_computation('A',global_pos,source_index,xHI,xHII,xHeI,xHeII,boundary)

 endif

elseif (do_LTE .eqv. .false.) then

              xHI=max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHII=max(intermediate_average_xHII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI=max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII=max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

              call photo_rate_computation('A',global_pos,source_index,xHI,xHII,xHeI,xHeII,boundary)

endif
          end select

        enddo
      enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

  end subroutine global_get_photo_rate

end module global_photo_rate
