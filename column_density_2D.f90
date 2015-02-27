module column_density_2D

  use precision, only: dp 
  use parameter, only: mesh, do_periodic
  use c2ray_parameters, only: epsilon_value
  use array, only: intermediate_average_xHI, &
                   intermediate_average_xHeI, &
                   intermediate_average_xHeII, &
                   xHI_array, xHeI_array, xHeII_array
  use sourceprops, only: srcpos
  use short, only: short_characteristic
  use array, only: global_LTE_array

contains

  subroutine column_density_2D_plane(source_index,i_plane,case_of_evolution,ADP_help_do_LTE_ray_tracing)

    implicit none

    integer,intent(in) :: source_index 
    integer,intent(in) :: i_plane 
    character,intent(in) :: case_of_evolution
    logical,intent(in) :: ADP_help_do_LTE_ray_tracing
    integer,dimension(1:3) :: ray_traced_pos 
    integer,dimension(1:3) :: global_pos
    integer,dimension(1:3) :: neg_dir_from_source_pos
    integer,dimension(1:3) :: pos_dir_from_source_pos
    real(kind=dp) :: xHI,xHeI,xHeII
    integer :: i,j,k

    ! The left/right, up/down, front/back position limit assuming infinite box
    if (do_periodic) then
      neg_dir_from_source_pos(:) = srcpos(:,source_index)-mesh/2
      pos_dir_from_source_pos(:) = srcpos(:,source_index)+mesh/2-1+mod(mesh,2)
    else
      neg_dir_from_source_pos(:) = 1
      pos_dir_from_source_pos(:) = mesh
    endif

    select case (i_plane)

    ! Sweep in +i,+j direction
    case(1)

      ray_traced_pos(3) = srcpos(3,source_index)
      global_pos(3) = modulo(ray_traced_pos(3)-1,mesh)+1

      do j = srcpos(2,source_index)+1,pos_dir_from_source_pos(2)
        do i = srcpos(1,source_index)+1,pos_dir_from_source_pos(1)

          ray_traced_pos(1) = i
          global_pos(1) = modulo(ray_traced_pos(1)-1,mesh)+1
          ray_traced_pos(2) = j
          global_pos(2) = modulo(ray_traced_pos(2)-1,mesh)+1

    if (case_of_evolution.eq.'L') then
      xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
    endif

        if (case_of_evolution.eq.'A') then

          if (ADP_help_do_LTE_ray_tracing .eqv. .false.) then
            xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
          else

            if (global_LTE_array(global_pos(1),global_pos(2),global_pos(3)) .ge. 0.0) then

              xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            else

              xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            endif

          endif

        endif

        call short_characteristic(ray_traced_pos,srcpos(:,source_index),global_pos,xHI,xHeI,xHeII)  

        enddo
      enddo

    ! Sweep in +i,-j direction
    case(2)

      ray_traced_pos(3) = srcpos(3,source_index)
      global_pos(3) = modulo(ray_traced_pos(3)-1,mesh)+1

      do j = srcpos(2,source_index)-1,neg_dir_from_source_pos(2),-1
        do i = srcpos(1,source_index)+1,pos_dir_from_source_pos(1)

          ray_traced_pos(1) = i
          global_pos(1) = modulo(ray_traced_pos(1)-1,mesh)+1
          ray_traced_pos(2) = j
          global_pos(2) = modulo(ray_traced_pos(2)-1,mesh)+1

    if (case_of_evolution.eq.'L') then
      xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
    endif

        if (case_of_evolution.eq.'A') then

          if (ADP_help_do_LTE_ray_tracing .eqv. .false.) then
            xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
          else

            if (global_LTE_array(global_pos(1),global_pos(2),global_pos(3)) .ge. 0.0) then

              xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            else

              xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            endif

          endif

        endif

        call short_characteristic(ray_traced_pos,srcpos(:,source_index),global_pos,xHI,xHeI,xHeII)  

        enddo
      enddo

    ! Sweep in -i,+j direction
    case(3)

      ray_traced_pos(3) = srcpos(3,source_index)
      global_pos(3) = modulo(ray_traced_pos(3)-1,mesh)+1

      do j = srcpos(2,source_index)+1,pos_dir_from_source_pos(2)
        do i = srcpos(1,source_index)-1,neg_dir_from_source_pos(1),-1

          ray_traced_pos(1) = i
          global_pos(1) = modulo(ray_traced_pos(1)-1,mesh)+1
          ray_traced_pos(2) = j
          global_pos(2) = modulo(ray_traced_pos(2)-1,mesh)+1

    if (case_of_evolution.eq.'L') then
      xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
    endif

        if (case_of_evolution.eq.'A') then

          if (ADP_help_do_LTE_ray_tracing .eqv. .false.) then
            xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
          else

            if (global_LTE_array(global_pos(1),global_pos(2),global_pos(3)) .ge. 0.0) then

              xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            else

              xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            endif

          endif

        endif

        call short_characteristic(ray_traced_pos,srcpos(:,source_index),global_pos,xHI,xHeI,xHeII)  

        enddo
      enddo

    ! Sweep in -i,-j direction
    case(4)

      ray_traced_pos(3) = srcpos(3,source_index)
      global_pos(3) = modulo(ray_traced_pos(3)-1,mesh)+1

      do j = srcpos(2,source_index)-1,neg_dir_from_source_pos(2),-1
        do i = srcpos(1,source_index)-1,neg_dir_from_source_pos(1),-1

          ray_traced_pos(1) = i
          global_pos(1) = modulo(ray_traced_pos(1)-1,mesh)+1
          ray_traced_pos(2) = j
          global_pos(2) = modulo(ray_traced_pos(2)-1,mesh)+1

    if (case_of_evolution.eq.'L') then
      xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
    endif

        if (case_of_evolution.eq.'A') then

          if (ADP_help_do_LTE_ray_tracing .eqv. .false.) then
            xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
          else

            if (global_LTE_array(global_pos(1),global_pos(2),global_pos(3)) .ge. 0.0) then

              xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            else

              xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            endif

          endif

        endif

        call short_characteristic(ray_traced_pos,srcpos(:,source_index),global_pos,xHI,xHeI,xHeII)  

        enddo
      enddo

    ! Sweep in +i,+k direction
    case(5)

      ray_traced_pos(2) = srcpos(2,source_index)
      global_pos(2) = modulo(ray_traced_pos(2)-1,mesh)+1

      do k = srcpos(3,source_index)+1,pos_dir_from_source_pos(3)
        do i = srcpos(1,source_index)+1,pos_dir_from_source_pos(1)

          ray_traced_pos(1) = i
          global_pos(1) = modulo(ray_traced_pos(1)-1,mesh)+1
          ray_traced_pos(3) = k
          global_pos(3) = modulo(ray_traced_pos(3)-1,mesh)+1

    if (case_of_evolution.eq.'L') then
      xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
    endif

        if (case_of_evolution.eq.'A') then

          if (ADP_help_do_LTE_ray_tracing .eqv. .false.) then
            xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
          else

            if (global_LTE_array(global_pos(1),global_pos(2),global_pos(3)) .ge. 0.0) then

              xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            else

              xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            endif

          endif

        endif

        call short_characteristic(ray_traced_pos,srcpos(:,source_index),global_pos,xHI,xHeI,xHeII)  

        enddo
      enddo

    ! Sweep in -i,+k direction
    case(6)

      ray_traced_pos(2) = srcpos(2,source_index)
      global_pos(2) = modulo(ray_traced_pos(2)-1,mesh)+1

      do k = srcpos(3,source_index)+1,pos_dir_from_source_pos(3)
        do i = srcpos(1,source_index)-1,neg_dir_from_source_pos(1),-1

          ray_traced_pos(1) = i
          global_pos(1) = modulo(ray_traced_pos(1)-1,mesh)+1
          ray_traced_pos(3) = k
          global_pos(3) = modulo(ray_traced_pos(3)-1,mesh)+1

    if (case_of_evolution.eq.'L') then
      xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
    endif

        if (case_of_evolution.eq.'A') then

          if (ADP_help_do_LTE_ray_tracing .eqv. .false.) then
            xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
          else

            if (global_LTE_array(global_pos(1),global_pos(2),global_pos(3)) .ge. 0.0) then

              xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            else

              xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            endif

          endif

        endif

        call short_characteristic(ray_traced_pos,srcpos(:,source_index),global_pos,xHI,xHeI,xHeII)  

        enddo
      enddo

    ! Sweep in -i,-k direction
    case(7)

      ray_traced_pos(2) = srcpos(2,source_index)
      global_pos(2) = modulo(ray_traced_pos(2)-1,mesh)+1

      do k = srcpos(3,source_index)-1,neg_dir_from_source_pos(3),-1
        do i = srcpos(1,source_index)-1,neg_dir_from_source_pos(1),-1

          ray_traced_pos(1) = i
          global_pos(1) = modulo(ray_traced_pos(1)-1,mesh)+1
          ray_traced_pos(3) = k
          global_pos(3) = modulo(ray_traced_pos(3)-1,mesh)+1

    if (case_of_evolution.eq.'L') then
      xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
    endif

        if (case_of_evolution.eq.'A') then

          if (ADP_help_do_LTE_ray_tracing .eqv. .false.) then
            xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
          else

            if (global_LTE_array(global_pos(1),global_pos(2),global_pos(3)) .ge. 0.0) then

              xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            else

              xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            endif

          endif

        endif

        call short_characteristic(ray_traced_pos,srcpos(:,source_index),global_pos,xHI,xHeI,xHeII)  

        enddo
      enddo

    ! Sweep in +i,-k direction
    case(8)

      ray_traced_pos(2) = srcpos(2,source_index)
      global_pos(2) = modulo(ray_traced_pos(2)-1,mesh)+1

      do k = srcpos(3,source_index)-1,neg_dir_from_source_pos(3),-1
        do i = srcpos(1,source_index)+1,pos_dir_from_source_pos(1)

          ray_traced_pos(1) = i
          global_pos(1) = modulo(ray_traced_pos(1)-1,mesh)+1
          ray_traced_pos(3) = k
          global_pos(3) = modulo(ray_traced_pos(3)-1,mesh)+1

    if (case_of_evolution.eq.'L') then
      xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
    endif

        if (case_of_evolution.eq.'A') then

          if (ADP_help_do_LTE_ray_tracing .eqv. .false.) then
            xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
          else

            if (global_LTE_array(global_pos(1),global_pos(2),global_pos(3)) .ge. 0.0) then

              xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            else

              xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            endif

          endif

        endif

        call short_characteristic(ray_traced_pos,srcpos(:,source_index),global_pos,xHI,xHeI,xHeII)  

        enddo
      enddo

    ! Sweep in +j,+k direction
    case(9) 

      ray_traced_pos(1) = srcpos(1,source_index)
      global_pos(1) = modulo(ray_traced_pos(1)-1,mesh)+1

      do k = srcpos(3,source_index)+1,pos_dir_from_source_pos(3)
        do j = srcpos(2,source_index)+1,pos_dir_from_source_pos(2)

          ray_traced_pos(2) = j
          global_pos(2) = modulo(ray_traced_pos(2)-1,mesh)+1
          ray_traced_pos(3) = k
          global_pos(3) = modulo(ray_traced_pos(3)-1,mesh)+1

    if (case_of_evolution.eq.'L') then
      xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
    endif

        if (case_of_evolution.eq.'A') then

          if (ADP_help_do_LTE_ray_tracing .eqv. .false.) then
            xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
          else

            if (global_LTE_array(global_pos(1),global_pos(2),global_pos(3)) .ge. 0.0) then

              xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            else

              xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            endif

          endif

        endif

        call short_characteristic(ray_traced_pos,srcpos(:,source_index),global_pos,xHI,xHeI,xHeII)  

        enddo
      enddo

    ! Sweep in -j,+k direction
    case(10) 

      ray_traced_pos(1) = srcpos(1,source_index)
      global_pos(1) = modulo(ray_traced_pos(1)-1,mesh)+1

      do k = srcpos(3,source_index)+1,pos_dir_from_source_pos(3)
        do j = srcpos(2,source_index)-1,neg_dir_from_source_pos(2),-1

          ray_traced_pos(2) = j
          global_pos(2) = modulo(ray_traced_pos(2)-1,mesh)+1
          ray_traced_pos(3) = k
          global_pos(3) = modulo(ray_traced_pos(3)-1,mesh)+1

    if (case_of_evolution.eq.'L') then
      xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
    endif

        if (case_of_evolution.eq.'A') then

          if (ADP_help_do_LTE_ray_tracing .eqv. .false.) then
            xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
          else

            if (global_LTE_array(global_pos(1),global_pos(2),global_pos(3)) .ge. 0.0) then

              xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            else

              xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            endif

          endif

        endif

        call short_characteristic(ray_traced_pos,srcpos(:,source_index),global_pos,xHI,xHeI,xHeII)  

        enddo
      enddo

    ! Sweep in +j,-k direction
    case(11) 

      ray_traced_pos(1) = srcpos(1,source_index)
      global_pos(1) = modulo(ray_traced_pos(1)-1,mesh)+1

      do k = srcpos(3,source_index)-1,neg_dir_from_source_pos(3),-1
        do j = srcpos(2,source_index)+1,pos_dir_from_source_pos(2)

          ray_traced_pos(2) = j
          global_pos(2) = modulo(ray_traced_pos(2)-1,mesh)+1
          ray_traced_pos(3) = k
          global_pos(3) = modulo(ray_traced_pos(3)-1,mesh)+1

    if (case_of_evolution.eq.'L') then
      xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
    endif

        if (case_of_evolution.eq.'A') then

          if (ADP_help_do_LTE_ray_tracing .eqv. .false.) then
            xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
          else

            if (global_LTE_array(global_pos(1),global_pos(2),global_pos(3)) .ge. 0.0) then

              xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            else

              xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            endif

          endif

        endif

        call short_characteristic(ray_traced_pos,srcpos(:,source_index),global_pos,xHI,xHeI,xHeII)  

        enddo
      enddo

    ! Sweep in -j,-k direction
    case(12) 

      ray_traced_pos(1) = srcpos(1,source_index)
      global_pos(1) = modulo(ray_traced_pos(1)-1,mesh)+1

      do k = srcpos(3,source_index)-1,neg_dir_from_source_pos(3),-1
        do j = srcpos(2,source_index)-1,neg_dir_from_source_pos(2),-1

          ray_traced_pos(2) = j
          global_pos(2) = modulo(ray_traced_pos(2)-1,mesh)+1
          ray_traced_pos(3) = k
          global_pos(3) = modulo(ray_traced_pos(3)-1,mesh)+1

    if (case_of_evolution.eq.'L') then
      xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
      xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
    endif

        if (case_of_evolution.eq.'A') then

          if (ADP_help_do_LTE_ray_tracing .eqv. .false.) then
            xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
            xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
          else

            if (global_LTE_array(global_pos(1),global_pos(2),global_pos(3)) .ge. 0.0) then

              xHI = max(xHI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(xHeI_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(xHeII_array(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            else

              xHI = max(intermediate_average_xHI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeI = max(intermediate_average_xHeI(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)
              xHeII = max(intermediate_average_xHeII(global_pos(1),global_pos(2),global_pos(3)),epsilon_value)

            endif

          endif

        endif

        call short_characteristic(ray_traced_pos,srcpos(:,source_index),global_pos,xHI,xHeI,xHeII)  

        enddo
      enddo

    end select
  
  end subroutine column_density_2D_plane

end module column_density_2D
