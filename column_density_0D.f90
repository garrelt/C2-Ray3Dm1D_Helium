module column_density_0D

  use precision, only: dp 
  use parameter, only: mesh
  use abundances, only: abu_he
  use c2ray_parameters, only: epsilon_value
  use grid, only: cell_size
  use doric_module, only: coldens
  use array, only: intermediate_average_xHI, &
                   intermediate_average_xHeI, &
                   intermediate_average_xHeII, &
                   xHI_array
  use sourceprops, only: srcpos
  use array, only: coldens_in_HI, &
                   coldens_in_HeI, &
                   coldens_in_HeII, &
                   coldens_out_HI, &
                   coldens_out_HeI, &
                   coldens_out_HeII, &
                   ring_vol, &
                   xHI_array, xHeI_array, xHeII_array, number_density_array
  use array, only: global_LTE_array

contains

  subroutine column_density_0D_point(source_index,case_of_evolution,ADP_help_do_LTE_ray_tracing)

    implicit none

    integer,intent(in) :: source_index  
    character,intent(in) :: case_of_evolution
    logical,intent(in) :: ADP_help_do_LTE_ray_tracing
    integer,dimension(1:3) :: ray_traced_pos 
    integer,dimension(1:3) :: global_pos
    real(kind=dp) :: path
    real(kind=dp) :: number_density_atom
    real(kind=dp) :: xHI,xHeI,xHeII
    integer :: i_dim
       
    do i_dim = 1,3
      ray_traced_pos(i_dim) = srcpos(i_dim,source_index)
      global_pos(i_dim) = ray_traced_pos(i_dim)
    enddo

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

    coldens_in_HI(global_pos(1),global_pos(2),global_pos(3)) = 0.0
    coldens_in_HeI(global_pos(1),global_pos(2),global_pos(3)) = 0.0
    coldens_in_HeII(global_pos(1),global_pos(2),global_pos(3)) = 0.0
    path = 0.5*cell_size
    number_density_atom = number_density_array(global_pos(1),global_pos(2),global_pos(3))
    ring_vol(global_pos(1),global_pos(2),global_pos(3)) = cell_size*cell_size*cell_size 
    coldens_out_HI(global_pos(1),global_pos(2),global_pos(3)) = coldens(path,xHI,number_density_atom,(1.0_dp-abu_he))
    coldens_out_HeI(global_pos(1),global_pos(2),global_pos(3)) = coldens(path,xHeI,number_density_atom,abu_he)
    coldens_out_HeII(global_pos(1),global_pos(2),global_pos(3)) = coldens(path,xHeII,number_density_atom,abu_he)

  end subroutine column_density_0D_point

end module column_density_0D

