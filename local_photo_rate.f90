module local_photo_rate

  use precision, only: dp 
  use my_mpi
  use parameter, only: max_coldensh, cell_volume, photo_conservation
  use array, only: number_density_array, &
                   LTE_photon_loss_array, ADP_photon_loss_array
  use sourceprops, only: srcpos
  use array, only: photoionization_HI_array
  use array, only: photoionization_HeI_array
  use array, only: photoionization_HeII_array
  use array, only: heating_array
  use array, only: coldens_in_HI
  use array, only: coldens_out_HI
  use array, only: coldens_in_HeI
  use array, only: coldens_out_HeI
  use array, only: coldens_in_HeII
  use array, only: coldens_out_HeII
  use array, only: ring_vol
  use abundances, only: abu_he
  use photo_rate, only: photoion_rates
  use type_definition, only: photrates, tablepos
  use OMP_LIB, only: omp_get_thread_num
#ifdef MY_OPENMP
  use OMP_LIB, only: omp_get_thread_num
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine photo_rate_computation(case_of_evolution,global_pos,source_index,xHI,xHII,xHeI,xHeII,boundary)

    implicit none

    character,intent(in) :: case_of_evolution
    integer, dimension(1:3), intent(in) :: global_pos
    integer, intent(in) :: source_index
    real(kind=dp), intent(in) :: xHI
    real(kind=dp), intent(in) :: xHII
    real(kind=dp), intent(in) :: xHeI
    real(kind=dp), intent(in) :: xHeII
    logical, intent(in) :: boundary
    type(photrates) :: phi
    real(kind=dp) :: ring_volume
    real(kind=dp) :: number_density_atom
    real(kind=dp) :: coldens_in_HI_sosad
    real(kind=dp) :: coldens_out_HI_sosad
    real(kind=dp) :: coldens_in_HeI_sosad
    real(kind=dp) :: coldens_out_HeI_sosad
    real(kind=dp) :: coldens_in_HeII_sosad
    real(kind=dp) :: coldens_out_HeII_sosad
    real(kind=dp) :: escaped_photon_rate
    integer :: thread_id


    coldens_in_HI_sosad = coldens_in_HI(global_pos(1),global_pos(2),global_pos(3))
    coldens_in_HeI_sosad = coldens_in_HeI(global_pos(1),global_pos(2),global_pos(3))
    coldens_in_HeII_sosad = coldens_in_HeII(global_pos(1),global_pos(2),global_pos(3))
    coldens_out_HI_sosad = coldens_out_HI(global_pos(1),global_pos(2),global_pos(3))
    coldens_out_HeI_sosad = coldens_out_HeI(global_pos(1),global_pos(2),global_pos(3))
    coldens_out_HeII_sosad = coldens_out_HeII(global_pos(1),global_pos(2),global_pos(3))
    number_density_atom = number_density_array(global_pos(1),global_pos(2),global_pos(3))  
    ring_volume = ring_vol(global_pos(1),global_pos(2),global_pos(3))

    if (coldens_in_HI_sosad .lt. max_coldensh .or. .true.) then 

      call photoion_rates(phi,coldens_in_HI_sosad,&
                          coldens_out_HI_sosad, &
                          coldens_in_HeI_sosad,&
                          coldens_out_HeI_sosad, &
                          coldens_in_HeII_sosad, &
                          coldens_out_HeII_sosad, &
                          ring_volume,source_index,xHII)


      phi%photo_cell_HI=phi%photo_cell_HI/(xHI*number_density_atom*(1.0_dp-abu_he))
      phi%photo_cell_HeI=phi%photo_cell_HeI/(xHeI*number_density_atom*abu_he)
      phi%photo_cell_HeII=phi%photo_cell_HeII/(xHeII*number_density_atom*abu_he)
      escaped_photon_rate = phi%photo_out

    else

      phi%photo_cell_HI = 0.0_dp
      phi%photo_cell_HeI = 0.0_dp
      phi%photo_cell_HeII = 0.0_dp
      phi%photo_out_HI = 0.0_dp
      phi%photo_out_HeI = 0.0_dp
      phi%photo_out_HeII = 0.0_dp
      phi%heat = 0.0_dp
      phi%photo_in = 0.0_dp
      phi%photo_out = 0.0_dp
      escaped_photon_rate = phi%photo_out

    endif


       photoionization_HI_array(global_pos(1),global_pos(2),global_pos(3))= &
            photoionization_HI_array(global_pos(1),global_pos(2),global_pos(3))+phi%photo_cell_HI
       photoionization_HeI_array(global_pos(1),global_pos(2),global_pos(3))=&
             photoionization_HeI_array(global_pos(1),global_pos(2),global_pos(3))+phi%photo_cell_HeI
       photoionization_HeII_array(global_pos(1),global_pos(2),global_pos(3))=&
             photoionization_HeII_array(global_pos(1),global_pos(2),global_pos(3))+phi%photo_cell_HeII   
       heating_array(global_pos(1),global_pos(2),global_pos(3))=&
             heating_array(global_pos(1),global_pos(2),global_pos(3))+phi%heat

!if (global_pos(1).eq.5 .and. global_pos(2).eq.5 .and. global_pos(3).eq.5) then
!write(*,*) 'photoionization_HI_array is ',photoionization_HI_array(global_pos(1),global_pos(2),global_pos(3))
!endif

    if (photo_conservation) then

      ! inclusion of background radiation
      call photoion_rates(phi,0.0_dp,coldens_out_HI_sosad-coldens_in_HI_sosad, &
                              0.0_dp,coldens_out_HeI_sosad-coldens_in_HeI_sosad, &
                              0.0_dp,coldens_out_HeII_sosad-coldens_in_HeII_sosad, &
                              cell_volume,0,xHII)
      phi%photo_cell_HI=phi%photo_cell_HI/(xHI*number_density_atom*(1.0_dp-abu_he))
      phi%photo_cell_HeI=phi%photo_cell_HeI/(xHeI*number_density_atom*abu_he)
      phi%photo_cell_HeII=phi%photo_cell_HeII/(xHeII*number_density_atom*abu_he)

       photoionization_HI_array(global_pos(1),global_pos(2),global_pos(3))= &
            photoionization_HI_array(global_pos(1),global_pos(2),global_pos(3))+phi%photo_cell_HI
       photoionization_HeI_array(global_pos(1),global_pos(2),global_pos(3))=&
             photoionization_HeI_array(global_pos(1),global_pos(2),global_pos(3))+phi%photo_cell_HeI
       photoionization_HeII_array(global_pos(1),global_pos(2),global_pos(3))=&
             photoionization_HeII_array(global_pos(1),global_pos(2),global_pos(3))+phi%photo_cell_HeII  

    endif

    thread_id = omp_get_thread_num()

    if (boundary) then

      select case (case_of_evolution)

      case ('L')

        LTE_photon_loss_array(source_index,thread_id) = &
           LTE_photon_loss_array(source_index,thread_id) + escaped_photon_rate*cell_volume/ring_volume

      case ('A')

        ADP_photon_loss_array(source_index,thread_id) = &
           ADP_photon_loss_array(source_index,thread_id) + escaped_photon_rate*cell_volume/ring_volume

      end select

    endif

  end subroutine photo_rate_computation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module local_photo_rate
