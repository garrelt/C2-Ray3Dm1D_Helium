module local_evolution

  use precision, only: dp
  use my_mpi 
  use c2ray_parameters, only: minimum_fractional_change
  use c2ray_parameters, only: minimum_fraction_of_atoms
  use c2ray_parameters, only: epsilon_value
  use array, only: xHI_array,xHII_array
  use array, only: xHeI_array,xHeII_array,xHeIII_array
  use array, only: temperature_array
  use array, only: photoionization_HI_array
  use array, only: photoionization_HeI_array
  use array, only: photoionization_HeII_array
  use array, only: heating_array
  use array, only: number_density_array
  use array, only: intermediate_average_xHI
  use array, only: intermediate_average_xHII
  use array, only: intermediate_average_xHeI
  use array, only: intermediate_average_xHeII
  use array, only: intermediate_average_xHeIII
  use array, only: intermediate_end_temperature
  use array, only: LTE_MPI_intermediate_average_xHI
  use array, only: LTE_MPI_intermediate_average_xHII
  use array, only: LTE_MPI_intermediate_average_xHeI
  use array, only: LTE_MPI_intermediate_average_xHeII
  use array, only: LTE_MPI_intermediate_average_xHeIII
  use array, only: LTE_MPI_intermediate_average_temperature
  use array, only: LTE_MPI_intermediate_end_xHI
  use array, only: LTE_MPI_intermediate_end_xHII
  use array, only: LTE_MPI_intermediate_end_xHeI
  use array, only: LTE_MPI_intermediate_end_xHeII
  use array, only: LTE_MPI_intermediate_end_xHeIII
  use array, only: LTE_MPI_intermediate_end_temperature
  use array, only: ADP_MPI_intermediate_average_xHI
  use array, only: ADP_MPI_intermediate_average_xHII
  use array, only: ADP_MPI_intermediate_average_xHeI
  use array, only: ADP_MPI_intermediate_average_xHeII
  use array, only: ADP_MPI_intermediate_average_xHeIII
  use array, only: ADP_MPI_intermediate_average_temperature
  use array, only: ADP_MPI_intermediate_end_xHI
  use array, only: ADP_MPI_intermediate_end_xHII
  use array, only: ADP_MPI_intermediate_end_xHeI
  use array, only: ADP_MPI_intermediate_end_xHeII
  use array, only: ADP_MPI_intermediate_end_xHeIII
  use array, only: ADP_MPI_intermediate_end_temperature
  use ionization_equation, only: solve_ionization_equation

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine evolution_of_cell(dt,pos,convergence_failure,case_of_evolution)

    implicit none

    real(kind=dp),intent(in) :: dt 
    integer,dimension(1:3),intent(in) :: pos
    integer,intent(inout) :: convergence_failure
    character,intent(in) :: case_of_evolution
    real(kind=dp) :: number_density_atom 
    real(kind=dp) :: prev_avg_xHI
    real(kind=dp) :: prev_avg_xHII
    real(kind=dp) :: prev_avg_xHeI
    real(kind=dp) :: prev_avg_xHeII
    real(kind=dp) :: prev_avg_xHeIII
    real(kind=dp) :: prev_end_temper
    real(kind=dp) :: current_begin_xHI
    real(kind=dp) :: current_begin_xHII
    real(kind=dp) :: current_begin_xHeI
    real(kind=dp) :: current_begin_xHeII
    real(kind=dp) :: current_begin_xHeIII
    real(kind=dp) :: current_avg_xHI
    real(kind=dp) :: current_avg_xHII
    real(kind=dp) :: current_avg_xHeI
    real(kind=dp) :: current_avg_xHeII
    real(kind=dp) :: current_avg_xHeIII
    real(kind=dp) :: current_end_xHI
    real(kind=dp) :: current_end_xHII
    real(kind=dp) :: current_end_xHeI
    real(kind=dp) :: current_end_xHeII
    real(kind=dp) :: current_end_xHeIII
    real(kind=dp) :: current_begin_temper
    real(kind=dp) :: current_avg_temper
    real(kind=dp) :: current_end_temper
    real(kind=dp) :: photoionization_HI_rate
    real(kind=dp) :: photoionization_HeI_rate
    real(kind=dp) :: photoionization_HeII_rate
    real(kind=dp) :: heating_rate

    current_begin_xHI = max(epsilon_value,xHI_array(pos(1),pos(2),pos(3)))
    current_begin_xHII = max(epsilon_value,xHII_array(pos(1),pos(2),pos(3)))
    current_begin_xHeI = max(epsilon_value,xHeI_array(pos(1),pos(2),pos(3)))
    current_begin_xHeII = max(epsilon_value,xHeII_array(pos(1),pos(2),pos(3)))
    current_begin_xHeIII = max(epsilon_value,xHeIII_array(pos(1),pos(2),pos(3)))
    current_begin_temper = temperature_array(pos(1),pos(2),pos(3))

    select case (case_of_evolution)

    case ('L')

    current_avg_xHI = xHI_array(pos(1),pos(2),pos(3))
    current_avg_xHII = xHII_array(pos(1),pos(2),pos(3))
    current_avg_xHeI = xHeI_array(pos(1),pos(2),pos(3))
    current_avg_xHeII = xHeII_array(pos(1),pos(2),pos(3))
    current_avg_xHeIII = xHeIII_array(pos(1),pos(2),pos(3))
    current_avg_temper = temperature_array(pos(1),pos(2),pos(3))
    current_end_xHI = xHI_array(pos(1),pos(2),pos(3))
    current_end_xHII = xHII_array(pos(1),pos(2),pos(3))
    current_end_xHeI = xHeI_array(pos(1),pos(2),pos(3))
    current_end_xHeII = xHeII_array(pos(1),pos(2),pos(3))
    current_end_xHeIII = xHeIII_array(pos(1),pos(2),pos(3))
    current_end_temper = temperature_array(pos(1),pos(2),pos(3))

    case ('A')

    current_avg_xHI = ADP_MPI_intermediate_average_xHI(pos(1),pos(2),pos(3))
    current_avg_xHII = ADP_MPI_intermediate_average_xHII(pos(1),pos(2),pos(3))
    current_avg_xHeI = ADP_MPI_intermediate_average_xHeI(pos(1),pos(2),pos(3))
    current_avg_xHeII = ADP_MPI_intermediate_average_xHeII(pos(1),pos(2),pos(3))
    current_avg_xHeIII = ADP_MPI_intermediate_average_xHeIII(pos(1),pos(2),pos(3))
    current_avg_temper = ADP_MPI_intermediate_average_temperature(pos(1),pos(2),pos(3))
    current_end_xHI = ADP_MPI_intermediate_end_xHI(pos(1),pos(2),pos(3))
    current_end_xHII = ADP_MPI_intermediate_end_xHII(pos(1),pos(2),pos(3))
    current_end_xHeI = ADP_MPI_intermediate_end_xHeI(pos(1),pos(2),pos(3))
    current_end_xHeII = ADP_MPI_intermediate_end_xHeII(pos(1),pos(2),pos(3))
    current_end_xHeIII = ADP_MPI_intermediate_end_xHeIII(pos(1),pos(2),pos(3))
    current_end_temper = ADP_MPI_intermediate_end_temperature(pos(1),pos(2),pos(3))

    end select

    number_density_atom = number_density_array(pos(1),pos(2),pos(3))
    photoionization_HI_rate = photoionization_HI_array(pos(1),pos(2),pos(3))
    photoionization_HeI_rate = photoionization_HeI_array(pos(1),pos(2),pos(3))
    photoionization_HeII_rate = photoionization_HeII_array(pos(1),pos(2),pos(3))
    heating_rate = heating_array(pos(1),pos(2),pos(3))
!if (pos(1).eq.5 .and. pos(2).eq.5 .and. pos(3).eq.5 ) then
!write(*,*) 'local evolution case ',case_of_evolution, 'photo HI rate is ',photoionization_HI_rate
!endif

    call solve_ionization_equation(dt,number_density_atom,&
        current_begin_xHI, current_begin_xHII, current_begin_xHeI, current_begin_xHeII, current_begin_xHeIII, &
        current_avg_xHI, current_avg_xHII, current_avg_xHeI, current_avg_xHeII, current_avg_xHeIII, &
        current_end_xHI, current_end_xHII, current_end_xHeI, current_end_xHeII, current_end_xHeIII, & 
        current_begin_temper, current_avg_temper, current_end_temper, &
        photoionization_HI_rate,photoionization_HeI_rate,photoionization_HeII_rate,heating_rate,pos)

    select case (case_of_evolution)

    case ('L')

    LTE_MPI_intermediate_average_xHI(pos(1),pos(2),pos(3)) = current_avg_xHI
    LTE_MPI_intermediate_average_xHII(pos(1),pos(2),pos(3)) = current_avg_xHII
    LTE_MPI_intermediate_average_xHeI(pos(1),pos(2),pos(3)) = current_avg_xHeI
    LTE_MPI_intermediate_average_xHeII(pos(1),pos(2),pos(3)) = current_avg_xHeII
    LTE_MPI_intermediate_average_xHeIII(pos(1),pos(2),pos(3)) = current_avg_xHeIII
    LTE_MPI_intermediate_average_temperature(pos(1),pos(2),pos(3)) = current_avg_temper
    LTE_MPI_intermediate_end_xHI(pos(1),pos(2),pos(3)) = current_end_xHI
    LTE_MPI_intermediate_end_xHII(pos(1),pos(2),pos(3)) = current_end_xHII
    LTE_MPI_intermediate_end_xHeI(pos(1),pos(2),pos(3)) = current_end_xHeI
    LTE_MPI_intermediate_end_xHeII(pos(1),pos(2),pos(3)) = current_end_xHeII
    LTE_MPI_intermediate_end_xHeIII(pos(1),pos(2),pos(3)) = current_end_xHeIII
    LTE_MPI_intermediate_end_temperature(pos(1),pos(2),pos(3)) = current_end_temper

    case ('A')

    prev_avg_xHI = intermediate_average_xHI(pos(1),pos(2),pos(3))
    prev_avg_xHII = intermediate_average_xHII(pos(1),pos(2),pos(3))
    prev_avg_xHeI = intermediate_average_xHeI(pos(1),pos(2),pos(3))
    prev_avg_xHeII = intermediate_average_xHeII(pos(1),pos(2),pos(3))
    prev_avg_xHeIII = intermediate_average_xHeIII(pos(1),pos(2),pos(3))
    prev_end_temper = intermediate_end_temperature(pos(1),pos(2),pos(3))

    if ( (abs((current_avg_xHI-prev_avg_xHI)) .gt. minimum_fractional_change .and. &
          abs((current_avg_xHI-prev_avg_xHI)/current_avg_xHI) .gt. minimum_fractional_change .and. &
              (current_avg_xHI .gt. minimum_fraction_of_atoms)) .or. &
         (abs((current_avg_xHeI-prev_avg_xHeI)) .gt. minimum_fractional_change .and. &
          abs((current_avg_xHeI-prev_avg_xHeI)/current_avg_xHeI) .gt. minimum_fractional_change .and. &
              (current_avg_xHeI .gt. minimum_fraction_of_atoms)) .or. &
         (abs((current_avg_xHeIII-prev_avg_xHeIII)) .gt. minimum_fractional_change .and. &
          abs((current_avg_xHeIII-prev_avg_xHeIII)/current_avg_xHeIII) .gt. minimum_fractional_change .and. &
              (current_avg_xHeIII .gt. minimum_fraction_of_atoms)) .or. & 
         (abs((current_end_temper-prev_end_temper)/prev_end_temper) .gt. minimum_fractional_change)) then
       convergence_failure = convergence_failure+1

    endif

    ADP_MPI_intermediate_average_xHI(pos(1),pos(2),pos(3)) = current_avg_xHI
    ADP_MPI_intermediate_average_xHII(pos(1),pos(2),pos(3)) = current_avg_xHII
    ADP_MPI_intermediate_average_xHeI(pos(1),pos(2),pos(3)) = current_avg_xHeI 
    ADP_MPI_intermediate_average_xHeII(pos(1),pos(2),pos(3)) = current_avg_xHeII
    ADP_MPI_intermediate_average_xHeIII(pos(1),pos(2),pos(3)) = current_avg_xHeIII
    ADP_MPI_intermediate_average_temperature(pos(1),pos(2),pos(3)) = current_avg_temper
    ADP_MPI_intermediate_end_xHI(pos(1),pos(2),pos(3)) = current_end_xHI
    ADP_MPI_intermediate_end_xHII(pos(1),pos(2),pos(3)) = current_end_xHII
    ADP_MPI_intermediate_end_xHeI(pos(1),pos(2),pos(3)) = current_end_xHeI
    ADP_MPI_intermediate_end_xHeII(pos(1),pos(2),pos(3)) = current_end_xHeII
    ADP_MPI_intermediate_end_xHeIII(pos(1),pos(2),pos(3)) = current_end_xHeIII
    ADP_MPI_intermediate_end_temperature(pos(1),pos(2),pos(3)) = current_end_temper

    end select


  end subroutine evolution_of_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module local_evolution
