module ionization_equation

  use precision, only: dp
  use my_mpi 
  !use recombination_collision_factor, only: recombination_collision_factor_initialization
  use c2ray_parameters, only: minimum_fractional_change
  use c2ray_parameters, only: minimum_fraction_of_atoms
  use thermalevolution, only: thermal
  use tped, only: electrondens
  use doric_module, only: doric, prepare_doric_factors, coldens

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine solve_ionization_equation(dt, number_density_atom, &
        current_begin_xHI, current_begin_xHII, current_begin_xHeI, current_begin_xHeII, current_begin_xHeIII, &
        current_avg_xHI, current_avg_xHII, current_avg_xHeI, current_avg_xHeII, current_avg_xHeIII, &
        current_end_xHI, current_end_xHII, current_end_xHeI, current_end_xHeII, current_end_xHeIII, & 
        current_begin_temper,current_avg_temper,current_end_temper, &
        photoionization_HI_rate,photoionization_HeI_rate,photoionization_HeII_rate,heating_rate,pos)

    implicit none

    real(kind=dp),intent(in) :: dt
    real(kind=dp),intent(in) :: number_density_atom
    real(kind=dp),intent(inout) :: current_begin_xHI
    real(kind=dp),intent(inout) :: current_begin_xHII
    real(kind=dp),intent(inout) :: current_begin_xHeI
    real(kind=dp),intent(inout) :: current_begin_xHeII
    real(kind=dp),intent(inout) :: current_begin_xHeIII
    real(kind=dp),intent(inout) :: current_avg_xHI
    real(kind=dp),intent(inout) :: current_avg_xHII
    real(kind=dp),intent(inout) :: current_avg_xHeI
    real(kind=dp),intent(inout) :: current_avg_xHeII
    real(kind=dp),intent(inout) :: current_avg_xHeIII
    real(kind=dp),intent(inout) :: current_end_xHI
    real(kind=dp),intent(inout) :: current_end_xHII
    real(kind=dp),intent(inout) :: current_end_xHeI
    real(kind=dp),intent(inout) :: current_end_xHeII
    real(kind=dp),intent(inout) :: current_end_xHeIII
    real(kind=dp), intent(in) :: current_begin_temper
    real(kind=dp), intent(inout) :: current_avg_temper
    real(kind=dp), intent(inout) :: current_end_temper
    real(kind=dp), intent(in) :: photoionization_HI_rate
    real(kind=dp), intent(in) :: photoionization_HeI_rate
    real(kind=dp), intent(in) :: photoionization_HeII_rate
    real(kind=dp), intent(in) :: heating_rate
    integer,dimension(1:3),intent(in) :: pos
    real(kind=dp) :: prev_avg_xHI
    real(kind=dp) :: prev_avg_xHII
    real(kind=dp) :: prev_avg_xHeI
    real(kind=dp) :: prev_avg_xHeII
    real(kind=dp) :: prev_avg_xHeIII
    real(kind=dp) :: prev_temper
    real(kind=dp) :: number_density_electron
    real(kind=dp) :: yfrac 
    real(kind=dp) :: zfrac
    real(kind=dp) :: y2afrac 
    real(kind=dp) :: y2bfrac 
    real(kind=dp) :: intermediate_ion_avg_HeII
    real(kind=dp) :: intermediate_ion_avg_HeI
    real(kind=dp) :: intermediate_ion_avg_HI
    real(kind=dp) :: intermediate_ion_end_HI
    real(kind=dp) :: intermediate_ion_end_HII
    real(kind=dp) :: intermediate_ion_end_HeI
    real(kind=dp) :: intermediate_ion_end_HeII
    real(kind=dp) :: intermediate_ion_end_HeIII
    integer :: i_iteration


    i_iteration = 0

    do 

      i_iteration = i_iteration+1

      prev_temper = current_end_temper
      current_end_temper = current_begin_temper

      ! Save the values of yh_av found in the previous iteration
      prev_avg_xHI = current_avg_xHI
      prev_avg_xHII = current_avg_xHII
      prev_avg_xHeI = current_avg_xHeI
      prev_avg_xHeII = current_avg_xHeII
      prev_avg_xHeIII = current_avg_xHeIII
          
      ! Calculate (mean) electron density
      number_density_electron = electrondens(number_density_atom,current_avg_xHII,&
                        current_avg_xHeII,current_avg_xHeIII)

      !call recombination_collision_factor_initialization(current_avg_temper) 

      call prepare_doric_factors(current_end_xHI,current_end_xHeI,current_end_xHeII,&
                                 number_density_atom,yfrac,zfrac,y2afrac,y2bfrac)
       
      call doric(dt,number_density_electron,&
        current_begin_xHI, current_begin_xHII, current_begin_xHeI, current_begin_xHeII, current_begin_xHeIII, &
        current_avg_xHI, current_avg_xHII, current_avg_xHeI, current_avg_xHeII, current_avg_xHeIII, &
        current_end_xHI, current_end_xHII, current_end_xHeI, current_end_xHeII, current_end_xHeIII, & 
        current_avg_temper, photoionization_HI_rate,photoionization_HeI_rate,photoionization_HeII_rate, &
        yfrac,zfrac,y2afrac,y2bfrac)

      number_density_electron = electrondens(number_density_atom,current_avg_xHII,&
                        current_avg_xHeII,current_avg_xHeIII)

      call prepare_doric_factors(current_end_xHI,current_end_xHeI,current_end_xHeII,&
                                 number_density_atom,yfrac,zfrac,y2afrac,y2bfrac)

      ! Save results of first pass over doric
      intermediate_ion_end_HI = current_end_xHI
      intermediate_ion_end_HII = current_end_xHII
      intermediate_ion_end_HeI = current_end_xHeI
      intermediate_ion_end_HeII = current_end_xHeII
      intermediate_ion_end_HeIII = current_end_xHeIII
      intermediate_ion_avg_HI = current_avg_xHI
      intermediate_ion_avg_HeI = current_avg_xHeI
      intermediate_ion_avg_HeII = current_avg_xHeII 

      call doric(dt,number_density_electron,&
        current_begin_xHI, current_begin_xHII, current_begin_xHeI, current_begin_xHeII, current_begin_xHeIII, &
        current_avg_xHI, current_avg_xHII, current_avg_xHeI, current_avg_xHeII, current_avg_xHeIII, &
        current_end_xHI, current_end_xHII, current_end_xHeI, current_end_xHeII, current_end_xHeIII, & 
        current_avg_temper, photoionization_HI_rate,photoionization_HeI_rate,photoionization_HeII_rate, &
        yfrac,zfrac,y2afrac,y2bfrac)

      ! Average the answers from the two passes over doric
      current_end_xHI = (current_end_xHI+intermediate_ion_end_HI)/2.0_dp
      current_end_xHII = (current_end_xHII+intermediate_ion_end_HII)/2.0_dp
      current_end_xHeI = (current_end_xHeI+intermediate_ion_end_HeI)/2.0_dp
      current_end_xHeII = (current_end_xHeII+intermediate_ion_end_HeII)/2.0_dp
      current_end_xHeIII = (current_end_xHeIII+intermediate_ion_end_HeIII)/2.0_dp
      current_avg_xHI = (current_avg_xHI+intermediate_ion_avg_HI)/2.0_dp
      current_avg_xHeI = (current_avg_xHeI+intermediate_ion_avg_HeI)/2.0_dp
      current_avg_xHeII = (current_avg_xHeII+intermediate_ion_avg_HeII)/2.0_dp

      number_density_electron = electrondens(number_density_atom,current_avg_xHII,&
                        current_avg_xHeII,current_avg_xHeIII)

      call thermal(dt,current_end_temper,current_avg_temper,number_density_electron,number_density_atom,&
        current_begin_xHI, current_begin_xHII, current_begin_xHeI, current_begin_xHeII, current_begin_xHeIII, &
        current_avg_xHI, current_avg_xHII, current_avg_xHeI, current_avg_xHeII, current_avg_xHeIII, &
        current_end_xHI, current_end_xHII, current_end_xHeI, current_end_xHeII, current_end_xHeIII, & 
        heating_rate,pos)

       if ((abs((current_avg_xHI-prev_avg_xHI)/current_avg_xHI) .lt. &
            minimum_fractional_change .or. &
            (current_avg_xHI < minimum_fraction_of_atoms)).and. &         
            (abs((current_avg_xHeI-prev_avg_xHeI)/current_avg_xHeI) .lt. &
            minimum_fractional_change .or. &
            (current_avg_xHeI < minimum_fraction_of_atoms)).and. &
            (abs((current_avg_xHeIII-prev_avg_xHeIII)/current_avg_xHeIII) .lt. & 
            minimum_fractional_change .or. &
            (current_avg_xHeIII < minimum_fraction_of_atoms)) .and. &           
            (abs((current_end_temper-prev_temper)/current_end_temper) .lt. minimum_fractional_change)) then  
          exit
       endif
       
       ! Warn about non-convergence and terminate iteration
       if (i_iteration .gt. 500) then
          write(*,*) 'in ionization equation not converge'
          exit
       endif
    enddo

  end subroutine solve_ionization_equation
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ionization_equation
