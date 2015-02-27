module photonstatistics
  
  use precision, only: dp
  use my_mpi
  use parameter, only: do_photonstatistics, number_of_source, next_output_time, sim_time, &
                       actual_dt, output_time_dt
  use array, only: LTE_photon_loss_array, ADP_photon_loss_array, global_LTE_array, &
                   global_evolved_LTE_array
  use sourceprops, only: BB_NormFlux, PL_NormFlux
  use c2ray_parameters, only: bb_S_star_nominal,pl_S_star_nominal
  use array, only: number_density_array, xHI_array, xHII_array, xHeI_array, xHeII_array, xHeIII_array
  use parameter, only: cell_volume, mesh
  use abundances, only: abu_h, abu_he
  use cgsconstants, only: temph0,temphe,colh0,colhe
  use array, only: intermediate_average_xHI, intermediate_average_xHII, intermediate_average_temperature, &
                   intermediate_average_xHeI, intermediate_average_xHeII, intermediate_average_xHeIII
  use clumping, only: clumping_factor
  use tped, only: electrondens

  implicit none

  ! Total number of recombinations
  real(kind=dp) :: LTE_recombination_number
  real(kind=dp) :: ADP_recombination_number
  ! Total number of collisional ionizations
  real(kind=dp) :: LTE_collisional_number
  real(kind=dp) :: ADP_collisional_number
  ! Total number of ionization due to recombination of Helium
  real(kind=dp) :: LTE_recombation_ionization_He
  real(kind=dp) :: ADP_recombation_ionization_He

  ! Grand total number of ionizing photons used
  real(kind=dp) :: grand_total_ionization_photons
  ! Grand total number of ionizing photons used
  real(kind=dp) :: grand_total_recombination_photons
  ! Grand total number of ionizing photons used
  real(kind=dp) :: grand_total_collisional_photons
  ! Grand total number of ionizing photons used
  real(kind=dp) :: grand_total_recombation_ionization_He
  ! Grand total number of ionizing photons produced
  real(kind=dp) :: grand_total_source_photons
  ! Grand total number of ionizing photons escaped
  real(kind=dp) :: grand_total_escaped_photons

  real(kind=dp) :: LTE_HI_before ! number of H atoms at start of time step
  real(kind=dp) :: LTE_HI_after ! number of H atoms at end of time step
  real(kind=dp) :: LTE_HII_before ! number of H ions at start of time step
  real(kind=dp) :: LTE_HII_after ! number of H ions at end of time step
  real(kind=dp) :: LTE_HeI_before ! number of He atoms at start of time step
  real(kind=dp) :: LTE_HeI_after ! number of He atoms at end of time step
  real(kind=dp) :: LTE_HeII_before ! number of He atoms at start of time step
  real(kind=dp) :: LTE_HeII_after ! number of He atoms at end of time step
  real(kind=dp) :: LTE_HeIII_before ! number of He atoms at start of time step
  real(kind=dp) :: LTE_HeIII_after ! number of He atoms at end of time step

  real(kind=dp) :: ADP_HI_before ! number of H atoms at start of time step
  real(kind=dp) :: ADP_HI_after ! number of H atoms at end of time step
  real(kind=dp) :: ADP_HII_before ! number of H ions at start of time step
  real(kind=dp) :: ADP_HII_after ! number of H ions at end of time step
  real(kind=dp) :: ADP_HeI_before ! number of He atoms at start of time step
  real(kind=dp) :: ADP_HeI_after ! number of He atoms at end of time step
  real(kind=dp) :: ADP_HeII_before ! number of He atoms at start of time step
  real(kind=dp) :: ADP_HeII_after ! number of He atoms at end of time step
  real(kind=dp) :: ADP_HeIII_before ! number of He atoms at start of time step
  real(kind=dp) :: ADP_HeIII_after ! number of He atoms at end of time step

  integer :: i ! mesh loop index (x)
  integer :: j ! mesh loop index (y)
  integer :: k ! mesh loop index (z)

  real(kind=dp) :: LTE_photon_number_in_ionization
  real(kind=dp) :: ADP_photon_number_in_ionization
  real(kind=dp) :: LTE_total_escaped_photons
  real(kind=dp) :: ADP_total_escaped_photons

  real(kind=dp) :: LTE_background
  real(kind=dp) :: ADP_background
  integer :: LTE_number_of_cell
  integer :: ADP_number_of_cell

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize the photon statistics
  subroutine photonstatistics_initialization()

    implicit none

    ! set total number of ionizing photons used to zero
    grand_total_ionization_photons = 0.0
    grand_total_recombination_photons = 0.0
    grand_total_collisional_photons = 0.0
    grand_total_recombation_ionization_He = 0.0
    grand_total_source_photons = 0.0
    grand_total_escaped_photons = 0.0

  end subroutine photonstatistics_initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize the photon statistics
  subroutine PS_array_initialization()

    implicit none

    allocate(LTE_photon_loss_array(1:number_of_source,0:nthreads-1))
    allocate(ADP_photon_loss_array(1:number_of_source,0:nthreads-1))
    LTE_photon_loss_array = 0.0
    ADP_photon_loss_array = 0.0
    LTE_background = 0.0
    ADP_background = 0.0

  end subroutine PS_array_initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine PS_array_destruction()

    implicit none

    deallocate(LTE_photon_loss_array)
    deallocate(ADP_photon_loss_array)

  end subroutine PS_array_destruction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine PS_photon_loss_total(case_of_evolution)

    implicit none

    character,intent(in) :: case_of_evolution
    real(kind=dp), dimension(1:number_of_source,0:nthreads-1) :: buffer
    real(kind=dp) :: delta_t
    integer :: mympierror

    select case (case_of_evolution)

    case ('L')

      LTE_total_escaped_photons = 0.0
      delta_t = next_output_time-sim_time
      call MPI_ALLREDUCE(LTE_photon_loss_array, buffer, number_of_source*nthreads, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_NEW, mympierror)        
      LTE_photon_loss_array(:,:) = buffer(:,:)       

      LTE_total_escaped_photons = sum(LTE_photon_loss_array)*delta_t
      LTE_background = LTE_background+sum(LTE_photon_loss_array)

      LTE_number_of_cell = count(global_LTE_array .eq. 1)-count(global_evolved_LTE_array .eq. 1)

      if (LTE_number_of_cell.eq.0) then
        BB_NormFlux(0) = 0.0
      else
        BB_NormFlux(0) = LTE_background/(LTE_number_of_cell*bb_S_star_nominal)
      endif

    case ('A')
 
      ADP_total_escaped_photons = 0.0 
      delta_t = actual_dt
      call MPI_ALLREDUCE(ADP_photon_loss_array, buffer, number_of_source*nthreads, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_NEW, mympierror)        
      ADP_photon_loss_array(:,:) = buffer(:,:)       
!write(*,*) ADP_photon_loss_array
      ADP_total_escaped_photons = sum(ADP_photon_loss_array)*delta_t 
      ADP_background = sum(ADP_photon_loss_array)

      ADP_number_of_cell = count(global_LTE_array .eq. 0)

      if (ADP_number_of_cell.eq.0) then
        BB_NormFlux(0) = 0.0
      else
        BB_NormFlux(0) = (ADP_background+LTE_background)/(ADP_number_of_cell*bb_S_star_nominal)
!write(*,*) BB_NormFlux(0)
      endif

    end select
 
  end subroutine PS_photon_loss_total

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine PS_reset()

    implicit none

    LTE_background = 0.0

  end subroutine PS_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine PS_ionization_state_before(case_of_evolution)

    implicit none

    character,intent(in) :: case_of_evolution

    select case (case_of_evolution)

    case ('L')

      ! Photon statistics: calculate the number of neutrals before integration
      LTE_HI_before = 0.0
      LTE_HII_before = 0.0
      LTE_HeI_before = 0.0
      LTE_HeII_before = 0.0
      LTE_HeIII_before = 0.0

      do k = 1,mesh
         do j = 1,mesh
            do i = 1,mesh
              if (global_LTE_array(i,j,k).eq.1 .and. &
                  global_evolved_LTE_array(i,j,k).eq.0) then

                LTE_HI_before = LTE_HI_before+number_density_array(i,j,k)*abu_h*xHI_array(i,j,k)
                LTE_HII_before = LTE_HII_before+number_density_array(i,j,k)*abu_h*xHII_array(i,j,k)
                LTE_HeI_before = LTE_HeI_before+number_density_array(i,j,k)*abu_he*xHeI_array(i,j,k)
                LTE_HeII_before = LTE_HeII_before+number_density_array(i,j,k)*abu_he*xHeII_array(i,j,k)
                LTE_HeIII_before = LTE_HeIII_before+number_density_array(i,j,k)*abu_he*xHeIII_array(i,j,k)

              endif
            enddo
         enddo
      enddo

      LTE_HI_before = LTE_HI_before*cell_volume
      LTE_HII_before = LTE_HII_before*cell_volume
      LTE_HeI_before = LTE_HeI_before*cell_volume
      LTE_HeII_before = LTE_HeII_before*cell_volume
      LTE_HeIII_before = LTE_HeIII_before*cell_volume

    case ('A')

      ! Photon statistics: calculate the number of neutrals before integration
      ADP_HI_before = 0.0
      ADP_HII_before = 0.0
      ADP_HeI_before = 0.0
      ADP_HeII_before = 0.0
      ADP_HeIII_before = 0.0

      do k = 1,mesh
         do j = 1,mesh
            do i = 1,mesh
              if (global_LTE_array(i,j,k).eq.0) then

                ADP_HI_before = ADP_HI_before+number_density_array(i,j,k)*xHI_array(i,j,k)
                ADP_HII_before = ADP_HII_before+number_density_array(i,j,k)*xHII_array(i,j,k)
                ADP_HeI_before = ADP_HeI_before+number_density_array(i,j,k)*xHeI_array(i,j,k)
                ADP_HeII_before = ADP_HeII_before+number_density_array(i,j,k)*xHeII_array(i,j,k)
                ADP_HeIII_before = ADP_HeIII_before+number_density_array(i,j,k)*xHeIII_array(i,j,k)

              endif
            enddo
         enddo
      enddo

      ADP_HI_before = ADP_HI_before*cell_volume*abu_h
      ADP_HII_before = ADP_HII_before*cell_volume*abu_h
      ADP_HeI_before = ADP_HeI_before*cell_volume*abu_he
      ADP_HeII_before = ADP_HeII_before*cell_volume*abu_he
      ADP_HeIII_before = ADP_HeIII_before*cell_volume*abu_he

    end select

  end subroutine PS_ionization_state_before

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine PS_ionization_state_after(case_of_evolution)

    implicit none

    character,intent(in) :: case_of_evolution

    select case (case_of_evolution)

    case ('L')

      ! Photon statistics: calculate the number of neutrals after integration
      LTE_HI_after = 0.0
      LTE_HII_after = 0.0
      LTE_HeI_after = 0.0
      LTE_HeII_after = 0.0
      LTE_HeIII_after = 0.0

      do k = 1,mesh
         do j = 1,mesh
            do i = 1,mesh
              if (global_LTE_array(i,j,k).eq.1 .and. &
                  global_evolved_LTE_array(i,j,k).eq.0) then

                LTE_HI_after = LTE_HI_after+number_density_array(i,j,k)*xHI_array(i,j,k)
                LTE_HII_after = LTE_HII_after+number_density_array(i,j,k)*xHII_array(i,j,k)
                LTE_HeI_after = LTE_HeI_after+number_density_array(i,j,k)*xHeI_array(i,j,k)
                LTE_HeII_after = LTE_HeII_after+number_density_array(i,j,k)*xHeII_array(i,j,k)
                LTE_HeIII_after = LTE_HeIII_after+number_density_array(i,j,k)*xHeIII_array(i,j,k)

              endif
            enddo
         enddo
      enddo

      LTE_HI_after = LTE_HI_after*cell_volume*abu_h
      LTE_HII_after = LTE_HII_after*cell_volume*abu_h
      LTE_HeI_after = LTE_HeI_after*cell_volume*abu_he
      LTE_HeII_after = LTE_HeII_after*cell_volume*abu_he
      LTE_HeIII_after = LTE_HeIII_after*cell_volume*abu_he

    case ('A')

      ! Photon statistics: calculate the number of neutrals after integration
      ADP_HI_after = 0.0
      ADP_HII_after = 0.0
      ADP_HeI_after = 0.0
      ADP_HeII_after = 0.0
      ADP_HeIII_after = 0.0

      do k = 1,mesh
         do j = 1,mesh
            do i = 1,mesh
              if (global_LTE_array(i,j,k).eq.0) then

                ADP_HI_after = ADP_HI_after+number_density_array(i,j,k)*abu_h*xHI_array(i,j,k)
                ADP_HII_after = ADP_HII_after+number_density_array(i,j,k)*abu_h*xHII_array(i,j,k)
                ADP_HeI_after = ADP_HeI_after+number_density_array(i,j,k)*abu_he*xHeI_array(i,j,k)
                ADP_HeII_after = ADP_HeII_after+number_density_array(i,j,k)*abu_he*xHeII_array(i,j,k)
                ADP_HeIII_after = ADP_HeIII_after+number_density_array(i,j,k)*abu_he*xHeIII_array(i,j,k)

              endif
            enddo
         enddo
      enddo

      ADP_HI_after = ADP_HI_after*cell_volume
      ADP_HII_after = ADP_HII_after*cell_volume
      ADP_HeI_after = ADP_HeI_after*cell_volume
      ADP_HeII_after = ADP_HeII_after*cell_volume
      ADP_HeIII_after = ADP_HeIII_after*cell_volume

    end select

  end subroutine PS_ionization_state_after

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine PS_recombination_collisional(case_of_evolution)

    implicit none

    character,intent(in) :: case_of_evolution

    real(kind=dp) :: number_density ! needed because ndens may be single precision
    real(kind=dp) :: average_xHI,average_xHII,average_temperature
    real(kind=dp) :: average_xHeI,average_xHeII,average_xHeIII
    real(kind=dp) :: delta_t
    real(kind=dp) :: sqrtt0, lambda_HI,lambda_HeI,lambda_HeII,dielectronic
    real(kind=dp) :: brech0,breche0,breche1,colli_HI,colli_HeI,colli_HeII


    select case (case_of_evolution)

    case ('L')

      LTE_recombination_number = 0.0
      LTE_collisional_number = 0.0
      LTE_recombation_ionization_He = 0.0
      delta_t = next_output_time-sim_time

      do k = 1,mesh 
        do j = 1,mesh
          do i = 1,mesh
           
            average_xHI = intermediate_average_xHI(i,j,k)
            average_xHII = intermediate_average_xHII(i,j,k)
            average_xHeI = intermediate_average_xHeI(i,j,k)
            average_xHeII = intermediate_average_xHeII(i,j,k)
            average_xHeIII = intermediate_average_xHeIII(i,j,k)
            average_temperature = intermediate_average_temperature(i,j,k)
            number_density = number_density_array(i,j,k)

            if (global_LTE_array(i,j,k).eq.1 .and. &
                global_evolved_LTE_array(i,j,k).eq.0) then

            lambda_HI =  2.0_dp*(temph0/average_temperature)                                             
            brech0=2.753e-14*lambda_HI**1.500_dp/(1.0_dp+(lambda_HI/2.740)**0.407)**2.242  

            if (average_temperature < 9.e3_dp) then
            lambda_HeI=2.0_dp*(temph0/average_temperature)
            breche0=2.753e-14_dp*lambda_HeI**1.500_dp/(1.0_dp+(lambda_HeI/2.740)**0.407)**2.242
            else
            lambda_HeI=   2.0_dp*(temphe(0)/average_temperature)
            dielectronic = 1.9e-3_dp*average_temperature**(-1.5_dp)* &
                           exp(-4.7e5_dp/average_temperature)*(1.0_dp+0.3_dp*exp(-9.4e4_dp/average_temperature))
            breche0= 1.260e-14_dp*lambda_HeI**0.750_dp + dielectronic
            endif

            lambda_HeII=2.0_dp*(temphe(1)/average_temperature)                                                  
            breche1=5.5060e-14_dp*lambda_HeII**1.5_dp/(1.0_dp+(lambda_HeII/2.740_dp)**0.407_dp)**2.242_dp 

            LTE_recombination_number = LTE_recombination_number+number_density* &
                       electrondens(number_density,average_xHII,average_xHeII,average_xHeIII)* &
                       clumping_factor* &
                       (average_xHII*brech0*abu_h+average_xHeII*breche0*abu_he*0.04_dp)

            sqrtt0 =sqrt(average_temperature)
            colli_HI =colh0*sqrtt0*exp(-temph0/average_temperature)
            colli_HeI=colhe(0)*sqrtt0*exp(-temphe(0)/average_temperature)
            colli_HeII=colhe(1)*sqrtt0*exp(-temphe(1)/average_temperature)

            LTE_collisional_number = LTE_collisional_number+number_density* &
                 electrondens(number_density,average_xHII,average_xHeII,average_xHeIII)* &
                (average_xHI*colli_HI+average_xHeI*colli_HeI+average_xHeII*colli_HeII)

            LTE_recombation_ionization_He = LTE_recombation_ionization_He + &
            number_density*abu_he*clumping_factor*electrondens(number_density,average_xHII,average_xHeII,average_xHeIII) * &
            (average_xHeII*breche0*0.96_dp+average_xHeIII*breche1*1.121_dp)

            endif

          enddo
        enddo
      enddo

      LTE_recombination_number = LTE_recombination_number*cell_volume*delta_t
      LTE_collisional_number = LTE_collisional_number*cell_volume*delta_t
      LTE_recombation_ionization_He = LTE_recombation_ionization_He*cell_volume*delta_t

    case ('A')

      ADP_recombination_number = 0.0
      ADP_collisional_number = 0.0
      delta_t = actual_dt

      do k = 1,mesh 
        do j = 1,mesh
          do i = 1,mesh

            average_xHI = intermediate_average_xHI(i,j,k)
            average_xHII = intermediate_average_xHII(i,j,k)
            average_xHeI = intermediate_average_xHeI(i,j,k)
            average_xHeII = intermediate_average_xHeII(i,j,k)
            average_xHeIII = intermediate_average_xHeIII(i,j,k)
            average_temperature = intermediate_average_temperature(i,j,k)
            number_density = number_density_array(i,j,k)

            if (global_LTE_array(i,j,k).eq.0) then

            lambda_HI =  2.0_dp*(temph0/average_temperature)                                             
            brech0=2.753e-14*lambda_HI**1.500_dp/(1.0_dp+(lambda_HI/2.740)**0.407)**2.242  

            if (average_temperature < 9.e3_dp) then
            lambda_HeI=2.0_dp*(temph0/average_temperature)
            breche0=2.753e-14_dp*lambda_HeI**1.500_dp/(1.0_dp+(lambda_HeI/2.740)**0.407)**2.242
            else
            lambda_HeI=   2.0_dp*(temphe(0)/average_temperature)
            dielectronic = 1.9e-3_dp*average_temperature**(-1.5_dp)* &
                           exp(-4.7e5_dp/average_temperature)*(1.0_dp+0.3_dp*exp(-9.4e4_dp/average_temperature))
            breche0= 1.260e-14_dp*lambda_HeI**0.750_dp + dielectronic
            endif

            lambda_HeII=2.0_dp*(temphe(1)/average_temperature)                                                  
            breche1=5.5060e-14_dp*lambda_HeII**1.5_dp/(1.0_dp+(lambda_HeII/2.740_dp)**0.407_dp)**2.242_dp 

            ADP_recombination_number = ADP_recombination_number+number_density* &
                       electrondens(number_density,average_xHII,average_xHeII,average_xHeIII)* &
                       clumping_factor* &
                       (average_xHII*brech0*abu_h+average_xHeII*breche0*abu_he*0.04_dp)

            sqrtt0 =sqrt(average_temperature)
            colli_HI =colh0*sqrtt0*exp(-temph0/average_temperature)
            colli_HeI=colhe(0)*sqrtt0*exp(-temphe(0)/average_temperature)
            colli_HeII=colhe(1)*sqrtt0*exp(-temphe(1)/average_temperature)

            ADP_collisional_number = ADP_collisional_number+number_density* &
                 electrondens(number_density,average_xHII,average_xHeII,average_xHeIII)* &
                (average_xHI*colli_HI+average_xHeI*colli_HeI+average_xHeII*colli_HeII)

            ADP_recombation_ionization_He = ADP_recombation_ionization_He + &
            number_density*abu_he*clumping_factor*electrondens(number_density,average_xHII,average_xHeII,average_xHeIII) * &
            (average_xHeII*breche0*0.96_dp+average_xHeIII*breche1*1.121_dp)

            endif

          enddo
        enddo
      enddo

      ADP_recombination_number = ADP_recombination_number*cell_volume*delta_t
      ADP_collisional_number = ADP_collisional_number*cell_volume*delta_t
      ADP_recombation_ionization_He = ADP_recombation_ionization_He*cell_volume*delta_t

    end select

  end subroutine PS_recombination_collisional
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine PS_summary(case_of_evolution)
   
    implicit none

    character,intent(in) :: case_of_evolution

    select case (case_of_evolution)

    case ('L')

      LTE_photon_number_in_ionization = (LTE_HI_before-LTE_HI_after)+(LTE_HeI_before-LTE_HeI_after)+&
                                        (LTE_HeIII_after-LTE_HeIII_before)
      grand_total_ionization_photons = grand_total_ionization_photons+LTE_photon_number_in_ionization
      grand_total_recombination_photons = grand_total_recombination_photons+LTE_recombination_number
      grand_total_collisional_photons = grand_total_collisional_photons+LTE_collisional_number
      grand_total_recombation_ionization_He = grand_total_recombation_ionization_He+LTE_recombation_ionization_He
      grand_total_escaped_photons = grand_total_escaped_photons+LTE_total_escaped_photons

    case ('A')

      ADP_photon_number_in_ionization = (ADP_HI_before-ADP_HI_after)+(ADP_HeI_before-ADP_HeI_after)+&
                                        (ADP_HeIII_after-ADP_HeIII_before)
      grand_total_ionization_photons = grand_total_ionization_photons+ADP_photon_number_in_ionization
      grand_total_recombination_photons = grand_total_recombination_photons+ADP_recombination_number
      grand_total_collisional_photons = grand_total_collisional_photons+ADP_collisional_number
      grand_total_recombation_ionization_He = grand_total_recombation_ionization_He+ADP_recombation_ionization_He
      grand_total_escaped_photons = grand_total_escaped_photons+ADP_total_escaped_photons

    end select
 
  end subroutine PS_summary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine PS_update_source_photons()

    implicit none

    grand_total_source_photons = grand_total_source_photons+&
                                 (sum(BB_NormFlux(1:number_of_source))*bb_S_star_nominal+&
                                  sum(PL_NormFlux(1:number_of_source))*pl_S_star_nominal)* &
                                 output_time_dt

  end subroutine PS_update_source_photons

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module photonstatistics
