module redshift

  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use astroconstants, only: YEAR
  use array, only: redshift_array
  use parameter, only: number_of_redshift, total_simulation_time
  use cosmology_parameters, only: Omega0,H0

  implicit none

  real(kind=dp) :: initial_redshift

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Construct redshift sequence
  subroutine redshift_initialization()

    implicit none

    real(kind=dp) :: t0,timestep
    integer :: nz 

    if (rank .eq. 0) then
      write(logf,*) "Beginning of redshift initialization" 
    endif

    ! Set the number of redshift slices
    !number_of_redshift =  5
    number_of_redshift =  9 ! Test 1, Test 4
    !number_of_redshift =  9 ! Test 1.5
    allocate(redshift_array(1:number_of_redshift))

    ! Time step
    !timestep = 1e7*YEAR
    !timestep = 1e7*YEAR ! Test 1
    !timestep = 1e7*YEAR ! Test 1.5
    timestep = 5e4*YEAR ! Test 4

    total_simulation_time = timestep*(number_of_redshift-1)    

    ! Starting redshift
    redshift_array(1) = 9.0 ! Test 1
    !redshift_array(1) = 6.0 ! Test 1.5

    initial_redshift = redshift_array(1)

    if (rank .eq. 0) then
      write(logf,*) "Initial redshift = ",initial_redshift
      write(logf,*) "Number of redshift slice = ",number_of_redshift
      write(logf,*) "Time between redshift slices = ",timestep/YEAR,"years"
    endif

    ! Cosmological time corresponding to (initial) redshift
    ! NOTE: Good only for high-z!!!
    t0 = 2.*(1.+redshift_array(1))**(-1.5)/(3.*H0*sqrt(Omega0))
    do nz = 2,number_of_redshift
      redshift_array(nz) = -1+(1.+redshift_array(1))*(t0/(t0+real(nz-1)*timestep))**(2./3.)
    enddo

    if (rank .eq. 0) then
      write(logf,*) "End of redshift initialization" 
      write(logf,*)      
      flush(logf) 
    endif

  end subroutine redshift_initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module redshift
