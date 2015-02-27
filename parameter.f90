module parameter

  use precision, only: dp

  implicit none

  ! size of a cell
  real(kind=dp) :: cell_size 

  ! volume of a cell
  real(kind=dp) :: cell_volume 

  ! cell_volume of entire simulation box 
  real(kind=dp) :: grid_volume

  ! Box size in Mpc/h comoving
  !real(kind=dp), parameter :: boxsize = 30
  !real(kind=dp), parameter :: boxsize = 10 ! Test 1
  !real(kind=dp), parameter :: boxsize = 30 ! Test 1.5
  real(kind=dp), parameter :: boxsize = 0.5 !Test 4

  ! size of the mesh for spatial coordinate.
  !integer, parameter :: mesh = 30
  integer, parameter :: mesh = 128 ! Test 4
  !integer, parameter :: mesh = 30 ! Test 1
  !integer, parameter :: mesh = 30 ! Test 1.5
  integer, parameter :: mesh3 = mesh*mesh*mesh

  !logical, parameter :: do_LTE = .true.
  logical, parameter :: do_LTE = .false.

  ! I-front takes this number of timestep to proceed
  integer, parameter :: time_partition = 10

  !> true if checking photonstatistics
  logical,parameter :: do_photonstatistics = .true.

  !> true if doing secondary ionization
  logical,parameter :: do_secondary = .true.
  !logical,parameter :: do_secondary = .false.

  !logical, parameter :: photo_conservation = .true.
  logical, parameter :: photo_conservation = .false.

  !> do thermal ! not yet implemented in heating equation
  logical, parameter :: do_thermal = .true.
  !logical, parameter :: do_thermal = .false.

  !> do periodic
  !logical, parameter :: do_periodic = .true.
  logical, parameter :: do_periodic = .false. ! Test 4

  ! coupled photo equation
  !logical, parameter :: do_couple = .true.
  logical, parameter :: do_couple = .false. 

  ! use constant timestep or adaptive timestep?
  !logical, parameter :: do_constant = .true.
  logical, parameter :: do_constant = .false.

  ! At time estimation method
  !logical, parameter :: do_AT_accurate_timestep = .true.
  logical, parameter :: do_AT_accurate_timestep =  .false.

!!!!!!!!!!!!!!!!!!!!!! radiation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer,parameter :: NumFreq = 512      ! Number of integration points in each of the frequency bins
  integer,parameter :: NumTau = 2000      ! Number of table points for the optical depth
  integer,parameter :: NumBndin1 = 1      ! Number of frequency sub-bins in interval 1 
  integer,parameter :: NumBndin2 = 26     ! Number of frequency sub-bins in interval 2
  integer,parameter :: NumBndin3 = 20     ! Number of frequency sub-bins in interval 3
  integer,parameter :: NumFreqBnd=NumBndin1+NumBndin2+NumBndin3       ! Total number of frequency bins
  integer,parameter :: NumheatBin=NumBndin1+NumBndin2*2+NumBndin3*3   ! Total number of heating bins 

  ! Optical depths at the entrance of the grid.
  ! It can be used if radiation enters the simulation volume from the outside.
  real(kind=dp) :: boundary_tauHI = 0.0
  real(kind=dp) :: boundary_tauHeI = 0.0
  real(kind=dp) :: boundary_tauHeII = 0.0

  ! Parameters defining the optical depth entries in the table.
  real(kind=dp),parameter :: minlogtau = -20.0                             ! Table position starts at log10(minlogtau)
  real(kind=dp),parameter :: maxlogtau = 4.0                               ! Table position ends at log10(maxlogtau) 
  real(kind=dp),parameter :: dlogtau = (maxlogtau-minlogtau)/real(NumTau)  ! dlogtau is the step size in log10(tau)

  ! Some boring variables  
  ! in general, I'm following Ricotti et al 2002
  real(kind=dp), dimension(1:3) :: CR1 = (/0.3908_dp, 0.0554_dp, 1.0_dp/)
  real(kind=dp), dimension(1:3) :: CR2 = (/0.6941_dp,0.0984_dp,3.9811_dp/)
  real(kind=dp), dimension(1:3) :: bR1 = (/0.4092_dp, 0.4614_dp, 0.2663_dp/)
  real(kind=dp), dimension(1:3) :: dR1 = (/1.7592_dp, 1.6660_dp, 1.3163_dp/)
  real(kind=dp), dimension(1:3) :: aR2 = (/0.2_dp,0.2_dp,0.4_dp/)
  real(kind=dp), dimension(1:3) :: bR2 = (/0.38_dp,0.38_dp,0.34_dp/)


  ! Stellar properties
  real(kind=dp) :: R_star_nominal       ! Black body radius
  real(kind=dp) :: L_star_nominal       ! Black body luminosity

  ! Power law source properties
  real(kind=dp) :: pl_minfreq           ! Minimum frequency for integration of total power
  real(kind=dp) :: pl_maxfreq           ! Maximum frequency for integration of total power
  real(kind=dp) :: pl_scaling           ! The scaling of the flux
  real(kind=dp) :: bb_source_ionzing_photon_rate  ! The rate of ionizing photon generated from the source
  real(kind=dp) :: pl_source_ionzing_photon_rate  ! The rate of ionizing photon generated from the source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=dp) :: input_source_temperature 
  real(kind=dp) :: input_IGM_temperature

  integer ::number_of_redshift
  real(kind=dp) :: total_simulation_time
  integer :: number_timesteps !< Number of time steps between redshift slices
  integer :: number_outputs !< Number of outputs between redshift slices
  integer, parameter ::  min_numproc_master_slave=10

  ! it records the average redshift at a time slice
  real(kind=dp) :: current_redshift 
  integer :: number_of_source

  ! choice of clumping model
  integer, parameter :: type_of_clumping = 1

  ! clumping factor if constant
  real, parameter :: constant_clumping_factor = 1.0

  integer :: iteration_counter
  integer :: convergence_failure

  integer :: LTE_iteration_counter
  integer :: LTE_convergence_failure


    integer :: conv_criterion
  real(kind=dp) :: actual_dt !< actual time step (s)

  real(kind=dp) :: output_time_dt !< time interval between outputs (s)

  real(kind=dp) :: end_time !< end time

  ! column density for stopping chemisty
  real(kind=dp),parameter :: max_coldensh = 1.998e21_dp

  ! number of cells needed to be evolved adaptively
  integer :: ADP_size_of_m

  ! number of cells needed to be evolved LTEly
  integer :: LTE_size_of_m 

  integer :: i_partition

  real(kind=dp) :: sim_time !< actual time (s)
  real(kind=dp) :: next_output_time !< time of next output (s)

  ! adaptive time step
  real(kind=dp) :: adaptive_dt

  ! assume the source is put at (0,0,0), the following are
  ! the max and min values for which a cell at (x,y,z), 
  ! it follows that
  ! length_neg_x <= x <= length_pos_x
  ! length_neg_y <= y <= length_pos_y
  ! length_neg_z <= z <= length_pos_z

  integer :: length_pos_x
  integer :: length_neg_x
  integer :: length_pos_y
  integer :: length_neg_y
  integer :: length_pos_z
  integer :: length_neg_z

  real :: real_length_pos_x
  real :: real_length_neg_x
  real :: real_length_pos_y
  real :: real_length_neg_y
  real :: real_length_pos_z
  real :: real_length_neg_z

  ! this number represents 1/2
  real :: half = 0.5
  logical :: LTE_non_evolved_cells_exist

end module parameter
