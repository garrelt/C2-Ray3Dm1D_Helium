module array

  use precision, only: si, dp
  use parameter, only: time_partition, number_of_source,mesh



  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! material !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

  !real(kind=dp),dimension(:,:,:),allocatable :: number_density_array
  real,dimension(:,:,:),allocatable :: number_density_array ! Test 4
  real(kind=dp),dimension(:,:,:),allocatable :: temperature_array
  real(kind=dp),dimension(:,:,:),allocatable :: xHI_array
  real(kind=dp),dimension(:,:,:),allocatable :: xHII_array
  real(kind=dp),dimension(:,:,:),allocatable :: xHeI_array
  real(kind=dp),dimension(:,:,:),allocatable :: xHeII_array
  real(kind=dp),dimension(:,:,:),allocatable :: xHeIII_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Evolution !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Time-averaged H,He ionization fraction during evolution
  real(kind=dp),dimension(:,:,:),allocatable :: intermediate_average_xHI
  real(kind=dp),dimension(:,:,:),allocatable :: intermediate_average_xHII
  real(kind=dp),dimension(:,:,:),allocatable :: intermediate_average_xHeI
  real(kind=dp),dimension(:,:,:),allocatable :: intermediate_average_xHeII
  real(kind=dp),dimension(:,:,:),allocatable :: intermediate_average_xHeIII

  ! Intermediate result for H,He ionization fraction during evolution
  real(kind=dp),dimension(:,:,:),allocatable,public :: intermediate_end_xHI
  real(kind=dp),dimension(:,:,:),allocatable,public :: intermediate_end_xHII
  real(kind=dp),dimension(:,:,:),allocatable,public :: intermediate_end_xHeI
  real(kind=dp),dimension(:,:,:),allocatable,public :: intermediate_end_xHeII
  real(kind=dp),dimension(:,:,:),allocatable,public :: intermediate_end_xHeIII

  ! Intermediate result for Temperature during evolution
  real(kind=dp),dimension(:,:,:),allocatable :: intermediate_end_temperature
  real(kind=dp),dimension(:,:,:),allocatable :: intermediate_average_temperature

  real(kind=dp),dimension(:,:,:),allocatable :: photoionization_HI_array
  real(kind=dp),dimension(:,:,:),allocatable :: photoionization_HeI_array
  real(kind=dp),dimension(:,:,:),allocatable :: photoionization_HeII_array
  real(kind=dp),dimension(:,:,:),allocatable :: heating_array

  real(kind=dp),dimension(:,:,:),allocatable :: coldens_in_HI
  real(kind=dp),dimension(:,:,:),allocatable :: coldens_out_HI
  real(kind=dp),dimension(:,:,:),allocatable :: coldens_in_HeI
  real(kind=dp),dimension(:,:,:),allocatable :: coldens_out_HeI
  real(kind=dp),dimension(:,:,:),allocatable :: coldens_in_HeII
  real(kind=dp),dimension(:,:,:),allocatable :: coldens_out_HeII

  real(kind=dp),dimension(:,:,:),allocatable :: ring_vol
  ! Buffer for MPI communication
  real(kind=dp),dimension(:,:,:), allocatable :: buffer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Radiation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=dp), dimension(:), allocatable :: delta_freq      ! Frequency width of integration 
  real(kind=dp), dimension(:), allocatable :: freq_max        ! Maximum freqeucny of integration 
  real(kind=dp), dimension(:), allocatable :: freq_min        ! Minimum freqeucny of integration

  ! Power law fit parameter for frequency range 1:3
  real(kind=dp), dimension(:), allocatable :: pl_index_cross_section_HI    ! Power law index of cross section of HI
  real(kind=dp), dimension(:), allocatable :: pl_index_cross_section_HeI   ! Power law index of cross section of HeI
  real(kind=dp), dimension(:), allocatable :: pl_index_cross_section_HeII  ! Power law index of cross section of HeII

  ! Cross section of atoms
  real(kind=dp), dimension(:), allocatable :: sigma_HI       ! Cross section of HI
  real(kind=dp), dimension(:), allocatable :: sigma_HeI      ! Cross section of HeI
  real(kind=dp), dimension(:), allocatable :: sigma_HeII     ! Cross section of HeII

  ! Parameters related to fraction of ionization and heating from different species
  real(kind=dp), dimension(:), allocatable :: f1ion_HI       ! Parameters related to ionization
  real(kind=dp), dimension(:), allocatable :: f1ion_HeI      ! Parameters related to ionization
  real(kind=dp), dimension(:), allocatable :: f1ion_HeII     ! Parameters related to ionization
  real(kind=dp), dimension(:), allocatable :: f2ion_HI       ! Parameters related to ionization
  real(kind=dp), dimension(:), allocatable :: f2ion_HeI      ! Parameters related to ionization
  real(kind=dp), dimension(:), allocatable :: f2ion_HeII     ! Parameters related to ionization
  real(kind=dp), dimension(:), allocatable :: f2heat_HI      ! Parameters related to heating
  real(kind=dp), dimension(:), allocatable :: f2heat_HeI     ! Parameters related to heating
  real(kind=dp), dimension(:), allocatable :: f2heat_HeII    ! Parameters related to heating
  real(kind=dp), dimension(:), allocatable :: f1heat_HI      ! Parameters related to heating
  real(kind=dp), dimension(:), allocatable :: f1heat_HeI     ! Parameters related to heating
  real(kind=dp), dimension(:), allocatable :: f1heat_HeII    ! Parameters related to heating
  
  ! Integrands ( frequency, optical depth )
  real(kind=dp),dimension(:,:), allocatable :: bb_photo_thick_integrand
  real(kind=dp),dimension(:,:), allocatable :: bb_photo_thin_integrand
  real(kind=dp),dimension(:,:), allocatable :: pl_photo_thick_integrand
  real(kind=dp),dimension(:,:), allocatable :: pl_photo_thin_integrand
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thick_integrand_HI
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thick_integrand_HeI
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thick_integrand_HeII
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thin_integrand_HI
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thin_integrand_HeI
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thin_integrand_HeII
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thick_integrand_HI
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thick_integrand_HeI
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thick_integrand_HeII
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thin_integrand_HI
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thin_integrand_HeI
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thin_integrand_HeII

  ! Integration table ( optical depth, sub-bin )
  real(kind=dp),dimension(:,:), target, allocatable :: bb_photo_thick_table  
  real(kind=dp),dimension(:,:), target, allocatable :: bb_photo_thin_table
  real(kind=dp),dimension(:,:), target, allocatable :: pl_photo_thick_table 
  real(kind=dp),dimension(:,:), target, allocatable :: pl_photo_thin_table
  real(kind=dp),dimension(:,:), target, allocatable :: bb_heat_thick_table  
  real(kind=dp),dimension(:,:), target, allocatable :: bb_heat_thin_table
  real(kind=dp),dimension(:,:), target, allocatable :: pl_heat_thick_table  
  real(kind=dp),dimension(:,:), target, allocatable :: pl_heat_thin_table


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=dp),dimension(:),allocatable :: redshift_array 

  integer,dimension(:,:),allocatable :: srcpos
  real(kind=dp),dimension(:),allocatable :: BB_NormFlux !< normalized ionizing flux of sources
  real(kind=dp),dimension(:),allocatable :: PL_NormFlux !< normalized ionizing flux of sources

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! ADP MPI Time-averaged H ionization fraction
  real(kind=dp),dimension(:,:,:),allocatable :: ADP_MPI_intermediate_average_xHI
  real(kind=dp),dimension(:,:,:),allocatable :: ADP_MPI_intermediate_average_xHII
  real(kind=dp),dimension(:,:,:),allocatable :: ADP_MPI_intermediate_average_xHeI
  real(kind=dp),dimension(:,:,:),allocatable :: ADP_MPI_intermediate_average_xHeII
  real(kind=dp),dimension(:,:,:),allocatable :: ADP_MPI_intermediate_average_xHeIII
  ! ADP MPI Intermediate result for H ionization fraction
  real(kind=dp),dimension(:,:,:),allocatable :: ADP_MPI_intermediate_end_xHI
  real(kind=dp),dimension(:,:,:),allocatable :: ADP_MPI_intermediate_end_xHII
  real(kind=dp),dimension(:,:,:),allocatable :: ADP_MPI_intermediate_end_xHeI
  real(kind=dp),dimension(:,:,:),allocatable :: ADP_MPI_intermediate_end_xHeII
  real(kind=dp),dimension(:,:,:),allocatable :: ADP_MPI_intermediate_end_xHeIII
  ! ADP MPI Intermediate result for Temperature
  real(kind=dp),dimension(:,:,:),allocatable :: ADP_MPI_intermediate_end_temperature
  real(kind=dp),dimension(:,:,:),allocatable :: ADP_MPI_intermediate_average_temperature
  ! transforming position from R to RxRxR for ADP evolution
  integer,dimension(:),allocatable :: ADP_m_to_i,ADP_m_to_j,ADP_m_to_k

  ! LTE MPI Time-averaged H ionization fraction
  real(kind=dp),dimension(:,:,:),allocatable :: LTE_MPI_intermediate_average_xHI
  real(kind=dp),dimension(:,:,:),allocatable :: LTE_MPI_intermediate_average_xHII
  real(kind=dp),dimension(:,:,:),allocatable :: LTE_MPI_intermediate_average_xHeI
  real(kind=dp),dimension(:,:,:),allocatable :: LTE_MPI_intermediate_average_xHeII
  real(kind=dp),dimension(:,:,:),allocatable :: LTE_MPI_intermediate_average_xHeIII
  ! LTE MPI Intermediate result for H ionization fraction
  real(kind=dp),dimension(:,:,:),allocatable :: LTE_MPI_intermediate_end_xHI
  real(kind=dp),dimension(:,:,:),allocatable :: LTE_MPI_intermediate_end_xHII
  real(kind=dp),dimension(:,:,:),allocatable :: LTE_MPI_intermediate_end_xHeI
  real(kind=dp),dimension(:,:,:),allocatable :: LTE_MPI_intermediate_end_xHeII
  real(kind=dp),dimension(:,:,:),allocatable :: LTE_MPI_intermediate_end_xHeIII
  ! LTE MPI Intermediate result for Temperature
  real(kind=dp),dimension(:,:,:),allocatable :: LTE_MPI_intermediate_end_temperature
  real(kind=dp),dimension(:,:,:),allocatable :: LTE_MPI_intermediate_average_temperature
  ! transforming position from R to RxRxR for LTE evolution
  integer,dimension(:),allocatable :: LTE_m_to_i,LTE_m_to_j,LTE_m_to_k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Adaptive Timestep !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! targeting ionized fractions to be achieved 
  real, dimension(1:time_partition) :: xfinal
  ! adaptive timestep derived from all sources
  real(kind=dp), dimension(:), allocatable :: dt_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! LTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer,dimension(1:mesh,1:mesh,1:mesh) :: global_Ifront_array
  real(kind=dp),dimension(1:mesh,1:mesh,1:mesh) :: global_LTE_array = -1.0
  integer,dimension(1:mesh,1:mesh,1:mesh) :: equivalent_class_index_array
  integer,dimension(1:mesh,1:mesh,1:mesh) :: global_evolved_LTE_array = 0
  integer,dimension(:,:,:),allocatable :: source_Ifront_array
  integer,dimension(:,:,:),allocatable :: ionized_array
  integer,dimension(:,:,:),allocatable :: ifront_plus_ionized_array
  integer,dimension(:),allocatable :: source_class_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! AT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer,dimension(:,:,:),allocatable :: non_screened_Ifront_array
  real(kind=dp),dimension(:,:,:),allocatable :: global_photoionization_HI_array
  real(kind=dp),dimension(:,:,:),allocatable :: global_photoionization_HeI_array
  real(kind=dp),dimension(:,:,:),allocatable :: global_photoionization_HeII_array
  real(kind=dp),dimension(:,:,:),allocatable :: global_timestep_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!AT geometry arrays !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer, dimension(:,:,:), allocatable :: pos_x_pos_y_pos_z_octant
  integer, dimension(:,:,:), allocatable :: pos_x_pos_y_neg_z_octant
  integer, dimension(:,:,:), allocatable :: pos_x_neg_y_pos_z_octant
  integer, dimension(:,:,:), allocatable :: pos_x_neg_y_neg_z_octant
  integer, dimension(:,:,:), allocatable :: neg_x_pos_y_pos_z_octant
  integer, dimension(:,:,:), allocatable :: neg_x_pos_y_neg_z_octant
  integer, dimension(:,:,:), allocatable :: neg_x_neg_y_pos_z_octant
  integer, dimension(:,:,:), allocatable :: neg_x_neg_y_neg_z_octant

  integer, dimension(:), allocatable, target :: pos_x_pos_y_pos_z_dir_x_triangular_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_pos_z_dir_x_triangular_x_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_pos_z_dir_x_triangular_y_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_pos_z_dir_x_triangular_z_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_pos_z_dir_y_triangular_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_pos_z_dir_y_triangular_x_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_pos_z_dir_y_triangular_y_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_pos_z_dir_y_triangular_z_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_pos_z_dir_z_triangular_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_pos_z_dir_z_triangular_x_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_pos_z_dir_z_triangular_y_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_pos_z_dir_z_triangular_z_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_neg_z_dir_x_triangular_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_neg_z_dir_x_triangular_x_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_neg_z_dir_x_triangular_y_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_neg_z_dir_x_triangular_z_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_neg_z_dir_y_triangular_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_neg_z_dir_y_triangular_x_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_neg_z_dir_y_triangular_y_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_neg_z_dir_y_triangular_z_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_neg_z_dir_z_triangular_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_neg_z_dir_z_triangular_x_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_neg_z_dir_z_triangular_y_array
  integer, dimension(:), allocatable, target :: pos_x_pos_y_neg_z_dir_z_triangular_z_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_pos_z_dir_x_triangular_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_pos_z_dir_x_triangular_x_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_pos_z_dir_x_triangular_y_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_pos_z_dir_x_triangular_z_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_pos_z_dir_y_triangular_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_pos_z_dir_y_triangular_x_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_pos_z_dir_y_triangular_y_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_pos_z_dir_y_triangular_z_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_pos_z_dir_z_triangular_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_pos_z_dir_z_triangular_x_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_pos_z_dir_z_triangular_y_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_pos_z_dir_z_triangular_z_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_neg_z_dir_x_triangular_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_neg_z_dir_x_triangular_x_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_neg_z_dir_x_triangular_y_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_neg_z_dir_x_triangular_z_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_neg_z_dir_y_triangular_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_neg_z_dir_y_triangular_x_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_neg_z_dir_y_triangular_y_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_neg_z_dir_y_triangular_z_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_neg_z_dir_z_triangular_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_neg_z_dir_z_triangular_x_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_neg_z_dir_z_triangular_y_array
  integer, dimension(:), allocatable, target :: pos_x_neg_y_neg_z_dir_z_triangular_z_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_pos_z_dir_x_triangular_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_pos_z_dir_x_triangular_x_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_pos_z_dir_x_triangular_y_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_pos_z_dir_x_triangular_z_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_pos_z_dir_y_triangular_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_pos_z_dir_y_triangular_x_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_pos_z_dir_y_triangular_y_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_pos_z_dir_y_triangular_z_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_pos_z_dir_z_triangular_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_pos_z_dir_z_triangular_x_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_pos_z_dir_z_triangular_y_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_pos_z_dir_z_triangular_z_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_neg_z_dir_x_triangular_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_neg_z_dir_x_triangular_x_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_neg_z_dir_x_triangular_y_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_neg_z_dir_x_triangular_z_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_neg_z_dir_y_triangular_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_neg_z_dir_y_triangular_x_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_neg_z_dir_y_triangular_y_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_neg_z_dir_y_triangular_z_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_neg_z_dir_z_triangular_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_neg_z_dir_z_triangular_x_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_neg_z_dir_z_triangular_y_array
  integer, dimension(:), allocatable, target :: neg_x_pos_y_neg_z_dir_z_triangular_z_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_pos_z_dir_x_triangular_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_pos_z_dir_x_triangular_x_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_pos_z_dir_x_triangular_y_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_pos_z_dir_x_triangular_z_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_pos_z_dir_y_triangular_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_pos_z_dir_y_triangular_x_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_pos_z_dir_y_triangular_y_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_pos_z_dir_y_triangular_z_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_pos_z_dir_z_triangular_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_pos_z_dir_z_triangular_x_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_pos_z_dir_z_triangular_y_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_pos_z_dir_z_triangular_z_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_neg_z_dir_x_triangular_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_neg_z_dir_x_triangular_x_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_neg_z_dir_x_triangular_y_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_neg_z_dir_x_triangular_z_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_neg_z_dir_y_triangular_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_neg_z_dir_y_triangular_x_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_neg_z_dir_y_triangular_y_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_neg_z_dir_y_triangular_z_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_neg_z_dir_z_triangular_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_neg_z_dir_z_triangular_x_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_neg_z_dir_z_triangular_y_array
  integer, dimension(:), allocatable, target :: neg_x_neg_y_neg_z_dir_z_triangular_z_array
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! photonstatistics !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(kind=dp), dimension(:,:), allocatable :: LTE_photon_loss_array
  real(kind=dp), dimension(:,:), allocatable :: ADP_photon_loss_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module array
