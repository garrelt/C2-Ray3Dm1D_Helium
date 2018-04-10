!>
!! \brief This module contains the data for the evolve modules which
!! calculate the ionization and temperature evolution of the entire grid (3D).
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date:
!!
!! \b Version: 3D, MPI & OpenMP

module evolve_data

  ! This version has been adapted for efficiency in order to be able
  ! to calculate large meshes.
    
  use precision, only: dp
  use my_mpi
  use sizes, only: mesh, Ndim
  use radiation_sizes, only: NumFreqBnd
  
  implicit none

  save

  !> Periodic boundary conditions, has to be true for this version
  logical,parameter :: periodic_bc = .true.

  !> Minimum number of MPI processes for using the master-slave setup 
  integer, parameter ::  min_numproc_master_slave=10

  ! GM/121127: This variable should always be set. If not running OpenMP
  ! it should be equal to 1. We initialize it to 1 here.
  integer :: tn=1 !< thread number

  ! Parameters for multiphase
  real(kind=dp),parameter :: epsilon_dx=1e-5 ! limit for considering special
  real(kind=dp),parameter :: limit_partial_cells=0.99 ! limit for turning off special
    
  ! Grid variables

  !> H Photo-ionization rate on the entire grid
  real(kind=dp),dimension(:,:,:),allocatable :: bb_phih_grid
  !> He Photo-ionization rate on the entire grid
  real(kind=dp),dimension(:,:,:,:),allocatable :: bb_phihe_grid
  !> Heating  rate on the entire grid
  real(kind=dp),dimension(:,:,:),allocatable :: bb_phiheat_grid
#ifdef PL
  !> H Photo-ionization rate on the entire grid
  real(kind=dp),dimension(:,:,:),allocatable :: pl_phih_grid
  !> He Photo-ionization rate on the entire grid
  real(kind=dp),dimension(:,:,:,:),allocatable :: pl_phihe_grid
  !> Heating  rate on the entire grid
  real(kind=dp),dimension(:,:,:),allocatable :: pl_phiheat_grid
#endif

#ifdef QUASARS
  !> H Photo-ionization rate on the entire grid
  real(kind=dp),dimension(:,:,:),allocatable :: qpl_phih_grid
  !> He Photo-ionization rate on the entire grid
  real(kind=dp),dimension(:,:,:,:),allocatable :: qpl_phihe_grid
  !> Heating  rate on the entire grid
  real(kind=dp),dimension(:,:,:),allocatable :: qpl_phiheat_grid
#endif
  
  !> Time-averaged H ionization fraction
  real(kind=dp),dimension(:,:,:,:),allocatable :: xh_av
  !> Time-averaged fraction of cell which is hot and ionized
  real(kind=dp),dimension(:,:,:),allocatable,target :: xh_hot_av
  !> Time-averaged He ionization fraction
  real(kind=dp),dimension(:,:,:,:),allocatable :: xhe_av
  !> Intermediate result for H ionization fraction
  real(kind=dp),dimension(:,:,:,:),allocatable :: xh_intermed
  !> Intermediate result for fraction of cell which is hot and ionized
  real(kind=dp),dimension(:,:,:),allocatable :: xh_hot_intermed
  !> Intermediate result for He ionization fraction
  real(kind=dp),dimension(:,:,:,:),allocatable :: xhe_intermed
  !> H0 Column density (outgoing)
  real(kind=dp),dimension(:,:,:),allocatable :: coldensh_out
  !> He Column density (outgoing)
  real(kind=dp),dimension(:,:,:,:),allocatable :: coldenshe_out
  !> Buffer for MPI communication
  real(kind=dp),dimension(:,:,:),allocatable :: buffer
  !> Photon loss from the grid
  real(kind=dp) :: photon_loss_all(1:NumFreqBnd)
  !> Photon loss from one source
  real(kind=dp),dimension(:),allocatable :: photon_loss_src_thread
  real(kind=dp) :: photon_loss_src

  integer,dimension(Ndim) :: last_l !< mesh position of left end point for RT
  integer,dimension(Ndim) :: last_r !< mesh position of right end point for RT

contains

  ! =======================================================================

  !> Allocate the arrays needed for evolve
  subroutine evolve_ini ()
    
    allocate(bb_phih_grid(mesh(1),mesh(2),mesh(3)))
    bb_phih_grid=0.0 ! Needs value for initial output
    allocate(bb_phihe_grid(mesh(1),mesh(2),mesh(3),0:1))
    bb_phihe_grid=0.0 ! Needs value for initial output
    allocate(bb_phiheat_grid(mesh(1),mesh(2),mesh(3)))
    bb_phiheat_grid=0.0 ! Needs value for initial output

#ifdef PL
    allocate(pl_phih_grid(mesh(1),mesh(2),mesh(3)))
    pl_phih_grid=0.0 ! Needs value for initial output
    allocate(pl_phihe_grid(mesh(1),mesh(2),mesh(3),0:1))
    pl_phihe_grid=0.0 ! Needs value for initial output
    allocate(pl_phiheat_grid(mesh(1),mesh(2),mesh(3)))
    pl_phiheat_grid=0.0 ! Needs value for initial output
#endif

#ifdef QUASARS
    allocate(qpl_phih_grid(mesh(1),mesh(2),mesh(3)))
    qpl_phih_grid=0.0 ! Needs value for initial output
    allocate(qpl_phihe_grid(mesh(1),mesh(2),mesh(3),0:1))
    qpl_phihe_grid=0.0 ! Needs value for initial output
    allocate(qpl_phiheat_grid(mesh(1),mesh(2),mesh(3)))
    qpl_phiheat_grid=0.0 ! Needs value for initial output
#endif
    
    allocate(xh_av(mesh(1),mesh(2),mesh(3),0:1))
    allocate(xh_hot_av(mesh(1),mesh(2),mesh(3)))
    allocate(xhe_av(mesh(1),mesh(2),mesh(3),0:2))

    allocate(xh_intermed(mesh(1),mesh(2),mesh(3),0:1))
    allocate(xh_hot_intermed(mesh(1),mesh(2),mesh(3)))
    allocate(xhe_intermed(mesh(1),mesh(2),mesh(3),0:2))

    allocate(coldensh_out(mesh(1),mesh(2),mesh(3)))
    allocate(coldenshe_out(mesh(1),mesh(2),mesh(3),0:1))

#ifdef MPI
    allocate(buffer(mesh(1),mesh(2),mesh(3)))
#endif
    allocate(photon_loss_src_thread(nthreads))
  !   allocate(photon_loss_src_thread(1:NumFreqBnd,1))   

  end subroutine evolve_ini

end module evolve_data
