!>
!! \brief This module contains data and routines for handling the data from
!! the Nbody simulations which provide the basis for C2Ray simulations
!! of Reionization.
!! 
!! \b Author: Garrelt Mellema, Ilian Iliev
!!
!! \b Date: 22-May-2008 (previous version was not dated)
!!
!! \b Version: test simulations
!!
!! The test module has a number of parameters hard coded. These are
!!
!! Starting redshift: 9
!!
!! Number of redshift slices: 8 (separated by 0.05 million years).
!!
!! Box size: 0.5/h Mpc

module nbody

  ! This file contains routine having to do with the input N-body simulation.
  ! The routine in here is
  ! - nbody_ini (called by main program)
  ! It reads the list of redshifts for which source lists and
  ! density fields are available.
  ! It calculates the mass (in particles) of the simulation box
  ! It sets an identifier (id_str) for the resolution (used when
  ! reading in source list files and density fields.

  ! Authors: Ilian Iliev, Garrelt Mellema

  ! Date: 22-May-2008 (previous version was not dated)

  use precision, only: dp
  use sizes, only: mesh
  use file_admin, only: stdinput, logf, file_input
  use astroconstants, only: Mpc, M_SOLAR, YEAR
  use my_mpi
  use cosmology_parameters, only: rho_crit_0, Omega0, h, H0

  implicit none

  character(len=10),parameter :: nbody_type="test4" !< ID of Nbody type

  real(kind=dp),parameter :: boxsize=0.5  !< Box size in Mpc/h comoving

  ! redshift sequence information
  real,parameter :: zred_start=8.8492 !< starting redshift  ! was 9.0
  integer, parameter,public :: NumZred=9    !< number of redshifts
  real,public :: timestep=0.05e6*YEAR !< time step

  real(kind=dp),dimension(NumZred),public :: zred_array !< array of redshifts
  integer,dimension(:),allocatable,public :: snap !< array of snapshot numbers (for compatibility)
  character(len=9),public :: id_str='test4 res'       !< resolution dependent string
  
  character(len=480),public :: dir_dens !< Path to directory with density files
  character(len=480),public :: dir_clump !< Path to directory with density files
  character(len=480),public :: dir_LLS !< Path to directory with LLS files  
  character(len=480),public :: dir_src !< Path to directory with source files

  !> Path to directory containing directory with density files:
  character(len=*),parameter,private :: dir_dens_path = "../TEST4/" 
  !> Name of directory with density files
  character(len=180),parameter,private :: dir_dens_name= ""

  !> Path to directory containing directory with clumping files:
  character(len=*),parameter,private :: dir_clump_path = "./" 
  !> Name of directory with files used for clumping
  character(len=*),parameter,private :: dir_clump_name= "./"

  !> Path to directory containing directory with source files:
  character(len=*),parameter,private :: dir_src_path = "../TEST4/" 
  !> Name of directory with source files
  character(len=*),parameter,private :: dir_src_name= ""

  !> Format of density file (unformatted or binary)
#ifdef IFORT
  ! ifort standard for "binary"
  character(len=*),parameter :: densityformat="unformatted"
  character(len=*),parameter :: densityaccess="sequential"
#else
  ! Fortran2003 standard for "binary"
!  character(len=*),parameter :: densityformat="unformatted"
!  character(len=*),parameter :: densityaccess="sequential"
#endif
  !> Format of clumping file (unformatted or binary)
#ifdef IFORT
  character(len=*),parameter :: clumpingformat="binary"
  character(len=*),parameter :: clumpingaccess="sequential"
#else
  character(len=15),parameter :: clumpingformat="unformatted"
  character(len=*),parameter :: clumpingaccess="sequential"
#endif
  !> Format of LLS file (unformatted or binary)
#ifdef IFORT
  character(len=*),parameter :: LLSformat="binary"
  character(len=*),parameter :: LLSaccess="sequential"
#else
  character(len=15),parameter :: LLSformat="unformatted"
  character(len=*),parameter :: LLSaccess="stream"
#endif
  !> density file with header?
  logical,parameter :: densityheader=.true.
  !> clumping file with header?
  logical,parameter :: clumpingheader=.true.
  !> LLS file with header?
  logical,parameter :: LLSheader=.true.    
  !> unit of density in density file
  !! can be "grid", "particle", "M0Mpc3"
  character(len=*),parameter :: density_unit="none"

#ifdef MPI
  integer,private :: mympierror !< MPI error flag variable
#endif

contains

  ! ===========================================================================

  subroutine nbody_ini ()
    
    real(kind=dp) :: t0!,timestep
    integer :: nz ! loop counter
    character(len=256) :: value
    integer :: len, status

    ! Set the various input directories
    dir_dens=trim(adjustl(dir_dens_path))//trim(adjustl(dir_dens_name))
    dir_clump=trim(adjustl(dir_clump_path))//trim(adjustl(dir_clump_name))
    dir_src=trim(adjustl(dir_src_path))//trim(adjustl(dir_src_name))

    write(*,*) dir_dens
    write(*,*) dir_src
    pause


    ! Construct redshift sequence

    ! Starting redshift
    zred_array(1)=zred_start
    
           ! Time step
       !timestep=1e7*YEAR

    ! Cosmological time corresponding to (initial) redshift
    ! NOTE: Good only for high-z!!!
    t0 = 2.*(1.+zred_array(1))**(-1.5)/(3.*H0*sqrt(Omega0))
    do nz=2,NumZred
       zred_array(nz)=-1+(1.+zred_array(1))* &
            (t0/(t0+real(nz-1)*timestep))**(2./3.)
    enddo

  end subroutine nbody_ini

end module nbody
