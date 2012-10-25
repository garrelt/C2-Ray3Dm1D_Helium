!>
!! \brief This module contains data and routines for handling the physical grid
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2012-10-13 (but older)
!!
!! \b Version: 1D, radial (spherically symmetric) grid.

module grid

  ! Handles grid properties
  ! Data: grid quantities: cell size, spatial coordinates, cell volume
  ! Routines: grid initialization: grid_ini

  ! Notes:
  ! Cell size is called dr, spatial coordinate is called x (but is really
  ! radial coordinate r). The latter is for consistency with the 3D version
  ! of the C2-Ray code.
  ! Also for consistency, the length of the arrays is mesh(1).

  use precision, only: dp
  use sizes, only: Ndim, mesh
  use astroconstants, only: Mpc
  use my_mpi
  use file_admin, only: stdinput, file_input

  implicit none

  ! Contains grid data
  ! dr - (radial) cell size
  ! x - spatial coordinates (radial)
  ! vol - volume of one cell
  real(kind=dp) :: dr(1) !< cell size
  real(kind=dp),dimension(mesh(1)) :: x !< spatial cooridnate   (was r)
  real(kind=dp),dimension(mesh(1)) :: y !< spatial cooridnate (compatibility)
  real(kind=dp),dimension(mesh(1)) :: z !< spatial cooridnate (compatibility)
  real(kind=dp),dimension(mesh(1)) :: vol !< volume of grid cell
  
contains

  ! =======================================================================

  !> Initializes grid properties\n
  !! \b Author: Garrelt Mellema\n
  !! \b Date: 20-Aug-2006 (f77 version: 15-Apr-2004)\n
  !! \b Version: One-dimensional spherical grid
    
  subroutine grid_ini()
    
    ! Initializes grid properties
    
    ! Author: Garrelt Mellema
    
    ! Date: 20-Aug-2006 (f77 version: 15-Apr-2004)
    
    ! Version: One-dimensional spherical grid
    
    use mathconstants, only: pi
    use string_manipulation, only: convert_case
    use astroconstants, only: pc,kpc,Mpc,AU

    integer :: i
    real(kind=dp) :: r_in,r_out
    character(len=10) :: str_length_unit
    real(kind=dp) :: conversion_factor
    
#ifdef MPI
    integer :: ierror
#endif

    ! Ask for grid size
    if (rank == 0) then
       if (.not.file_input) then
          write(*,*) 'Note: for cosmological applications, specify'
          write(*,*) 'comoving values below.'
          write(*,'(A,$)') 'Enter inner and outer radius of grid (specify units): '
       endif
       read(stdinput,*) r_in,r_out,str_length_unit
       
       ! Convert to cms
       call convert_case(str_length_unit,0) ! conversion to lower case
       select case (trim(adjustl(str_length_unit)))
       case ('cm','centimeter','cms','centimeters')
          conversion_factor=1.0
       case ('m','meter','ms','meters')
          conversion_factor=100.0
       case ('km','kilometer','kms','kilometers','clicks')
          conversion_factor=1.0e5
       case ('pc','parsec','parsecs')
          conversion_factor=pc
       case ('kpc','kiloparsec','kiloparsecs')
          conversion_factor=kpc
       case ('mpc','megaparsec','megaparsecs')
          conversion_factor=Mpc
       case default
          write(*,*) 'Length unit not recognized, assuming cm'
          conversion_factor=1.0
       end select
       r_in=r_in*conversion_factor
       r_out=r_out*conversion_factor
    endif
    
#ifdef MPI
    ! Distribute the total grid size over all processors
    call MPI_BCAST(r_in,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(r_out,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
#endif
    
    dr(1)=(r_out-r_in)/real(mesh(1))
      
    ! Radial coordinate of a cell
    do i=1,mesh(1)
       x(i)=(real(i)-0.5)*dr(1)+r_in
    enddo

    ! Volume of a cell.
    do i=1,mesh(1)
       vol(i)=4.0*pi/3.0*((x(i)+0.5*dr(1))**3-(x(i)-0.5*dr(1))**3)
       !vol(i)=4.0*pi*r(i)*r(i)*dr
    enddo

  end subroutine grid_ini
  
end module grid
