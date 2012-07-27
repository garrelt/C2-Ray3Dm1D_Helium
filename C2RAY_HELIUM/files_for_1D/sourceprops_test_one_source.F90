module sourceprops

  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use cgsconstants, only: m_p
  use astroconstants, only: M_SOLAR, YEAR
  use cosmology_parameters, only: Omega_B, Omega0
  !use nbody, only:  dir_src
  !use material, only: xh
  use grid, only: x,y,z
  use c2ray_parameters, only: phot_per_atom, lifetime, &
       S_star_nominal, StillNeutral, Number_Sourcetypes
  real(kind=dp),dimension(:,:),allocatable :: srcMass !< masses of sources 
  integer :: NumSrc 
  integer,dimension(:,:),allocatable :: srcpos
  real(kind=dp),dimension(:,:),allocatable :: rsrcpos
  real(kind=dp),dimension(:),allocatable :: NormFlux,NormFluxPL
  integer,dimension(:),allocatable :: srcSeries

contains
  
  ! =======================================================================

  ! =======================================================================

  !> Initialization routine: dummy
  !! Author: Garrelt Mellema
  
  subroutine source_properties_ini (sourcetype)
    use file_admin, only: stdinput
    character :: sourcetype
    allocate(NormFlux(1))
    allocate(NormFluxPL(1))
    allocate(srcMass(1,1:3))
    NormFlux(1)=1.0_dp
    NormFluxPL(1)=0.0_dp
    if (sourcetype.eq.'P')    NormFluxPL(1)=1.0_dp

    srcMass(1,:)=(/1.0_dp, 1.0_dp, 1.0_dp/)  !not needed?
  end subroutine source_properties_ini

end module sourceprops
