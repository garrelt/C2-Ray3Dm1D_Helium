!>
!! \brief This module contains data and routines for handling the source
!! luminosities.
!! 
!! \b Author: Garrelt Mellema, Ilian Iliev
!!
!! \b Date: 30-Jan-2008
!!
!! \b Version: Test simulation with reading in file with source locations
!! \b and ionizing photon rates

module sourceprops
  
  use precision, only: dp
  use my_mpi
  use file_admin, only: logf,sourcefile
  use cgsconstants, only: m_p
  use astroconstants, only: M_SOLAR, YEAR
  use cosmology_parameters, only: Omega_B, Omega0
  use nbody, only:  dir_src
  !use material, only: xh
  use grid, only: x,y,z
#ifdef PL
  use c2ray_parameters, only: xray_phot_per_atom, pl_S_star_nominal,&
       EddLeff_nominal,mass_nominal,EddLum
#endif
#ifdef QUASARS
  use c2ray_parameters, only: qpl_S_star_nominal,qEddLeff_nominal,&
      qmass_nominal,qEddLum, qpl_index_nominal, qpl_MinFreq_nominal,&
      qpl_MaxFreq_nominal
#endif
  use c2ray_parameters, only: phot_per_atom, lifetime, &
       S_star_nominal, StillNeutral, Number_Sourcetypes
  
  implicit none

  !> base name of source list files
  character(len=100),parameter,private :: &
       sourcelistfile_base="test_sources.dat"

  integer :: NumSrc=0 !< Number of sources 
  integer,dimension(:,:),allocatable :: srcpos
  real(kind=dp),dimension(:),allocatable :: NormFlux !< normalized ionizing flux of sources
  real(kind=dp),dimension(:),allocatable :: NormFluxPL !< normalized ionizing flux of sources
  integer,dimension(:),allocatable :: srcSeries  !< a randomized list of sources

  character(len=512),private :: sourcelistfile,sourcelistfilesuppress

contains
  
  ! =======================================================================

  !> Set the source properties for this redshift
  subroutine source_properties(zred_now,nz,lifetime2,restart)

    ! Input routine: establish the source properties

    ! For random permutation of sources
    use  m_ctrper

    real(kind=dp),intent(in) :: zred_now ! current redshift
    real(kind=dp),intent(in) :: lifetime2 ! time step
    integer,intent(in) :: nz
    integer,intent(in) :: restart

    integer :: ns,ns0

#ifdef MPI
    integer :: mympierror
#endif

#ifdef MPILOG     
    write(logf,*) "Check sourceprops: ",zred_now,nz,lifetime2,restart
#endif 
    
    ! Deallocate source arrays
    if (allocated(srcpos)) deallocate(srcpos)
    if (allocated(NormFlux)) deallocate(NormFlux)
    if (allocated(srcSeries)) deallocate(srcSeries)
#ifdef PL
    if (allocated(NormFluxPL)) deallocate(NormFluxPL)
#endif
    
    ! Rank 0 reads in sources
    if (rank == 0) then
       ! Construct the file names
       sourcelistfile=trim(adjustl(dir_src))//"test_sources.dat"
       open(unit=sourcefile,file=sourcelistfile,status="old")

       ! Establish number of sources
       read(sourcefile,*) NumSrc

       write(*,*)  "Number of sources: ",NumSrc

    endif ! end of rank 0 test
    
#ifdef MPI
    ! Distribute source number to all other nodes
    call MPI_BCAST(NumSrc,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

#ifdef MPILOG
    if (rank /=0) write(logf,*) "Number of sources, with suppression: ",NumSrc
#endif

    ! Allocate arrays for this NumSrc
    if (NumSrc > 0) then
       write(*,*)  "Reading number of sources: ",NumSrc
       allocate(srcpos(3,NumSrc))
       allocate(NormFlux(0:NumSrc)) ! 0 will hold lost photons
       allocate(SrcSeries(NumSrc))
#ifdef PL
       allocate(NormFluxPL(NumSrc))
#endif

       ! Fill in the source arrays
       if (rank == 0) then
          do ns=1,NumSrc
#ifdef PL
             read(sourcefile,*) srcpos(1,ns),srcpos(2,ns),srcpos(3,ns), &
                  NormFlux(ns),NormFluxPL(ns)
             NormFluxPL(ns)=NormFluxPL(ns)/pl_S_star_nominal
#else
             read(sourcefile,*) srcpos(1,ns),srcpos(2,ns),srcpos(3,ns), &
                  NormFlux(ns)
#endif
             NormFlux(ns)=NormFlux(ns)/S_star_nominal
          enddo
          
          ! STANDARD TEST
          ! Source positions in file start at 1!
          !srcpos(1:3,1)=(/ 50, 50, 50 /)
          !srcpos(1:3,2)=(/ 51, 50, 50 /)
          !srcpos(1:3,3)=(/ 52, 50, 50 /)
          !srcpos(1:3,4)=(/ 53, 50, 50 /)
          !NormFlux(1:4)=1e55_dp/S_star_nominal
          !NormFluxPL(1:4)=0.0_dp
          
          !srcpos(1:3,5)=(/ 20, 10, 10 /)
          !NormFlux(5)=1e57_dp/S_star_nominal
          !NormFluxPL(5)=0.0_dp
          
          !srcpos(1:3,6)=(/ 70, 70, 50 /)
          !srcpos(1:3,7)=(/ 72, 70, 50 /)
          !srcpos(1:3,8)=(/ 70, 72, 50 /)
          !srcpos(1:3,9)=(/ 72, 72, 50 /)
          !NormFlux(6:8)=1e55_dp/S_star_nominal
          !NormFlux(9)=1e56_dp/S_star_nominal
          !NormFluxPL(6:9)=0.0
          
          !srcpos(1:3,10)=(/ 20, 10, 90 /)
          !NormFlux(10)=1e54_dp/S_star_nominal
          ! Note, no normalization is used for PL source
          !NormFluxPL(10)=3.0_dp
          
          write(logf,*) 'Total photon rate (BB)= ', &
               sum(NormFlux)*S_star_nominal,' s^-1'
#ifdef PL
          write(logf,*) 'Total photon rate (PL)= ', &
               sum(NormFluxPL)*pl_S_star_nominal,' s^-1'
#endif

          ! Create array of source numbers for generating random order
          do ns=1,NumSrc
             SrcSeries(ns)=ns
          enddo
          ! Make a random order
          call ctrper(SrcSeries(1:NumSrc),1.0)
       endif ! of rank 0 test

#ifdef MPI
       ! Distribute the source parameters to the other nodes
       call MPI_BCAST(srcpos,3*NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(NormFlux,NumSrc+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#ifdef PL
       call MPI_BCAST(NormFluxPL,NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif
       ! Distribute the source series to the other nodes
       call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
       
    endif

    ! close the source file if you are rank 0
    if (rank == 0) close(sourcefile)

  end subroutine source_properties
  
  ! =======================================================================
  
  !> Initialization routine: dummy
  !! Author: Garrelt Mellema
  !! I'm using it now for the 1D-code
  subroutine source_properties_ini ()
  end subroutine source_properties_ini

end module sourceprops
