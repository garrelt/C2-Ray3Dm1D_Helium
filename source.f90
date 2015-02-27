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
  use file_admin, only: logf
  use cgsconstants, only: m_p
  use astroconstants, only: M_SOLAR, YEAR
  use cosmology_parameters, only: Omega_B, Omega0
  use c2ray_parameters, only: phot_per_atom, lifetime, &
       bb_S_star_nominal, pl_S_star_nominal,StillNeutral, Number_Sourcetypes
  use parameter, only: number_of_source
  use array, only: srcpos,BB_NormFlux,PL_NormFlux

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine source_initialization()

    implicit none

    character(len=512) :: sourcelistfile
    integer :: ns,ns0
    integer :: mympierror

    ! Deallocate source arrays
    if (allocated(srcpos)) deallocate(srcpos)
    if (allocated(BB_NormFlux)) deallocate(BB_NormFlux)
    if (allocated(PL_NormFlux)) deallocate(PL_NormFlux)
    
    ! Rank 0 reads in sources
    if (rank .eq. 0) then

       ! Construct the file names
       !sourcelistfile="test_sources.dat"
       sourcelistfile="sources_test4.dat" ! Test 4

       write(logf,*) "Reading source file from ",trim(adjustl(sourcelistfile))

       open(unit=50,file=sourcelistfile,status="old")

       ! Establish number of sources
       read(50,*) number_of_source
       write(logf,*) "Number of source = ",number_of_source
       flush(logf) 

    endif 
    
#ifdef MPI
    ! Distribute source number to all other nodes
    call MPI_BCAST(number_of_source,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

    ! Allocate arrays for this number_of_source
    if (number_of_source > 0) then
      allocate(srcpos(3,1:number_of_source))
      allocate(BB_NormFlux(0:number_of_source)) 
      allocate(PL_NormFlux(1:number_of_source))

      ! Fill in the source arrays
      if (rank .eq. 0) then
        BB_NormFlux(0) = 0.0

        do ns = 1,number_of_source
          read(50,*) srcpos(1,ns),srcpos(2,ns),srcpos(3,ns),BB_NormFlux(ns),PL_NormFlux(ns)
          BB_NormFlux(ns)=BB_NormFlux(ns)/bb_S_star_nominal
          PL_NormFlux(ns)=PL_NormFlux(ns)/pl_S_star_nominal
        enddo

        close(50)

      endif ! of rank 0 test

#ifdef MPI
       ! Distribute the source parameters to the other nodes
       call MPI_BCAST(srcpos,3*number_of_source,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(BB_NormFlux,number_of_source+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(PL_NormFlux,number_of_source,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif
       
    endif

  end subroutine source_initialization
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end module sourceprops
