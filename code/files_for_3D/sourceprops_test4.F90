!>
!! \brief This module contains data and routines for handling the sources
!! of reionization.
!!
!! 
!! \b Author: Garrelt Mellema, Ilian Iliev
!!
!! \b Date: 30-Jan-2008
!!
!! \b Version: CUBEP3M Nbody simulation, with source suppression, and different source models

module sourceprops

  use precision, only: dp
  use my_mpi
  use file_admin, only: stdinput, logf, file_input, sourcefile
  use cgsconstants, only: m_p
  use astroconstants, only: M_SOLAR, YEAR
  use cosmology_parameters, only: Omega_B, Omega0
  use nbody, only: id_str, dir_src, NumZred
  use material, only: xh
  use grid, only: x,y,z
  use c2ray_parameters, only: phot_per_atom, lifetime, &
       S_star_nominal, StillNeutral, Number_Sourcetypes

  implicit none

  !> base name of source list files
  character(len=100),private :: sourcelistfile_base="_sources.dat"
  !character(len=100),private :: sourcelistfile_base="_wsubgrid_sources.dat"
  character(len=100),private :: sourcelistfilesuppress_base="_sources_used_wfgamma.dat"

  !> maximum increase in uv to use up cumulated photons
  real(kind=dp),parameter,private :: cumfrac_max=0.15 

  integer :: NumSrc=0 !< Number of sources
  integer :: Prev_NumSrc !< Previous number of sources
  integer,dimension(:,:),allocatable :: srcpos !< mesh position of sources
  real(kind=dp),dimension(:,:),allocatable :: rsrcpos !< grid position of sources
  real(kind=dp),dimension(:,:),allocatable :: srcMass !< masses of sources 
  real(kind=dp),dimension(:),allocatable :: NormFlux !< normalized ionizing flux of sources
   real(kind=dp),dimension(:),allocatable :: NormFluxPL !< normalized ionizing flux of sources 
  integer,dimension(:),allocatable :: srcSeries  !< a randomized list of sources
  real(kind=dp),dimension(:),allocatable :: uv_array  !< list of UV flux evolution (for some sources models)
  !> The cumulative number of uv photons. We save this number so we can add it
  !! to the uv luminosity the first time sources appear.
  real(kind=dp) :: cumulative_uv=0.0

  character(len=30) :: UV_Model !< type of UV model
  integer :: NumZred_uv !< Number of redshift points in UV model
  integer,private :: NumSrc0=0 !< intermediate source count
  integer,dimension(3),private :: srcpos0
  real(kind=dp),private :: srcMass00,srcMass01,total_SrcMass
  character(len=6) :: z_str !< string value of redshift
  integer,private :: NumMassiveSrc !< counter: number of massive sources
  integer,private :: NumSupprbleSrc !< counter: number of suppressible sources
  integer,private :: NumSupprsdSrc !< counter: number of suppressed sources
  real(kind=dp),private :: cumfrac

  character(len=512),private :: sourcelistfile,sourcelistfilesuppress

contains
  
  ! =======================================================================

  !> Set the source properties for this redshift
  !! Authors: Garrelt Mellema, Ilian Iliev
  !! Update: 30-Jan-2008 (20-Sep-2006 (3-jan-2005, 15-Apr-2004))

  subroutine source_properties(zred_now,nz,lifetime2,restart)


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
    if (allocated(rsrcpos)) deallocate(rsrcpos)
    if (allocated(srcMass)) deallocate(srcMass)
    if (allocated(NormFlux)) deallocate(NormFlux)
    if (allocated(NormFluxPL)) deallocate(NormFluxPL)    
    if (allocated(srcSeries)) deallocate(srcSeries)
    
    ! Rank 0 reads in sources
    if (rank == 0) then
       ! Construct the file names
       sourcelistfile=trim(adjustl(dir_src))//"sources.dat"

       call establish_number_of_active_sources (restart)

    endif ! end of rank 0 test

#ifdef MPI
    ! Distribute source number to all other nodes
    call MPI_BCAST(NumSrc,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
             
#ifdef MPILOG
    if (rank /=0) write(logf,*) "Number of sources: ",NumSrc
#endif
    
    ! Allocate arrays for this NumSrc
    if (NumSrc > 0) then
       allocate(srcpos(3,NumSrc))
       allocate(rsrcpos(3,NumSrc))
       allocate(NormFlux(0:NumSrc)) ! 0 will hold lost photons
       allocate(NormFluxPL(0:NumSrc)) ! 0 will hold lost photons       
       allocate(SrcSeries(NumSrc))

       ! Fill in the source arrays
       if (rank == 0) then
          call read_in_sources (restart)
       endif ! of rank 0 test
    
#ifdef MPI
       ! Distribute the source parameters to the other nodes
       call MPI_BCAST(srcpos,3*NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(rsrcpos,3*NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(NormFlux,1+NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(NormFluxPL,NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)      
       call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)       
#endif

       if (rank == 0) then
          write(logf,*) 'Total flux= ',sum(NormFlux(1:NumSrc))*S_star_nominal,' s^-1'
          write(logf,*) 'Total flux PL relative= ',sum(NormFluxPL),' s^-1'
          ! Create array of source numbers for generating random order
          do ns=1,NumSrc
             SrcSeries(ns)=ns
          enddo
          
          ! Make a random order
          call ctrper(SrcSeries(1:NumSrc),1.0)
       endif

#ifdef MPI
       ! Distribute the source series to the other nodes
       call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
      
    else

       if (rank == 0) then
          write(logf,*) "WARNING (sourceprops_test4): No sources found."
       endif

    endif

  end subroutine source_properties

  ! =======================================================================

  subroutine establish_number_of_active_sources (restart)

    integer,intent(in) :: restart

    integer :: srcpos1,srcpos2,srcpos3,sumnumhalo,ioflag
    real :: flux

    open(unit=sourcefile,file=sourcelistfile,status='old')
    
    ! Count number of sources in the source file by looking for EOF
    NumSrc=0
    ioflag=0
    sumnumhalo=0
    do while (ioflag >= 0) ! EOF gives a -1 here. Some of the files have
       ! other errors, which should be ignored.
       read(sourcefile,*,iostat=ioflag) srcpos1,srcpos2,srcpos3,flux
       NumSrc=NumSrc+1
    enddo
    NumSrc=NumSrc-1 ! we counted one too many with the above procedure
    
    ! Report number of sources
    write(logf,"(A,I6)") "Number of sources: ", NumSrc
    write(logf,*) "Flux of last source: ", flux   
    ! Now that we have counted the number of halos, close the file
    close(unit=sourcefile)

  end subroutine establish_number_of_active_sources
  
  ! =======================================================================

  subroutine read_in_sources (restart)

    integer,intent(in) :: restart

    integer :: ns

    open(unit=sourcefile,file=sourcelistfile,status='old')

    ! Read in source positions and mass
    write(logf,*) 'NumSrc', NumSrc
    do ns=1,NumSrc
       read(sourcefile,*) srcpos(1,ns),srcpos(2,ns),srcpos(3,ns),NormFlux(ns),NormFluxPL(ns)

       ! Source is always at cell centre!!
       rsrcpos(1,ns)=x(srcpos(1,ns))
       rsrcpos(2,ns)=y(srcpos(2,ns))
       rsrcpos(3,ns)=z(srcpos(3,ns))
    enddo
    close(sourcefile)
    write(logf,*) "Fluxes: ", NormFlux
  end subroutine read_in_sources

  ! =======================================================================

  subroutine assign_uv_luminosities (lifetime2,nz)
    
    real(kind=dp),intent(in) :: lifetime2 ! time step
    integer,intent(in) :: nz

    integer :: ns

    ! Scale read in photon rates
    do ns=1,NumSrc
       !note that now photons/atom are included in SrcMass
       NormFlux(ns)=NormFlux(ns)/S_star_nominal
    enddo
    
  end subroutine assign_uv_luminosities
  
  ! =======================================================================

  !> Initialization routine: determine the source model and optionally read 
  !! in source properties
  !! Author: Garrelt Mellema
  
  !! Version: empty dummy for compatibility reasons

  subroutine source_properties_ini ()
    
    ! Empty
    
  end subroutine source_properties_ini
  
end module sourceprops
