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
  use file_admin, only: stdinput, logf, file_input
  use cgsconstants, only: m_p
  use astroconstants, only: M_SOLAR, YEAR
  use cosmology_parameters, only: Omega_B, Omega0
  use nbody, only: id_str, M_grid, dir_src, NumZred
  use material, only: xh
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
  use c2ray_parameters, only: phot_per_atom, lifetime, S_star_nominal,&
       StillNeutral, Number_Sourcetypes

  implicit none

  ! Definitions for the source list
  ! These integer variables can be used to make accessing the source
  ! properties from the source list more abstract.
  integer,parameter :: HMACH=4
  integer,parameter :: LMACH=5

#ifdef QUASARS
  integer,parameter :: QSO=6
  !> number of columns in source list with quasars
  integer,parameter,private :: ncolumns_srcfile=6

  !> base name of source list files
  character(len=100),parameter,private :: &
       sourcelistfile_base="_wsubgrid_wquasar_sources.dat"
  character(len=100),parameter,private :: &
       sourcelistfilesuppress_base="_wquasar_sources_used_wfgamma.dat"
#else
  !> number of columns in source list without quasars
  integer,parameter,private :: ncolumns_srcfile=6

  !> base name of source list files
  character(len=100),parameter,private :: &
       sourcelistfile_base="_wsubgrid_sources.dat"
  character(len=100),parameter,private :: &
       sourcelistfilesuppress_base="_sources_used_wfgamma.dat"
#endif

  real(kind=dp),dimension(ncolumns_srcfile),private :: srclist

  !Temp array added by Hannah Ross
  real(kind=dp),dimension(:),allocatable :: temparray

  !> maximum increase in uv to use up cumulated photons
  real(kind=dp),parameter,private :: cumfrac_max=0.15 

  integer :: NumSrc=0 !< Number of sources
  integer :: Prev_NumSrc !< Previous number of sources
  integer,dimension(:,:),allocatable :: srcpos !< mesh position of sources
  real(kind=dp),dimension(:,:),allocatable :: srcMass !< masses of sources 
  real(kind=dp),dimension(:),allocatable :: NormFlux !< normalized ionizing flux of sources
#ifdef PL
  real(kind=dp),dimension(:),allocatable :: NormFluxPL !< normalized ionizing flux of sources
#endif
#ifdef QUASARS
 real(kind=dp),dimension(:),allocatable :: NormFluxQPL !< normalized ionizing flux of sources
 integer :: NumQsr !< counter: number of quasar sources
#endif
  integer,dimension(:),allocatable :: srcSeries  !< a randomized list of sources
  real(kind=dp),dimension(:),allocatable :: uv_array  !< list of UV flux evolution (for some sources models)
  !> The cumulative number of uv photons. We save this number so we can add it
  !! to the uv luminosity the first time sources appear.
  real(kind=dp) :: cumulative_uv=0.0

  character(len=30) :: UV_Model !< type of UV model
  integer :: NumZred_uv !< Number of redshift points in UV model
  integer,private :: NumSrc0=0 !< intermediate source count
  integer,dimension(3),private :: srcpos0
  real(kind=dp),private :: total_SrcMass,NormFluxQPL0!,srcPL
  character(len=6),private :: z_str !< string value of redshift
  integer,private :: NumMassiveSrc !< counter: number of massive sources
  integer,private :: NumSupprbleSrc !< counter: number of suppressible sources
  integer,private :: NumSupprsdSrc !< counter: number of suppressed sources
  integer,private :: NumPlSrc !< counter: number of power law sources
  real(kind=dp),private :: cumfrac

  character(len=512),private :: sourcelistfile,sourcelistfilesuppress

#ifdef MPI
    integer :: mympierror
#endif
 
contains
  
  ! =======================================================================

  !> Set the source properties for this redshift
  subroutine source_properties(zred_now,nz,lifetime2,restart)
    
    ! Input routine: establish the source properties
    
    real(kind=dp),intent(in) :: zred_now ! current redshift
    real(kind=dp),intent(in) :: lifetime2 ! time step
    integer,intent(in) :: nz
    integer,intent(in) :: restart
    
    integer :: ns,ns0
    
#ifdef MPILOG     
    write(logf,*) "Check sourceprops: ",zred_now,nz,lifetime2,restart
#endif 
    
    ! Deallocate source arrays
    if (allocated(srcpos)) deallocate(srcpos)
    if (allocated(srcMass)) deallocate(srcMass)
    if (allocated(NormFlux)) deallocate(NormFlux)
#ifdef PL
    if (allocated(NormFluxPL)) deallocate(NormFluxPL)
#endif
#ifdef QUASARS
    if (allocated(NormFluxQPL)) deallocate(NormFluxQPL)
#endif

#ifdef MPILOG     
    write(logf,*) "deallocated all"
#endif 

    ! Save previous number of sources
    Prev_NumSrc=NumSrc
    
    ! Rank 0 reads in sources
    if (rank == 0) then
       
       ! Sources are read from files with redshift in the file name:
       ! construct redshift string
       write(z_str,'(f6.3)') zred_now
       
       ! Construct the file names
       sourcelistfile=trim(adjustl(dir_src))//&       
            trim(adjustl(z_str))//"-"//trim(adjustl(id_str))// &
            trim(adjustl(sourcelistfile_base))
       sourcelistfilesuppress=trim(adjustl(dir_src))//&
            trim(adjustl(z_str))//"-"//trim(adjustl(id_str))// &
            trim(adjustl(sourcelistfilesuppress_base))
       
       call establish_number_of_active_sources (restart)
       
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
       allocate(srcpos(3,NumSrc))
       allocate(SrcMass(NumSrc,0:Number_Sourcetypes))
       allocate(NormFlux(0:NumSrc)) ! 0 will hold lost photons
#ifdef PL
       allocate(NormFluxPL(NumSrc))
#endif
#ifdef QUASARS
       allocate(NormFluxQPL(NumSrc))
#endif

#ifdef MPILOG     
       write(logf,*) "allocated all" 
#endif

       ! Fill in the source arrays
       if (rank == 0) then
          call read_in_sources (restart)
          
          ! Set cumulative number of uv photons to zero if this is not the
          ! first redshift for which sources are active (this way cumulative_uv
          ! can be used in all cases).
          ! New version: cumulative_uv is slowly reduced
          !if (Prev_NumSrc /= 0) cumulative_uv=0.0
          call assign_uv_luminosities (lifetime2,nz)
          
          write(logf,*) 'Source lifetime=', lifetime2/(1e6*YEAR),' Myr'
          write(logf,*) 'Total photon rate (BB)= ', &
               sum(NormFlux)*S_star_nominal,' s^-1'
#ifdef PL
          write(logf,*) 'Total photon rate (PL)= ', &
               sum(NormFluxPL)*pl_S_star_nominal,' s^-1'
#endif
#ifdef QUASARS
          write(logf,*) 'Total photon rate (Q)= ', &
               sum(NormFluxQPL)*qpl_S_star_nominal,' s^-1'
#endif
          
       endif ! of rank 0 test

#ifdef MPI
       ! Distribute the source parameters to the other nodes
       call MPI_BCAST(srcpos,3*NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(SrcMass,(1+Number_Sourcetypes)*NumSrc, &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(NormFlux,NumSrc+1, &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#ifdef PL
       call MPI_BCAST(NormFluxPL,NumSrc, &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)     
#endif
#ifdef QUASARS
       call MPI_BCAST(NormFluxQPL,NumSrc, &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)    

#endif
#endif
    
       ! Randomize the source list
       !call create_random_source_list()

    else
    
       ! For the prescribed uv evolution model, accumulate the photons
       ! (to be used once the first sources appear)
       if (UV_Model == "Fixed N_gamma") &
            cumulative_uv=cumulative_uv + uv_array(nz)
       
    endif
    
  end subroutine source_properties

  ! =======================================================================

  subroutine establish_number_of_active_sources (restart)

    integer,intent(in) :: restart
    
    real :: odens	
    
    integer :: ns0
    integer :: ns


    if (restart == 0 .or. restart == 1) then
       open(unit=50,file=sourcelistfile,status='old')
       ! Number of sources
       read(50,*) NumSrc0
       
       ! Report
       write(logf,*) "Total number of source locations, no suppression: ", &
            NumSrc0
       
       ! Read in source positions and mass to count the number
       ! of non-suppressed sources
       NumSrc = 0
       NumMassiveSrc = 0
       NumSupprbleSrc = 0
       NumSupprsdSrc = 0
#ifdef QUASARS
       NumQsr = 0
#endif
       do ns0=1,NumSrc0
          ! If you change the following lines, also change it below in
          ! read_in_sources
          read(50,*) srclist(1:ncolumns_srcfile)
          srcpos0(1:3)=int(srclist(1:3))

          ! Massive sources are never suppressed.
#ifdef QUASARS
          if (srclist(HMACH) /= 0.0 .or. srclist(QSO) /= 0.0) then
#else
          if (srclist(HMACH) /= 0.0) then
#endif
             NumSrc=NumSrc+1
          ! if the cell is still neutral, no suppression (if we use the Iliev
          ! et al source model)   
          elseif (xh(srcpos0(1),srcpos0(2),srcpos0(3),1) < StillNeutral .and. &
               UV_Model == "Iliev et al") then
             NumSrc=NumSrc+1
          endif

          ! Count different types of sources
          if (srclist(HMACH) /= 0.0) NumMassiveSrc=NumMassiveSrc+1
          if (srclist(LMACH) /= 0.0) NumSupprbleSrc=NumSupprbleSrc+1
#ifdef QUASARS
          if (srclist(QSO) /= 0.0) NumQsr = NumQsr+1 
#endif
          ! How many suppressed?
          if (srclist(LMACH) /= 0.0) then
             if (xh(srcpos0(1),srcpos0(2),srcpos0(3),1) > StillNeutral .or. &
               UV_Model /= "Iliev et al") NumSupprsdSrc=NumSupprsdSrc+1
          endif
       enddo
       close(50)

       ! Report to log file
       write(logf,*) "Number of suppressable sources: ",NumSupprbleSrc
       write(logf,*) "Number of suppressed sources: ",NumSupprsdSrc
       write(logf,*) "Number of massive sources: ",NumMassiveSrc
#ifdef QUASARS
       write(logf,*) "Number of quasar sources: ",NumQsr
#endif
       if (NumSupprbleSrc > 0) write(logf,*) "Suppressed fraction: ", &
            real(NumSupprsdSrc)/real(NumSupprbleSrc)
    else
       ! Upon restart from intermediate redshift use the previously 
       ! calculated suppressed source list
       open(unit=49,file=sourcelistfilesuppress,status='unknown')
       ! Number of sources
       read(49,*) NumSrc
       close(49)
    endif
    write(logf,*) "Number of sources, with suppression: ",NumSrc

  end subroutine establish_number_of_active_sources

  ! =======================================================================

  subroutine read_in_sources (restart)

    integer,intent(in) :: restart

    integer :: ns
    integer :: ns0
    character(len=24) :: format_str  !string for formatting output 

    if (restart == 0 .or. restart == 1) then
       open(unit=50,file=sourcelistfile, status='old')
       ! Number of sources
       read(50,*) NumSrc0
       ! Read in source positions and mass
       ns=0
       do ns0=1,NumSrc0
          read(50,*) srclist(1:ncolumns_srcfile)
          srcpos0(1:3)=int(srclist(1:3))

          ! Process the source list entry through the suppression
          ! algorithm
          if (xh(srcpos0(1),srcpos0(2),srcpos0(3),1) < StillNeutral) then
#ifdef QUASARS
             if (UV_Model == "Iliev et al" .or. srclist(HMACH) > 0.0d0 .or. &
                 srclist(QSO) > 0.0d0) then
#else
             if (UV_Model == "Iliev et al" .or. srclist(HMACH) > 0.0d0) then
#endif
                ! the cell is still neutral, no suppression
                ns=ns+1
                ! Source positions in file start at 1!
                srcpos(1,ns)=srcpos0(1)
                srcpos(2,ns)=srcpos0(2)
                srcpos(3,ns)=srcpos0(3)
                SrcMass(ns,1)=srclist(HMACH)
                if (UV_Model == "Iliev et al") then
                   SrcMass(ns,2)=srclist(LMACH)
                else
                   SrcMass(ns,2)=0.0
                endif
#ifdef QUASARS
                NormFluxQPL(ns)=srclist(QSO) ! set QSO source flux 
                ! (here it is assumed they are always active)
#endif
             endif
#ifdef QUASARS
          elseif (srclist(HMACH) > 0.0d0 .or. srclist(QSO) > 0.0d0) then
#else
          elseif (srclist(HMACH) > 0.0d0) then
#endif
             !the cell is ionized but source is massive enough to survive
             !and is assumed Pop. II 
             ns=ns+1
             ! Source positions in file start at 1!
             srcpos(1,ns)=srcpos0(1)
             srcpos(2,ns)=srcpos0(2)
             srcpos(3,ns)=srcpos0(3)
             SrcMass(ns,1)=srclist(HMACH)
             SrcMass(ns,2)=0.0
#ifdef QUASARS
             NormFluxQPL(ns)=srclist(QSO) ! set QSO source flux 
             ! (here it is assumed they are always active)
#endif
          endif

       enddo
       
       ! Collect total source mass (weigthed with efficiency factor
       ! in case of the Iliev et al source model).
       if (UV_Model == "Iliev et al") then
          SrcMass(:,0)=SrcMass(:,1)*phot_per_atom(1)  & !massive sources
               +SrcMass(:,2)*phot_per_atom(2)      !small sources  
#ifdef PL
          NormFluxPL(:)=xray_phot_per_atom*(SrcMass(:,1)+SrcMass(:,2))
#endif
       else
          SrcMass(:,0)=SrcMass(:,1)!+SrcMass(:,2)
#ifdef PL
          NormFluxPL(:)=xray_phot_per_atom*SrcMass(:,1)
#endif
       endif

       ! Save new, processed, source list (without the suppressed sources)
       call save_source_list()

    else ! of restart test

       ! Read source list from file saved previously
       open(unit=49,file=sourcelistfilesuppress,status="old")
       write(logf,*) "Reading ",NumSrc," sources from ", &
            trim(adjustl(sourcelistfilesuppress))
       read(49,*) NumSrc

       ! Data in saved source list depends on what sources were
       ! recorded / are being used.
#if defined(QUASARS) && defined(PL)
         if (allocated(temparray)) deallocate(temparray)
          allocate(temparray(3))
#elif defined(QUASARS)
          if (allocated(temparray)) deallocate(temparray)
          allocate(temparray(2))   
#elif defined(PL)
          if (allocated(temparray)) deallocate(temparray)
          allocate(temparray(2))
#else
          if (allocated(temparray)) deallocate(temparray)
          allocate(temparray(1))
#endif

       do ns0=1,NumSrc
          
          read(49,*) srcpos(1,ns0),srcpos(2,ns0),srcpos(3,ns0),&
               temparray
          
          SrcMass(ns0,0) = temparray(1)
#if defined(QUASARS) && defined(PL)
          NormFluxPL(ns0) = temparray(2)
          NormFluxQPL(ns0) = temparray(3)
#elif defined(QUASARS)
          NormFluxQPL(ns0) = temparray(2)
#elif defined(PL)
          NormFluxPL(ns0) = temparray(2)
#endif
       enddo
       deallocate(temparray)
       close(49)
    endif ! of restart test
    
  end subroutine read_in_sources

  ! =======================================================================

  subroutine save_source_list ()

    integer :: ns0
    character(len=24) :: format_str  !string for formatting output 

    ! Open the processed source list file
    open(unit=49,file=sourcelistfilesuppress,status='unknown')
    ! Write number of active sources positions
    write(49,*) NumSrc

    ! Write out the source information
    do ns0=1,NumSrc
       
       ! Content depends on what sources are being used
#if defined(QUASARS) && defined(PL)
       format_str = "(3i4,2f10.3,1es15.3)"!,1e10.3)"
       if (allocated(temparray)) deallocate(temparray)
       allocate(temparray(3))
       temparray(1) = SrcMass(ns0,0)
       temparray(2) = NormFluxPL(ns0)
       temparray(3) = NormFluxQPL(ns0)
#elif defined(QUASARS)
       format_str = "(3i4,1f10.3,1es15.3)"
       if (allocated(temparray)) deallocate(temparray)
       allocate(temparray(2))   
       temparray(1) = SrcMass(ns0,0)
       temparray(2) = NormFluxQPL(ns0)
#elif defined(PL)
       format_str = "(3i4,2f10.3)"
       if (allocated(temparray)) deallocate(temparray)
       allocate(temparray(2))
       temparray(1) = SrcMass(ns0,0)
       temparray(2) = NormFluxPL(ns0)
#else
       format_str = "(3i4,1f10.3)"
       if (allocated(temparray)) deallocate(temparray)
       allocate(temparray(1))
       temparray(1) = SrcMass(ns0,0)
#endif
       
       ! Write the source position information
       write(49,format_str) srcpos(1,ns0),srcpos(2,ns0),srcpos(3,ns0), &
            temparray(1:size(temparray))
       
    enddo
    
    ! Close the file
    close(49)
    
  end subroutine save_source_list

  ! =======================================================================

  subroutine assign_uv_luminosities (lifetime2,nz)
    
    real(kind=dp),intent(in) :: lifetime2 ! time step
    integer,intent(in) :: nz

    integer :: ns

#ifdef MPILOG
    if (rank /=0) then
       write(logf,*) 'Source lifetime=', lifetime2/(1e6*YEAR),' Myr'
       write(logf,*) 'S_star_nominal=',S_star_nominal
#ifdef PL
       write(logf,*) 'pl_S_star_nominal=',pl_S_star_nominal
#endif
#ifdef QUASARS
       write(logf,*) 'qpl_S_star_nominal=',qpl_S_star_nominal
#endif
    endif
#endif
    
    ! Turn masses into luminosities
    select case (UV_Model)

    case ("Iliev et al")
       do ns=1,NumSrc
          ! Calculate normalized luminosity of stellar sources
          ! note that escaping photons/atom are included in SrcMass
          NormFlux(ns)=Luminosity_from_mass(SrcMass(ns,0))/lifetime2

          ! Calculate normalized luminosity of power law source
#ifdef PL
          NormFluxPL(ns)=PL_Luminosity_from_mass(NormFluxPL(ns))/lifetime2
#endif
#ifdef QUASARS
          NormFluxQPL(ns)=QPL_Luminosity_convert(NormFluxQPL(ns)) !NormFluxQPL(ns)/(qpl_S_star_nominal)
#endif
       enddo

    case ("Fixed N_gamma")
       if (nz <= NumZred_uv) then
          cumfrac=min(cumfrac_max,cumulative_uv/uv_array(nz))
          if (rank == 0) then 
             write(logf,*) 'Cumulative versus current photons: ', &
                  cumulative_uv,uv_array(nz),cumulative_uv/uv_array(nz)
             write(logf,*) 'Cumulative fraction used: ', cumfrac
          endif
          total_SrcMass=sum(SrcMass(:,0))
          ! Only set NormFlux when data is available!
          do ns=1,NumSrc
             NormFlux(ns)=(1.0+cumfrac)*uv_array(nz)/lifetime2*SrcMass(ns,0)/ &
                  (total_SrcMass*S_star_nominal)
             ! Calculate normalized luminosity of power law source
#ifdef PL
             NormFluxPL(ns)=PL_Luminosity_from_mass(NormFluxPL(ns))/lifetime2
#endif
#ifdef QUASARS
             NormFluxQPL(ns)=QPL_Luminosity_convert(NormFluxQPL(ns)) !NormFluxQPL(ns)/(qpl_S_star_nominal)
#endif
          enddo
          ! Subtract extra photons from cumulated photons
          cumulative_uv=max(0.0_dp,cumulative_uv-cumfrac*uv_array(nz))
          !write(logf,*) uv_array(nz),SrcMass(:,0),uv_array(nz)/lifetime2
       else
          ! For high redshifts there may not be a uv model.
          ! Set fluxes to zero.
          NormFlux(:)=0.0
          if (rank == 0) write(logf,*) &
               "No UV model available, setting fluxes to zero."
       endif

    case ("Fixed Ndot_gamma")
       if (nz <= NumZred_uv) then
          total_SrcMass=sum(SrcMass(:,0))
          ! Only set NormFlux when data is available!
          do ns=1,NumSrc
             NormFlux(ns)=uv_array(nz)* &
                  SrcMass(ns,0)/(total_SrcMass*S_star_nominal)
             ! Calculate normalized luminosity of power law source
#ifdef PL
             NormFluxPL(ns)=PL_Luminosity_from_mass(NormFluxPL(ns))/lifetime2
#endif
#ifdef QUASARS
             NormFluxQPL(ns)=QPL_Luminosity_convert(NormFluxQPL(ns)) !NormFluxQPL(ns)/(qpl_S_star_nominal)
#endif

          enddo
       else
          NormFlux(:)=0.0
          if (rank == 0) write(logf,*) "No UV model available, setting fluxes to zero."
       endif
    end select
    
#ifdef MPILOG
    if (rank /=0) then
       write(logf,*) 'Assigned UV luminosities to sources'
       write(logf,*) 'UV_Model=',UV_Model
    endif
#endif
    
  end subroutine assign_uv_luminosities
  
  ! =======================================================================
  
  function Luminosity_from_mass (Mass)

    ! The normalized flux (luminosity) for normal sources is expressed
    ! in terms of a standard ionizing photon rate, called S_star_nominal.
    ! In radiation the tables have been calculated for a spectrum
    ! with an ionizing photon rate of S_star_nominal.

    ! Mass is supposed to be total mass of the source MULTIPLIED with
    ! the efficiency factor (f) which is the product of the star formation
    ! fraction, the escape fraction and the number of photons produced
    ! per baryon. Because of the latter factor we need to convert the
    ! mass to number of baryons by dividing my the proton mass m_p.
    ! NOTE: number of baryons is the total number of nucleons, which
    ! is why we divide my m_p and NOT by mu*m_p (where mu is mean mass
    ! of atoms/ions).

    real(kind=dp) :: Mass !< mass in units of grid masses
    real(kind=dp) :: Luminosity_from_mass

    Luminosity_from_mass = Mass*M_grid*Omega_B/(Omega0*m_p)/S_star_nominal

  end function Luminosity_from_mass

  ! =======================================================================
#ifdef PL
  function PL_Luminosity_from_mass (Mass)

    ! The normalized flux (luminosity) for normal sources is expressed
    ! in terms of a standard ionizing photon rate, called S_star_nominal.
    ! In radiation the tables have been calculated for a spectrum
    ! with an ionizing photon rate of S_star_nominal.

    ! Mass is supposed to be total mass of the source MULTIPLIED with
    ! the efficiency factor (f) which is the product of the star formation
    ! fraction, the escape fraction and the number of photons produced
    ! per baryon. Because of the latter factor we need to convert the
    ! mass to number of baryons by dividing my the proton mass m_p.
    ! NOTE: number of baryons is the total number of nucleons, which
    ! is why we divide my m_p and NOT by mu*m_p (where mu is mean mass
    ! of atoms/ions).

    ! Note: this model assumes that the x-ray luminosity of a halo
    ! is proportional to its SFR.

    real(kind=dp) :: Mass !< mass in units of grid masses
    real(kind=dp) :: PL_Luminosity_from_mass

    PL_Luminosity_from_mass = Mass*M_grid*Omega_B/(Omega0*m_p)/pl_S_star_nominal

  end function PL_Luminosity_from_mass
#endif
  ! =======================================================================
#ifdef QUASARS
  function QPL_Luminosity_convert (Lum)

    ! The normalized flux (luminosity) for normal sources is expressed
    ! in terms of a standard ionizing photon rate, called qpl_S_star_nominal.
    ! In radiation the tables have been calculated for a spectrum
    ! with an ionizing photon rate of qpl_S_star_nominal.

    ! The luminosity is given in ergs/second at 2kev, needs to be converted
    ! to the total number of photons per second.

    use cgsconstants, only: ev2fr, ev2erg

    real(kind=dp) :: Lum
    real(kind=dp) :: Emax
    real(kind=dp) :: Emin
    real(kind=dp) :: delta_E
    real(kind=dp) :: QPL_Luminosity_convert
    real(kind=dp) :: alpha
    
    !Convert luminosity to total number of photons per second
    Emin = qpl_MinFreq_nominal/(ev2fr) !in ev
    Emax = qpl_MaxFreq_nominal/(ev2fr) !in ev
    delta_E = (Emax - Emin)*ev2erg !Needs to be in ergs as luminosity is given in ergs

    alpha = qpl_index_nominal - 1.0

    QPL_Luminosity_convert = -1.0/delta_E * Lum/(2000**(-alpha)) * &
                1.0/alpha*(Emax**(-alpha) - Emin**(-alpha))
!    write(*,*) "Emax, Emin: ", Emax, Emin    
!    write(*,*) "Lum: ", Lum
!    if (QPL_Luminosity_convert /= 0.0) write(*,*) "number of photons: ", QPL_Luminosity_convert

    !Normalise
    QPL_Luminosity_convert = QPL_Luminosity_convert/qpl_S_star_nominal

  end function QPL_Luminosity_convert
#endif
  ! =======================================================================

  !> Initialization routine: determine the source model and optionally read 
  !! in source properties
  !! Author: Garrelt Mellema
  
  
  !! This accomodates different source models
  !! 0: Iliev et al source, Ndot_gamma= f*M_halo/timestep
  !! 1: Fixed total N_gamma, still need to divide by time step
  !! 2: Fixed total Ndot_gamma.
  subroutine source_properties_ini ()
    

    integer :: uv_answer
    real(kind=dp) :: z_in, N_source_nosupp, N_source_supp, N_gamma_nosupp
    character(len=180) :: uv_file ! name of file with uv model for redshifts
    integer :: nz

    ! Ask for redshift file
    if (rank == 0) then
       if (.not.file_input) write(*,"(A,$)") "UV Luminosity recipe (0,1,2): "
       read(stdinput,*) uv_answer
       select case (uv_answer)
       case(0)
          UV_Model = "Iliev et al"
       case(1)
          UV_Model = "Fixed N_gamma"
       case(2)
          UV_Model = "Fixed Ndot_gamma"
       end select
       
       if (uv_answer > 0) then
          if (.not.file_input) write(*,"(A,$)") "File with UV data: "
          read(stdinput,*) uv_file
          
          ! Open and read redshift file
          open(unit=60,file=uv_file,form="formatted",status="old")
          read(unit=60,fmt=*) NumZred_uv
          if (NumZred_uv /= NumZred) then
             write(logf,*) "WARNING: Number of redshifts in UV luminosity file (", &
                  NumZred_uv,") does not match number of redshifts in ", &
                  "redshift file (",NumZred,")."
          endif
          allocate(uv_array(NumZred_uv))
          if (uv_answer == 1) then
             do nz=1,NumZred_uv
                read(unit=60,fmt=*) z_in, N_source_nosupp, N_source_supp, & 
                     N_gamma_nosupp, uv_array(nz)
             enddo
          else
             do nz=1,NumZred_uv
                read(unit=60,fmt=*) z_in, uv_array(nz)
             enddo
          endif
          close(60)
       endif
    endif
#ifdef MPI
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(uv_answer,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(UV_Model,30,MPI_CHARACTER,0,MPI_COMM_NEW,mympierror)
    if (uv_answer > 0) then
       call MPI_BCAST(NumZred_uv,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)    
       if (rank /= 0) allocate(uv_array(NumZred_uv))
       call MPI_BCAST(uv_array,NumZred_uv,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,&
            mympierror)
    endif
#endif
    
  end subroutine source_properties_ini

  ! =======================================================================
  
  subroutine create_random_source_list()

    ! For random permutation of sources
    use  m_ctrper
    
    integer :: ns

    ! First deallocate the array (if it was allocated)
    if (allocated(srcSeries)) deallocate(srcSeries)

    ! Allocate the array for the current number of sources
    allocate(SrcSeries(NumSrc))

    if (rank == 0) then
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

  end subroutine create_random_source_list

end module sourceprops
