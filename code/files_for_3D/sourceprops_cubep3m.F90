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
  use c2ray_parameters, only: phot_per_atom, xray_phot_per_atom, lifetime, &
       S_star_nominal, pl_S_star_nominal, StillNeutral, Number_Sourcetypes
  use c2ray_parameters, only: EddLeff_nominal,mass_nominal,EddLum

  implicit none

  !> base name of source list files
  character(len=100),parameter,private :: &
       sourcelistfile_base="_sources.dat" ! "_wsubgrid_sources.dat"
  character(len=100),parameter,private :: &
       sourcelistfilesuppress_base="_sources_used_wfgamma.dat"

  !> number of columns in source list
  integer,parameter,private :: ncolumns_srcfile=5
  real,dimension(ncolumns_srcfile),private :: srclist

  !> maximum increase in uv to use up cumulated photons
  real(kind=dp),parameter,private :: cumfrac_max=0.15 

  integer :: NumSrc=0 !< Number of sources
  integer :: Prev_NumSrc !< Previous number of sources
  integer,dimension(:,:),allocatable :: srcpos !< mesh position of sources
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
  real(kind=dp),private :: srcMass00,srcMass01,total_SrcMass,srcPL
  character(len=6),private :: z_str !< string value of redshift
  integer,private :: NumMassiveSrc !< counter: number of massive sources
  integer,private :: NumSupprbleSrc !< counter: number of suppressible sources
  integer,private :: NumSupprsdSrc !< counter: number of suppressed sources
  integer,private :: NumPlSrc !< counter: number of power law sources
  real(kind=dp),private :: cumfrac

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
    if (allocated(srcMass)) deallocate(srcMass)
    if (allocated(NormFlux)) deallocate(NormFlux)
    if (allocated(NormFluxPL)) deallocate(NormFluxPL)
    if (allocated(srcSeries)) deallocate(srcSeries)
    
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
       allocate(NormFluxPL(NumSrc))
       allocate(SrcSeries(NumSrc))
       
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
          write(logf,*) 'Total photon rate (PL)= ', &
               sum(NormFluxPL)*pl_S_star_nominal,' s^-1'
          
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
       call MPI_BCAST(SrcMass,(1+Number_Sourcetypes)*NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(NormFlux,NumSrc+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(NormFluxPL,NumSrc,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)       
       ! Distribute the source series to the other nodes
       call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
    
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
       do ns0=1,NumSrc0
          ! If you change the following lines, also change it below in          
          ! read_in_sources                                                     
          read(50,*) srclist(1:ncolumns_srcfile)
          srcpos0(1:3)=int(srclist(1:3))
          srcMass00=srclist(4) !massive sources (HMACHs)                        
          srcMass01=srclist(5) !low-mass sources (LMACHs)                       
          !read(50,*) srcpos0(1),srcpos0(2),srcpos0(3),SrcMass00,SrcMass01
          ! Massive sources are never suppressed.
          if (SrcMass00 /= 0.0) then
             NumSrc=NumSrc+1
          ! if the cell is still neutral, no suppression (if we use the Iliev
          ! et al source model)   
          elseif (xh(srcpos0(1),srcpos0(2),srcpos0(3),1) < StillNeutral .and. &
               UV_Model == "Iliev et al") then
             NumSrc=NumSrc+1
          endif

          ! Count different types of sources
          if (SrcMass00 /= 0.0) NumMassiveSrc=NumMassiveSrc+1
          if (SrcMass01 /= 0.0) NumSupprbleSrc=NumSupprbleSrc+1
          ! How many suppressed?
          if (SrcMass01 /= 0.0) then
             if (xh(srcpos0(1),srcpos0(2),srcpos0(3),1) > StillNeutral .or. &
               UV_Model /= "Iliev et al") NumSupprsdSrc=NumSupprsdSrc+1
          endif
       enddo
       close(50)
       write(logf,*) "Number of suppressable sources: ",NumSupprbleSrc
       write(logf,*) "Number of suppressed sources: ",NumSupprsdSrc
       write(logf,*) "Number of massive sources: ",NumMassiveSrc
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

    if (restart == 0 .or. restart == 1) then
       open(unit=50,file=sourcelistfile,status='old')
       ! Number of sources
       read(50,*) NumSrc0
       ! Read in source positions and mass
       ns=0
       do ns0=1,NumSrc0
          read(50,*) srclist(1:ncolumns_srcfile)
          srcpos0(1:3)=int(srclist(1:3))
          srcMass00=srclist(4) !massive sources (HMACHs)                       
          srcMass01=srclist(5) !low-mass sources (LMACHs)                       
          !read(50,*) srcpos0(1),srcpos0(2),srcpos0(3), &
          !     SrcMass00,SrcMass01
          
          if (xh(srcpos0(1),srcpos0(2),srcpos0(3),1) < StillNeutral) then
             if (UV_Model == "Iliev et al" .or. SrcMass00 > 0.0d0) then
                ! the cell is still neutral, no suppression
                ns=ns+1
                ! Source positions in file start at 1!
                srcpos(1,ns)=srcpos0(1)
                srcpos(2,ns)=srcpos0(2)
                srcpos(3,ns)=srcpos0(3)
                SrcMass(ns,1)=SrcMass00
                if (UV_Model == "Iliev et al") then
                   SrcMass(ns,2)=SrcMass01
                else
                   SrcMass(ns,2)=0.0
                endif
             endif
          elseif (SrcMass00 > 0.0d0) then
             !the cell is ionized but source is massive enough to survive
             !and is assumed Pop. II 
             ns=ns+1
             ! Source positions in file start at 1!
             srcpos(1,ns)=srcpos0(1)
             srcpos(2,ns)=srcpos0(2)
             srcpos(3,ns)=srcpos0(3)
             SrcMass(ns,1)=SrcMass00
             SrcMass(ns,2)=0.0
          endif
       enddo
       
       ! Collect total source mass (weigthed with efficiency factor
       ! in case of the Iliev et al source model).
       if (UV_Model == "Iliev et al") then
          SrcMass(:,0)=SrcMass(:,1)*phot_per_atom(1)  & !massive sources
               +SrcMass(:,2)*phot_per_atom(2)      !small sources  
          NormFluxPL(:)=xray_phot_per_atom*(SrcMass(:,1)+SrcMass(:,2))
       else
          SrcMass(:,0)=SrcMass(:,1)!+SrcMass(:,2)
          NormFluxPL(:)=xray_phot_per_atom*SrcMass(:,1)
       endif
       ! Save new source list, without the suppressed ones
       open(unit=49,file=sourcelistfilesuppress,status='unknown')
       write(49,*) NumSrc
       do ns0=1,NumSrc
          write(49,"(3i4,2f10.3)") srcpos(1,ns0),srcpos(2,ns0),srcpos(3,ns0), &
               SrcMass(ns0,0),NormFluxPL(ns0)
       enddo
       close(49)
    else ! of restart test
       ! Read source list from file saved previously
       open(unit=49,file=sourcelistfilesuppress,status="old")
       write(logf,*) "Reading ",NumSrc," sources from ", &
            trim(adjustl(sourcelistfilesuppress))
       read(49,*) NumSrc
       do ns0=1,NumSrc
          read(49,*) srcpos(1,ns0),srcpos(2,ns0),srcpos(3,ns0), &
               SrcMass(ns0,0),NormFluxPL(ns0)
       enddo
       close(49)
    endif ! of restart test
    
  end subroutine read_in_sources

  ! =======================================================================

  subroutine assign_uv_luminosities (lifetime2,nz)
    
    real(kind=dp),intent(in) :: lifetime2 ! time step
    integer,intent(in) :: nz

    integer :: ns

#ifdef MPILOG
    if (rank /=0) then
       write(logf,*) 'Source lifetime=', lifetime2/(1e6*YEAR),' Myr'
       write(logf,*) 'S_star_nominal=',S_star_nominal
       write(logf,*) 'pl_S_star_nominal=',pl_S_star_nominal
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
          NormFluxPL(ns)=PL_Luminosity_from_mass(NormFluxPL(ns))/lifetime2
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
             NormFluxPL(ns)=PL_Luminosity_from_mass(NormFluxPL(ns))/lifetime2
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
             NormFluxPL(ns)=PL_Luminosity_from_mass(NormFluxPL(ns))/lifetime2
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

#ifdef MPI
    integer :: mympierror
#endif

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
  
end module sourceprops
