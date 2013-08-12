!     This module contains data and routines which deal with radiative
!     effects. Its main part deal with photo-ionizing radiation, but it
!     also initializes other radiative properties, such as cooling (which
!     are contained in different modules).
!     It can be used in hydrodynamic or stand-alone radiative transfer 
!     calculations.

module radiation
  
  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use mathconstants, only: pi
  use cgsconstants, only: sigma_SB, &                    ! Stefan-Boltzmann constant
                          hplanck, &                     ! Planck constant
                          k_B, &                         ! Boltzmann constant
                          two_pi_over_c_square, &        ! two times pi over c aquare
                          ev2fr                          ! eV to Hz conversion
  use cgsphotoconstants, only: ion_freq_HI,&             ! HI ionization energy in frequency
                               ion_freq_HeI,&            ! HeI ionization energy in frequency
                               ion_freq_HeII,&           ! HeII ionization energy in frequency
                               sigma_HI_at_ion_freq,&    ! HI cross section at its ionzing frequency
                               sigma_HeI_at_ion_freq,&   ! HeI cross section at its ionzing frequency
                               sigma_HeII_at_ion_freq    ! HeII cross section at its ionzing frequency
  use astroconstants, only: R_SOLAR, &                   ! Solar radius
                            L_SOLAR                      ! Solar luminosity
  use romberg, only: scalar_romberg, &                   ! 1D integration function
                     vector_romberg, &                   ! 1D integration subroutine
                     romberg_initialisation              ! Romberg initialisation procedure
  use c2ray_parameters, only: T_eff_nominal,&            ! Black body  effective temperature for for nominal SED
                              S_star_nominal, &          ! Ionizing photon rate for for nominal SED
                              pl_index_nominal,&         ! Power law index for for nominal SED
                              EddLeff_nominal,&          ! Eddington efficiency for for nominal SED
                              EddLum                     ! Eddington luminosity for for nominal SED
  use material, only: isothermal

  implicit none

  integer,parameter :: NumFreq = 512      ! Number of integration points in each of the frequency bins
  integer,parameter :: NumTau = 2000      ! Number of table points for the optical depth
  integer,parameter :: NumBndin1 = 1      ! Number of frequency sub-bins in interval 1 
  integer,parameter :: NumBndin2 = 26     ! Number of frequency sub-bins in interval 2
  integer,parameter :: NumBndin3 = 20     ! Number of frequency sub-bins in interval 3
  integer,parameter :: NumFreqBnd=NumBndin1+NumBndin2+NumBndin3       ! Total number of frequency bins
  integer,parameter :: NumheatBin=NumBndin1+NumBndin2*2+NumBndin3*3   ! Total number of heating bins 

  ! Optical depths at the entrance of the grid.
  ! It can be used if radiation enters the simulation volume from the outside.
  real(kind=dp) :: boundary_tauHI = 0.0
  real(kind=dp) :: boundary_tauHeI = 0.0
  real(kind=dp) :: boundary_tauHeII = 0.0

  ! Parameters defining the optical depth entries in the table.
  real(kind=dp),parameter :: minlogtau = -20.0                             ! Table position starts at log10(minlogtau) 
  real(kind=dp),parameter :: maxlogtau = 4.0                               ! Table position ends at log10(maxlogtau) 
  real(kind=dp),parameter :: dlogtau = (maxlogtau-minlogtau)/real(NumTau)  ! dlogtau is the step size in log10(tau)

  ! Variables connected to secondary ionizations.
  ! in general, I'm following Ricotti et al 2002

  real(kind=dp), dimension(1:3) :: CR1=(/0.3908_dp, 0.0554_dp, 1.0_dp/)
  real(kind=dp), dimension(1:3) :: bR1=(/0.4092_dp, 0.4614_dp, 0.2663_dp/)
  real(kind=dp), dimension(1:3) :: dR1=(/1.7592_dp, 1.6660_dp, 1.3163_dp/)
    
  real(kind=dp), dimension(1:3) :: CR2=(/0.6941_dp,0.0984_dp,3.9811_dp/)
  real(kind=dp), dimension(1:3) :: aR2=(/0.2_dp,0.2_dp,0.4_dp/)
  real(kind=dp), dimension(1:3) :: bR2=(/0.38_dp,0.38_dp,0.34_dp/)
  !real(kind=dp), dimension(1:3) :: dR2=(/2.0_dp,2.0_dp,2.0_dp/) write explicitly ^2 -> introduce xeb

  ! Number of photon of the power law source
  real(kind=dp) :: pl_input_flux = 0.0

  ! Logical that determines the use of grey opacities
  logical,parameter :: grey = .false. 

  ! Type of source, B=black body, P=power law source
  Character :: sourcetype = " "

  ! Stellar properties
  real(kind=dp) :: T_eff  = 0.0      ! Black body effective temperature
  real(kind=dp) :: R_star = 0.0      ! Black body radius
  real(kind=dp) :: L_star = 0.0      ! Black body luminosity
  real(kind=dp) :: L_star_ion = 0.0 ! Black body ionizing luminosity
  real(kind=dp) :: S_star = 0.0      ! Black body ionizing photons rate
  real(kind=dp) :: R_star2 = 0.0     ! Square of R_star
  real(kind=dp) :: h_over_kT    ! Planck constant over k_B * T_eff

  ! The lowest and highest frequency subbands used for the bb and pl source
  integer :: bb_FreqBnd_UpperLimit=NumFreqBnd
  integer :: pl_FreqBnd_LowerLimit
  integer :: pl_FreqBnd_UpperLimit

  ! Power law source properties
  real(kind=dp) :: pl_index = 1.0            ! Power law index
  real(kind=dp) :: pl_minfreq           ! Minimum frequency for integration of total power
  real(kind=dp) :: pl_maxfreq           ! Maximum frequency for integration of total power
  real(kind=dp) :: pl_scaling = 1.0     ! The scaling of the flux (needs to be initialized)
  real(kind=dp) :: Edd_Efficiency = 0.0 ! Eddinton efficieny
  real(kind=dp) :: pl_S_star = 1.0

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

  ! Frequency dependence of cross section (subband dependent)
  real(kind=dp), dimension(0:NumFreq) :: cross_section_freq_dependence

  ! Weights
  real(kind=dp), dimension(0:NumFreq, 0:NumTau) :: vector_weight

  ! Frequency array
  real(kind=dp), dimension(0:NumFreq) :: frequency

  ! Optical depth array
  real(kind=dp), dimension(0:NumTau) :: tau

  ! Integration table ( optical depth, sub-bin )
  real(kind=dp),dimension(:,:), target, allocatable :: bb_photo_thick_table  
  real(kind=dp),dimension(:,:), target, allocatable :: bb_photo_thin_table
  real(kind=dp),dimension(:,:), target, allocatable :: pl_photo_thick_table 
  real(kind=dp),dimension(:,:), target, allocatable :: pl_photo_thin_table
  real(kind=dp),dimension(:,:), target, allocatable :: bb_heat_thick_table  
  real(kind=dp),dimension(:,:), target, allocatable :: bb_heat_thin_table
  real(kind=dp),dimension(:,:), target, allocatable :: pl_heat_thick_table  
  real(kind=dp),dimension(:,:), target, allocatable :: pl_heat_thin_table

  ! photrates contains all the photo-ionization rates and heating rates
  type photrates    
     real(kind=dp) :: photo_cell_HI          ! HI photoionization rate of the cell    
     real(kind=dp) :: photo_cell_HeI         ! HeI photoionization rate of the cell    
     real(kind=dp) :: photo_cell_HeII        ! HeII photoionization rate of the cell    
     real(kind=dp) :: heat_cell_HI           ! HI heating rate of the cell       
     real(kind=dp) :: heat_cell_HeI          ! HeI heating rate of the cell    
     real(kind=dp) :: heat_cell_HeII         ! HeII heating rate of the cell          
     real(kind=dp) :: photo_in_HI            ! HI photoionization rate incoming to the cell    
     real(kind=dp) :: photo_in_HeI           ! HeI photoionization rate incoming to the cell
     real(kind=dp) :: photo_in_HeII          ! HeII photoionization rate incoming to the cell
     real(kind=dp) :: heat_in_HI             ! HI heating rate incoming to the cell
     real(kind=dp) :: heat_in_HeI            ! HeI heating rate incoming to the cell
     real(kind=dp) :: heat_in_HeII           ! HeII heating rate incoming to the cell 
     real(kind=dp) :: photo_out_HI           ! HI photoionization rate outgoing from the cell
     real(kind=dp) :: photo_out_HeI          ! HeI photoionization rate outgoing from the cell
     real(kind=dp) :: photo_out_HeII         ! HeII photoionization rate outgoing from the cell 
     real(kind=dp) :: heat_out_HI            ! HI heating rate outgoing from the cell
     real(kind=dp) :: heat_out_HeI           ! HeI heating rate outgoing from the cell
     real(kind=dp) :: heat_out_HeII          ! HeII heating rate outgoing from the cell
     real(kind=dp) :: heat                   ! Total heating rate of the cell
     real(kind=dp) :: photo_in               ! Total photoionization rate incoming to the cell
  end type photrates

  ! This definition allows adding two variables of type photrates using the 
  ! + sign.
  ! The function photrates_add is defined below.
  interface operator (+)
     module procedure photrates_add
  end interface operator (+)

  ! tablepos helps to locate correct position of the photoionization and heating tables
  type tablepos
    real(kind=dp), dimension(NumFreqBnd) :: tau            
    real(kind=dp), dimension(NumFreqBnd) :: odpos          
    real(kind=dp), dimension(NumFreqBnd) :: residual       
    integer, dimension(NumFreqBnd)       :: ipos           
    integer, dimension(NumFreqBnd)       :: ipos_p1        
  end type tablepos 

#ifdef MPI       
    integer,private :: mympierror
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                             !
  subroutine rad_ini ()                                                                      !
                                                                                             !
    ! Ask for the parameters of the spectrum                                                 !
    call spectrum_parms ()                                                                   ! 
                                                                                             ! 
    ! Initializes constants and tables for radiation processes                               ! 
    ! (heating, cooling and ionization)                                                      !
    call setup_scalingfactors ()                                                             !
                                                                                             !
    ! Initialize integration routines                                                        !
    call romberg_initialisation (NumFreq)                                                    !
                                                                                             !
    ! Determine spectrum diagnostics                                                         !
    call spec_diag ()                                                                        !

#ifdef MPILOG
    write(logf,*) 'about to integrate'
#endif
                                                                                           !                                                                                             !
    ! Generate photoionization tables and heating tables                                     !
    call spec_integration ()                                                                 !                            

#ifdef MPILOG
    write(logf,*) 'end of radiation'
#endif
                                                                                           !
  end subroutine rad_ini                                                                     !                        
                                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Ask for the parameters of the spectrum
  subroutine spectrum_parms
    
    use file_admin, only: stdinput, file_input
    
    integer :: i_choice = 0                 ! option number
    real(kind=dp) :: bb_luminosity_unscaled ! black body total flux
    
    ! Ask for the input if you are processor 0 and the
    ! spectral parameters are not set in the c2ray_parameters
    ! Note that it is assumed that if teff_nominal is set, 
    ! S_star_nominal is ALSO set.
    if (T_eff_nominal == 0.0) then
       
       if (rank == 0) then
          do while ((sourcetype  /= 'B') .and. (sourcetype /= 'P'))
             
             if (.not.file_input) &
                  write(*,"(A,$)") "Blackbody source (B) or power law source (P)"
             read(stdinput,*) sourcetype ! Read source type, either blackbody or power law source
             
          enddo
          
          !In blackbody case, ask for effective temperature of the source. 
          !The temperature should be bounded below by 2000 and above by 1000000.
          !And then ask for some parameters of the blackbody.
          if (sourcetype == 'B') then
             
             write(logf,*) 'Black body source'
             
             do while (T_eff < 2000.0 .or. T_eff > 1000000.) 
                if (.not.file_input) write(*,'(A,$)') 'Give black body effective temperature (2000 <= T <= 1000000): '
                read(stdinput,*) T_eff      ! Read temperature of black body
                write(logf,*) 'Temperature of the black body : ', T_eff
                if (T_eff < 2000.0 .or. T_eff > 1000000.) then
                   write(*,*) 'Error: Effective temperature out of range. Try again'
                endif
             enddo
             
             ! Find total flux of blackbody (Stefan-Boltzmann law)
             bb_luminosity_unscaled = sigma_SB*T_eff*T_eff*T_eff*T_eff
          
             ! Ask for radius, luminosity, ionizing luminosity or ionizing photon rate?
             if (.not.file_input) then
                write(*,'(A)') 'You can specify' 
                write(*,'(A)') ' 1) a stellar radius'
                write(*,'(A)') ' 2) a total luminosity'
                write(*,'(A)') ' 3) Total ionizing luminosity'
                write(*,'(A)') ' 4) Total number of ionizing photons'
             endif
             
             ! Report error if options are not 1, 2, 3 and 4
             do while (i_choice <= 0 .or. i_choice > 4)
                if (.not.file_input) write(*,'(A,$)') 'Preferred option (1, 2, 3 or 4): '
                read(stdinput,*) i_choice       ! Read option from the input, 1 to 4
                if (i_choice <= 0 .or. i_choice > 4) then
                   write(*,*) 'Error: Choose between 1 2 3 or 4'
                endif
             enddo
             
             select case (i_choice)
                
             case (1)
                write(logf,*) 'A stellar radius is specified'
                if (.not.file_input) write(*,'(A,$)') 'Give radius in solar radius: '
                read(stdinput,*) R_star      ! Read radius of the black body
                write(logf,*) 'The radius is ', R_star, ' solar radius'
                R_star=R_star*r_solar
                L_star=4.0d0*pi*R_star*R_star*bb_luminosity_unscaled
                S_star=0.0  ! Number of photo-ionizing photons is set to zero here but will be reset in spec_diag routine
                
             case (2)
                write(logf,*) 'A total luminosity is specified'
                if (.not.file_input) write(*,'(A,$)') 'Give total luminosity in solar luminosity: '
                read(stdinput,*) L_star      ! Read luminosity of the black body
                write(logf,*) 'The luminosity is ', L_star, ' solar luminosity'
                L_star=L_star*l_solar
                R_star=dsqrt(L_star/(4.0d0*pi*bb_luminosity_unscaled))
                S_star=0.0   ! Number of photo-ionizing photons is set to zero here but will be reset in spec_diag routine
                
             case (3)
                write(logf,*) 'Total ionizing luminosity is specified'
                if (.not.file_input) write(*,'(A,$)') 'Give total ionizing luminosity in solar luminosity: '
                read(stdinput,*) L_star_ion   ! Read ionizing luminosity of the black body
                write(logf,*) 'Total ionizing luminosity is ', L_star_ion, ' solar luminosty'
                L_star_ion=L_star_ion*l_solar         	
                ! Assign some fiducial values, these are overwritten in routine spec_diag
                R_star=r_solar
                L_star=4.0d0*pi*R_star*R_star*bb_luminosity_unscaled
                S_star=0.0   ! Number of photo-ionizing photons is set to zero here but will be reset in spec_diag routine
                
             case (4)
                write(logf,*) 'Rate of ionzing photons (S_star) is specified'
                if (.not.file_input) write(*,'(A,$)') 'Give the number of ionizing photons per second: '
                read(stdinput,*) S_star
                write(logf,*) 'The number of photons per second is ', S_star
                ! Assign some fiducial values for R_star and L_star, 
                ! these are scaled to correspond to S_star in routine spec_diag
                R_star=r_solar
                L_star=4.0d0*pi*R_star*R_star*bb_luminosity_unscaled
                
             end select
             
             ! In power-law case, we ask for number of ionizing photons per second or Eddinton luminosity efficiency 
          elseif(sourcetype == 'P') then
             
             write(logf,*) 'Power law source'
             if (.not.file_input) then
                write(*,'(A)') 'You can specify'
                write(*,'(A)') ' 1)  Number of ionizing photons per second '
                write(*,'(A)') ' 2)  Efficiency parameter assuming a 1e6 solar mass BH' 
             endif
             
             ! Report error if options are not 1 and 2
             do while (i_choice <= 0 .or. i_choice > 2)
                if (.not.file_input) write(*,'(A,$)') 'Preferred option (1 or 2): '
                read(stdinput,*) i_choice
                if (i_choice <= 0 .or. i_choice > 2) then
                   write(*,*) 'Error: Choose between 1 or 2'
                endif
             enddo
             
             select case (i_choice)
                
             case (1)
                write(logf,*) 'Rate of ionizing photons is specified'
                if (.not.file_input) write(*,'(A,$)') 'give number of ionizing photons per second '                        	      
                read(stdinput,*) pl_S_star          ! Read ionizing photons per second
                write(logf,*) 'The rate is ', pl_S_star
                
                ! Set the Eddinton luminosity efficiency to the nominal value 
                Edd_Efficiency=EddLeff_nominal
                pl_scaling=1.0 ! fiducial value, to be updated in spec_diag
                
             case (2)
                write(logf,*) 'Efficiency parameter is specified'
                if (.not.file_input) write(*,'(A,$)') 'give efficiency parameter '
                read(stdinput,*) Edd_Efficiency         ! Read Eddington efficiency
                write(logf,*) 'The efficiency parameter is ', Edd_Efficiency
                ! Set some fiducial value, to be updated in spec_diag
                pl_S_star=0.0
                pl_scaling=1.0
             end select
             
             if (.not.file_input) write(*,'(A,$)') 'Specify power law index (for number of photons, not energy) '
             read(stdinput,*) pl_index      ! Read power law index, this number equal to one plus that of energy 
             write(logf,*) 'Power law index is ', pl_index
             if (.not.file_input) write(*,'(A,$)') 'give lower and upper frequency limits in eV '
             read(stdinput,*) pl_MinFreq,pl_MaxFreq     ! Read lower and upper frequency limits in eV	
             write(logf,*) 'The lower energy limit is ', pl_MinFreq, ' eV'
             write(logf,*) 'The upper energy limit is ', pl_MaxFreq, ' eV'
             if (.not.file_input) write(*,'(A)') 'However, this is not implemented right now '          	

             ! Convert eVs to Hz for the frequency limits
             pl_MinFreq = pl_MinFreq * ev2fr
             pl_MaxFreq = pl_MaxFreq * ev2fr

             ! set some fiducial values for the BB source here, though they are not useful
             R_star=r_solar
             S_star=0.0
             L_star=0.0
             T_eff=1.0e5
             
          endif
       endif
#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(T_eff,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(R_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(L_star,1,MPI_DOUBLE_PRECISION,0, & 
         MPI_COMM_NEW,mympierror)
    call MPI_BCAST(L_star_ion,1,MPI_DOUBLE_PRECISION,0, & 
         MPI_COMM_NEW,mympierror)
    call MPI_BCAST(S_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(Edd_Efficiency,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(pl_S_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(pl_MinFreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
    call MPI_BCAST(pl_MaxFreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         mympierror)
#endif
    
       ! In case neither blackbody nor power-law source
    else
       ! T_eff and S_star are assumed to have been set in the c2ray_parameter module
       T_eff=T_eff_nominal
       S_star=S_star_nominal
       ! Assign some fiducial values, these are scaled to correspond 
       ! to S_star in routine spec_diag
       R_star=r_solar
       L_star=R_star*R_star*4.0d0*pi*sigma_SB*T_eff**4
       
       ! Power law source properties set to nominal values from c2ray_parameter module
       Edd_Efficiency=EddLeff_nominal
       pl_index=pl_index_nominal
       ! Assign some fiducial values, these are scaled to correspond 
       ! to S_star in routine spec_diag
       pl_scaling=1.0
       pl_S_star=0.0
       pl_MinFreq=pl_MinFreq_nominal
       pl_MaxFreq=pl_MaxFreq_nominal
    endif
    
    ! This is h/kT
    h_over_kT = hplanck/(k_B*T_eff)
    
  end subroutine spectrum_parms
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Determine spectrum diagnostics
  ! This routine integrates the SEDs in order to achieve proper scaling
  ! of them. The spectrum_parms subroutine has set some of the SED properties
  ! but they are not necessarily consistent. Consistency is set here
  subroutine spec_diag ()
    
    real(kind=dp) :: bb_ionizing_luminosity_unscaled
    real(kind=dp) :: S_star_unscaled, S_scaling
    real(kind=dp) :: pl_S_star_unscaled
    real(kind=dp) :: pl_S_star_wanted
    real(kind=dp) :: pl_ionizing_luminosity_unscaled
    real(kind=dp) :: pl_ionizing_luminosity
    real(kind=dp) :: pl_ionizing_luminosity_wanted
    real(kind=dp) :: bb_S_star_band1, bb_S_star_band2, bb_S_star_band3
    real(kind=dp) :: pl_S_star_band1, pl_S_star_band2, pl_S_star_band3
    
    ! Case L_star_ion is provided. Find R_star, L_star and blackblody flux
    if (L_star_ion /= 0.0d0) then
       
       ! Black-body ionizing luminosity (energy sense)
       bb_ionizing_luminosity_unscaled = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"B","L")
       ! Find radius from the scaled and specified ionizing luminosities
       R_star = sqrt(L_star_ion/bb_ionizing_luminosity_unscaled)*R_star
       ! Find total luminosity from Stefan-Boltzmann law
       L_star = R_star*R_star*4.0_dp*pi*sigma_SB*T_eff**4
       ! Find the S_star
       S_star = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"B","S")
       
    else
       
       ! Case L_star_ion is not provided.  
       ! Now we know R_star and L_star (except in the case where S_star is provided).
       ! So we can continue to find blackbody flux 
       
       ! Black-body flux (photon sense)
       S_star_unscaled = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"B","S")
       
       ! If S_star is zero, it is set here.
       if (S_star == 0.0) then
          S_star=S_star_unscaled
       else
          ! If S_star is specified, then R_star and L_star have to be tuned accordingly.
          S_scaling = S_star/S_star_unscaled
          R_star = sqrt(S_scaling)*R_star
          L_star = S_scaling*L_star
          L_star_ion = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"B","L")
       endif
    endif
       
    ! Report back to the log file
    if (rank == 0) then
       if(sourcetype == 'B' .or. sourcetype == " ") then
          write(logf,'(/a)')           'Using a black body with'
          write(logf,'(a,1pe10.3,a)')   ' Teff =       ', T_eff, ' K'
          write(logf,'(a,1pe10.3,a)')   ' Radius =     ', R_star/r_solar, ' R_solar'
          write(logf,'(a,1pe10.3,a)')   ' Luminosity = ', L_star/l_solar, ' L_solar'
          write(logf,'(a,1pe10.3,a)')   ' Ionzing photon rate = ', S_star, ' s^-1'
       endif
    endif

    ! Determine the scaling factor pl_scaling
    ! needed to achieve either the specified ionizing photon rate or ionizing luminosity
    if (pl_S_star > 0.0) then
       ! Total power-law ionizing photon rate is specified (photon sense)
       pl_S_star_unscaled = integrate_sed(pl_MinFreq,pl_MaxFreq,"P","S")
       !pl_S_star_unscaled = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"P","S")
       pl_S_star_wanted = pl_S_star
       pl_scaling = pl_S_star_wanted/pl_S_star_unscaled
       pl_ionizing_luminosity = integrate_sed(pl_MinFreq,pl_MaxFreq,"P","L")
       !pl_ionizing_luminosity = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"P","L")
    else
       ! The power-law ionizing luminosity is specified (energy sense). 
       pl_ionizing_luminosity = EddLum*Edd_Efficiency
       pl_ionizing_luminosity_unscaled = integrate_sed(pl_MinFreq,pl_MaxFreq,"P","L")
       !pl_ionizing_luminosity_unscaled = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"P","L")
       pl_ionizing_luminosity_wanted = pl_ionizing_luminosity
       pl_scaling = pl_ionizing_luminosity_wanted/pl_ionizing_luminosity_unscaled
       pl_S_star = integrate_sed(pl_MinFreq,pl_MaxFreq,"P","S")
       !pl_S_star = integrate_sed(freq_min(1),freq_max(NumFreqBnd),"P","S")
    endif
    
    ! Report back to the log file
    if (rank == 0) then
       if (sourcetype == 'P' .or. sourcetype == " ") then
          write(logf,'(/a)')           'Using a power law source with'
          write(logf,'(a,1pe10.3)')   ' Power law index = ', pl_index
          write(logf,'(a,1pe10.3)')   ' Efficiency parameter = ', pl_ionizing_luminosity/EddLum
          write(logf,'(a,1pe10.3)')   ' Ionizing photon rate = ', pl_S_star
          write(logf,'(a,1pe10.3)')   ' Ionizing luminosity = ', pl_ionizing_luminosity
          write(logf,'(a,2f8.3,a)')   ' between energies ', &
               pl_MinFreq/(1e3*ev2fr),pl_MaxFreq/(1e3*ev2fr),' kEv'
       endif
    endif
    
    ! Find out the number of photons in band 1
    bb_S_star_band1 = integrate_sed(freq_min(NumBndin1),freq_max(NumBndin1),"B","S")
    pl_S_star_band1 = integrate_sed(freq_min(NumBndin1),freq_max(NumBndin1),"P","S")
    
    ! Find out the number of photons in band 2
    bb_S_star_band2 = integrate_sed(freq_min(NumBndin1+1),freq_max(NumBndin1+NumBndin2),"B","S")
    pl_S_star_band2 = integrate_sed(freq_min(NumBndin1+1),freq_max(NumBndin1+NumBndin2),"P","S")
    
    ! Find out the number of photons in band 3
    bb_S_star_band3 = integrate_sed(freq_min(NumBndin1+NumBndin2+1), &
         freq_max(NumBndin1+NumBndin2+NumBndin3),"B","S")
    pl_S_star_band3 = integrate_sed(freq_min(NumBndin1+NumBndin2+1), &
         freq_max(NumBndin1+NumBndin2+NumBndin3),"P","S")

    ! Report back to the log file
    if (rank == 0) then
       write(logf,'(/A,(1PE12.5),A//)') ' Number of BB photons in band 1: ', &
            bb_S_star_band1, ' s^-1'
       write(logf,'(A,(1PE12.5),A//)') ' Number of PL photons in band 1: ', &
            pl_S_star_band1, ' s^-1'
       write(logf,'(A,(1PE12.5),A//)') ' Number of BB photons in band 2: ', &
            bb_S_star_band2, ' s^-1'
       write(logf,'(A,(1PE12.5),A//)') ' Number of PL photons in band 2: ', &
            pl_S_star_band2, ' s^-1'
       write(logf,'(A,(1PE12.5),A//)') ' Number of BB photons in band 3: ', &
            bb_S_star_band3, ' s^-1'
       write(logf,'(A,(1PE12.5),A//)') ' Number of PL photons in band 3: ', &
            pl_S_star_band3, ' s^-1'
       write(logf,'(A,(1PE12.5),A//)') &
            ' Total number of ionizing photons (BB): ', &
            bb_S_star_band1+bb_S_star_band2+bb_S_star_band3, 's^-1'            
       write(logf,'(A,(1PE12.5),A//)') &
            ' Total number of ionizing photons (PL): ', &
            pl_S_star_band1+pl_S_star_band2+pl_S_star_band3, 's^-1'            
    endif
    
  end subroutine spec_diag

  function integrate_sed(freq_min, freq_max, sourcetype, sedtype)
    
    ! function type
    real(kind=dp) :: integrate_sed

    ! arguments
    real(kind=dp),intent(in) :: freq_min
    real(kind=dp),intent(in) :: freq_max
    character(len=1),intent(in) :: sourcetype ! P or B
    character(len=1),intent(in) :: sedtype ! L or S

    integer :: i_freq
    real(kind=dp) :: freq_step
    real(kind=dp),dimension(0:NumFreq) :: frequency
    real(kind=dp),dimension(0:NumFreq) :: weight
    real(kind=dp),dimension(0:NumFreq) :: integrand

    ! Set the frequency step for the integration
    freq_step=(freq_max-freq_min)/real(NumFreq)
    
    ! Fill the frequency array and the weight array
    do i_freq=0,NumFreq
       frequency(i_freq)=freq_min+freq_step*real(i_freq)
       weight(i_freq)=freq_step
    enddo

    select case (sourcetype)

    case("B")
       do i_freq=0,NumFreq
          if (frequency(i_freq)*h_over_kT .le. 709.0_dp) then
             ! this blackbody is in number of photon sense
             integrand(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                  (exp(frequency(i_freq)*h_over_kT)-1.0_dp)  
             ! when the argument of the exponential function gets too high
          else
             integrand(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                  (exp((frequency(i_freq)*h_over_kT)/2.0_dp))/ &
                  (exp((frequency(i_freq)*h_over_kT)/2.0_dp))
          endif
          if (sedtype == "L") integrand(i_freq) = hplanck*frequency(i_freq)*integrand(i_freq)
       enddo
       integrate_sed = 4.0*pi*R_star*R_star*scalar_romberg(integrand,weight,NumFreq,NumFreq,0)

    case("P")
       do i_freq=0,NumFreq
          ! this power-law is in number of photon sense
          integrand(i_freq)=frequency(i_freq)**(-pl_index)
          if (sedtype == "L") integrand(i_freq) = hplanck*frequency(i_freq)*integrand(i_freq)
       enddo
       integrate_sed = pl_scaling*scalar_romberg(integrand,weight,NumFreq,NumFreq,0)

    end select

  end function integrate_sed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Generate photoionization tables and heating tables   
  subroutine spec_integration ()

    integer :: i_tau
    integer :: i_subband

    ! Allocate integrands and table arrays
    call allocate_spec_integration_arrays

    ! GM/130729: These should be set somewhere else:
    ! This is h/kT
    h_over_kT=hplanck/(k_B*T_eff)
    ! This is R_star^2
    R_star2=R_star*R_star


    ! fill the optical depth array used to fill the tables 
    ! it is filled in NumTau logarithmic steps 
    ! from minlogtau to maxlogtau
    do i_tau = 1,NumTau
       tau(i_tau) = 10.0**(minlogtau+dlogtau*real(i_tau-1))
    enddo

    ! Position zero corresponds to zero optical depth
    tau(0)=0.0

    ! Warn about grey opacities:
    if (grey .and. rank == 0) write(logf,*) 'WARNING: Using grey opacities'

    ! Find out how far to follow the BB SED (used in lookuptable routines)
    do i_subband=1,NumFreqBnd
       if (freq_min(i_subband)*h_over_kT > 25.) then
          bb_FreqBnd_UpperLimit=i_subband-1
          exit
       endif
    enddo
    if (rank == 0) then
       write(logf,"(A,I3)") "Using BB up to frequency band ", &
            bb_FreqBnd_UpperLimit
       write(logf,"(A,F10.2,A)") "  this is energy ", &
            freq_min(bb_FreqBnd_UpperLimit+1)/ev2fr," eV"
    endif

    ! Find out how far to follow the PL SED (used in lookuptable routines)
    pl_FreqBnd_UpperLimit=NumFreqBnd
    do i_subband=1,NumFreqBnd
       if (freq_min(i_subband) > pl_MaxFreq) then
          pl_FreqBnd_UpperLimit=i_subband-1
          exit
       endif
    enddo
    pl_FreqBnd_LowerLimit=1
    do i_subband=NumFreqBnd,1,-1
       if (freq_min(i_subband) < pl_MinFreq) then
          pl_FreqBnd_LowerLimit=i_subband
          exit
       endif
    enddo
    if (rank == 0) then
       write(logf,"(2(A,I3))") "Using PL from frequency band ", &
            pl_FreqBnd_LowerLimit," to ",pl_FreqBnd_UpperLimit
       write(logf,"(A,F10.2,A)") "  these are energies ", &
            freq_min(pl_FreqBnd_LowerLimit)/ev2fr," and ", &
            (freq_min(pl_FreqBnd_UpperLimit)+ &
            delta_freq(pl_FreqBnd_UpperLimit)*real(NumFreq))/ev2fr," eV"
    endif

#ifdef MPILOG
    write(logf,*) 'Band 1'
#endif

    ! In frequency band 1, fill in integrands and make tables

    ! Go through all the sub-bin in band 1
    do i_subband=1,NumBndin1  
       
       ! Fill frequency array
       call set_frequency_array(i_subband)
       
       ! Set frequency dependence of cross section
       call set_cross_section_freq_dependence(i_subband, &
            pl_index_cross_section_HI(i_subband),grey)
       
       ! Make photo integrands
       call fill_photo_integrands(i_subband)
       
       ! Assign values to the heating integrands
       if (.not.isothermal) then
          call fill_heating_integrands_HI
       endif
       
       ! Set integration weights
       call set_integration_weights(i_subband)
       
       ! Make photo tables
       call make_photo_tables(i_subband)
       
       ! Make heating tables
       if (.not.isothermal) then
          call make_heat_tables_HI(1)
       endif
       
    enddo
    
#ifdef MPILOG
    write(logf,*) 'Band 2'
#endif
    
    ! In frequency band 2, fill in integrands and make tables
    
    ! Go through all the sub-bin in band 2
    do i_subband=NumBndin1+1,NumBndin1+NumBndin2
       
       ! Fill frequency array
       call set_frequency_array(i_subband)
       
       ! Set frequency dependence of cross section
       call set_cross_section_freq_dependence(i_subband, &
            pl_index_cross_section_HeI(i_subband),grey)
       
       ! Make photo integrands
       call fill_photo_integrands(i_subband)
       
       ! Assign values to the heating integrands
       if (.not.isothermal) then
          call fill_heating_integrands_HI
          call fill_heating_integrands_HeI
       endif
       
       ! Set integration weights
       call set_integration_weights(i_subband)
       
       ! Make photo tables
       call make_photo_tables(i_subband)
       
       ! Make heating tables
       if (.not.isothermal) then
          call make_heat_tables_HI(i_subband*2-2)
          call make_heat_tables_HeI(i_subband*2-1)
       endif
       
    enddo
    
#ifdef MPILOG
    write(logf,*) 'Band 3'
#endif
    !
    ! In frequency band 3, fill in integrands and make tables

    ! Go through all the sub-bin in band 3
    do i_subband=NumBndin1+NumBndin2+1,NumBndin1+NumBndin2+NumBndin3
       
#ifdef MPILOG
       write(logf,*) 'Sub band: ',i_subband
#endif
       
       ! Fill frequency array
       call set_frequency_array(i_subband)
       
       ! Set frequency dependence of cross section
       call set_cross_section_freq_dependence(i_subband, &
            pl_index_cross_section_HeII(i_subband),grey)
       
       ! Make photo integrands
       call fill_photo_integrands(i_subband)
       
       ! Assign values to the heating integrands
       if (.not.isothermal) then
          call fill_heating_integrands_HI
          call fill_heating_integrands_HeI
          call fill_heating_integrands_HeII
       endif
       
#ifdef MPILOG
       write(logf,*) 'Integrate photo tables'
#endif
       
       ! Set integration weights
       call set_integration_weights(i_subband)
       
       ! Make photo tables
       call make_photo_tables(i_subband)
       
#ifdef MPILOG
       write(logf,*) 'Integrate heating tables'
#endif
       ! Make heating tables
       if (.not.isothermal) then
          call make_heat_tables_HeII(i_subband*3-NumBndin2-4)
          call make_heat_tables_HeII(i_subband*3-NumBndin2-3)
          call make_heat_tables_HeII(i_subband*3-NumBndin2-2)
       endif
       
    enddo
    
#ifdef MPILOG
    write(logf,*) 'Tables made, cleaning up'
#endif
    
    ! deallocate the useless photo integrand
    call deallocate_integrand_arrays

#ifdef MPILOG
    write(logf,*) 'Clean!'
#endif

    ! report a table value
    write(logf,*) "bb_photo_thick_table: ",sum(bb_photo_thick_table(0,:))
    write(logf,*) "bb_photo_thin_table: ",sum(bb_photo_thin_table(0,:))
    write(logf,*) "bb_heat_thick_table: ",(bb_heat_thick_table(0,1))
    write(logf,*) "bb_heat_thin_table: ",(bb_heat_thin_table(0,1))
    write(logf,*) "Rstar2: ",R_star2
    
  end subroutine spec_integration
  

!---------------------------------------------------------------------------
  
  subroutine allocate_spec_integration_arrays
    
    call allocate_integrand_arrays
    call allocate_table_arrays

  end subroutine allocate_spec_integration_arrays

!---------------------------------------------------------------------------

  subroutine allocate_integrand_arrays

    ! Photoionization integrand as a function of frequency and tau
    allocate(bb_photo_thick_integrand(0:NumFreq, 0:NumTau))    
    allocate(bb_photo_thin_integrand(0:NumFreq, 0:NumTau)) 
    allocate(pl_photo_thick_integrand(0:NumFreq, 0:NumTau))    
    allocate(pl_photo_thin_integrand(0:NumFreq, 0:NumTau)) 

    ! Heating integrand as a function of frequency and tau
    if (.not.isothermal) then
       allocate(bb_heat_thick_integrand_HI(0:NumFreq, 0:NumTau))  
       allocate(bb_heat_thick_integrand_HeI(0:NumFreq, 0:NumTau))  
       allocate(bb_heat_thick_integrand_HeII(0:NumFreq, 0:NumTau))      
       allocate(bb_heat_thin_integrand_HI(0:NumFreq, 0:NumTau))   
       allocate(bb_heat_thin_integrand_HeI(0:NumFreq, 0:NumTau))   
       allocate(bb_heat_thin_integrand_HeII(0:NumFreq, 0:NumTau))   
       allocate(pl_heat_thick_integrand_HI(0:NumFreq, 0:NumTau))   
       allocate(pl_heat_thick_integrand_HeI(0:NumFreq, 0:NumTau))   
       allocate(pl_heat_thick_integrand_HeII(0:NumFreq, 0:NumTau))      
       allocate(pl_heat_thin_integrand_HI(0:NumFreq, 0:NumTau))   
       allocate(pl_heat_thin_integrand_HeI(0:NumFreq, 0:NumTau))   
       allocate(pl_heat_thin_integrand_HeII(0:NumFreq, 0:NumTau))   
    endif

  end subroutine allocate_integrand_arrays

!---------------------------------------------------------------------------

  subroutine deallocate_integrand_arrays
    
    deallocate(bb_photo_thick_integrand)
    deallocate(bb_photo_thin_integrand)
    deallocate(pl_photo_thick_integrand)
    deallocate(pl_photo_thin_integrand)
    
    ! deallocate the useless heating integrand
    if (.not.isothermal) then
       deallocate(bb_heat_thick_integrand_HI)
       deallocate(bb_heat_thick_integrand_HeI)
       deallocate(bb_heat_thick_integrand_HeII)
       deallocate(bb_heat_thin_integrand_HI)
       deallocate(bb_heat_thin_integrand_HeI)
       deallocate(bb_heat_thin_integrand_HeII)
       deallocate(pl_heat_thick_integrand_HI)
       deallocate(pl_heat_thick_integrand_HeI)
       deallocate(pl_heat_thick_integrand_HeII)
       deallocate(pl_heat_thin_integrand_HI)
       deallocate(pl_heat_thin_integrand_HeI)
       deallocate(pl_heat_thin_integrand_HeII)
    endif
    
  end subroutine deallocate_integrand_arrays

!---------------------------------------------------------------------------

  subroutine allocate_table_arrays

    ! Allocate the table arrays

    ! Photoionization table as a function of photo sub-bin and tau
    allocate(bb_photo_thick_table(0:NumTau, 1:NumFreqBnd))
    allocate(bb_photo_thin_table(0:NumTau, 1:NumFreqBnd))
    allocate(pl_photo_thick_table(0:NumTau, 1:NumFreqBnd))
    allocate(pl_photo_thin_table(0:NumTau, 1:NumFreqBnd))
    
    ! Heating table as a function of heating sub-bin and tau
    if (.not.isothermal) then
       allocate(bb_heat_thick_table(0:NumTau, 1:NumheatBin))
       allocate(bb_heat_thin_table(0:NumTau, 1:NumheatBin))
       allocate(pl_heat_thick_table(0:NumTau, 1:NumheatBin))
       allocate(pl_heat_thin_table(0:NumTau, 1:NumheatBin))
    endif
    
  end subroutine allocate_table_arrays

!---------------------------------------------------------------------------

  subroutine set_frequency_array(i_subband)
    
    integer,intent(in) :: i_subband
    
    integer :: i_freq
    
    do i_freq=0,NumFreq
       frequency(i_freq) = freq_min(i_subband)+ &
            delta_freq(i_subband)*real(i_freq)
    enddo
    
  end subroutine set_frequency_array

!---------------------------------------------------------------------------

  subroutine set_cross_section_freq_dependence(i_subband,pl_index,grey)

    integer,intent(in) :: i_subband
    real(kind=dp),intent(in) :: pl_index
    logical,intent(in) :: grey
    
    integer :: i_freq

    if (grey) then
       do i_freq=0,NumFreq
          cross_section_freq_dependence(i_freq) = 1.0
       enddo
    else
       do i_freq=0,NumFreq
          cross_section_freq_dependence(i_freq) = &
               (frequency(i_freq)/freq_min(i_subband))**(-pl_index)        
       enddo
    endif
    
  end subroutine set_cross_section_freq_dependence
    
!---------------------------------------------------------------------------

  subroutine fill_photo_integrands(i_subband)

    integer,intent(in) :: i_subband

    integer :: i_tau
    integer :: i_freq
    
    ! Loop through the tau partition
    do i_tau=0,NumTau 
       
       ! Loop through the frequency partition 
       do i_freq=0,NumFreq
          
          ! Assign values to the photo integrands
          if (tau(i_tau)*cross_section_freq_dependence(i_freq) < 700.0) then
             ! GM/130729 For these high frequencies this
             ! BB exponential term can overflow. Test for this.
             if (frequency(i_freq)*h_over_kT < 700.0) then  
                bb_photo_thick_integrand(i_freq,i_tau) = &
                     4.0_dp*pi*R_star2*two_pi_over_c_square* &
                     frequency(i_freq)*frequency(i_freq)* &
                     exp(-tau(i_tau)*cross_section_freq_dependence(i_freq))/ &
                     (exp(frequency(i_freq)*h_over_kT)-1.0)  
                bb_photo_thin_integrand(i_freq,i_tau) = &
                     4.0_dp*pi*R_star2*two_pi_over_c_square* &
                     frequency(i_freq)*frequency(i_freq)* &
                     cross_section_freq_dependence(i_freq)* &
                     exp(-tau(i_tau)*cross_section_freq_dependence(i_freq))/ &
                     (exp(frequency(i_freq)*h_over_kT)-1.0)  
             else
                bb_photo_thick_integrand(i_freq,i_tau) = 0.0
                bb_photo_thin_integrand(i_freq,i_tau) = 0.0
             endif
             pl_photo_thick_integrand(i_freq,i_tau) = &
                  pl_scaling*frequency(i_freq)**(-pl_index)* &
                  exp(-tau(i_tau)*cross_section_freq_dependence(i_freq))
             pl_photo_thin_integrand(i_freq,i_tau) = &
                  pl_scaling*frequency(i_freq)**(-pl_index)* &
                  cross_section_freq_dependence(i_freq)* &
                  exp(-tau(i_tau)*cross_section_freq_dependence(i_freq))
          else
             bb_photo_thick_integrand(i_freq,i_tau) = 0.0
             bb_photo_thin_integrand(i_freq,i_tau) = 0.0
             pl_photo_thick_integrand(i_freq,i_tau) = 0.0
             pl_photo_thin_integrand(i_freq,i_tau) = 0.0
          endif
          
       enddo
    enddo

  end subroutine fill_photo_integrands

!---------------------------------------------------------------------------

  subroutine fill_heating_integrands_HI

    integer :: i_tau
    integer :: i_freq
    
    ! Loop through the tau partition
    do i_tau=0,NumTau 
       
       ! Loop through the frequency partition
       do i_freq=0,NumFreq
          
          bb_heat_thick_integrand_HI(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HI)* &
               bb_photo_thick_integrand(i_freq,i_tau)
          bb_heat_thin_integrand_HI(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HI)* &
               bb_photo_thin_integrand(i_freq,i_tau)
          pl_heat_thick_integrand_HI(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HI)* &
               pl_photo_thick_integrand(i_freq,i_tau)
          pl_heat_thin_integrand_HI(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HI)* &
               pl_photo_thin_integrand(i_freq,i_tau)
          
       enddo
       
    enddo

  end subroutine fill_heating_integrands_HI

!---------------------------------------------------------------------------

  subroutine fill_heating_integrands_HeI
    
    integer :: i_tau
    integer :: i_freq
    
    ! Loop through the tau partition
    do i_tau=0,NumTau 
       
       ! Loop through the frequency partition
       do i_freq=0,NumFreq
          
          bb_heat_thick_integrand_HeI(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HeI)* &
               bb_photo_thick_integrand(i_freq,i_tau)
          bb_heat_thin_integrand_HeI(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HeI)* &
               bb_photo_thin_integrand(i_freq,i_tau)
          pl_heat_thick_integrand_HeI(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HeI)* &
               pl_photo_thick_integrand(i_freq,i_tau)
          pl_heat_thin_integrand_HeI(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HeI)* &
               pl_photo_thin_integrand(i_freq,i_tau)
       enddo
       
    enddo
    
  end subroutine fill_heating_integrands_HeI

!---------------------------------------------------------------------------

  subroutine fill_heating_integrands_HeII
    
    integer :: i_tau
    integer :: i_freq

    ! Loop through the tau partition
    do i_tau=0,NumTau 
       
       ! Loop through the frequency partition 
       do i_freq=0,NumFreq
          
          bb_heat_thick_integrand_HeII(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HeII)* &
               bb_photo_thick_integrand(i_freq,i_tau)
          bb_heat_thin_integrand_HeII(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HeII)* &
               bb_photo_thin_integrand(i_freq,i_tau)
          pl_heat_thick_integrand_HeII(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HeII)* &
               pl_photo_thick_integrand(i_freq,i_tau)
          pl_heat_thin_integrand_HeII(i_freq,i_tau) = &
               hplanck*(frequency(i_freq)-ion_freq_HeII)* &
               pl_photo_thin_integrand(i_freq,i_tau)
       enddo
       
    enddo

  end subroutine fill_heating_integrands_HeII

!---------------------------------------------------------------------------

  subroutine set_integration_weights(i_subband)

    ! Sets the weights for the call to vector_romberg
    integer,intent(in) :: i_subband

    vector_weight(:,:) = delta_freq(i_subband)  

  end subroutine set_integration_weights

!---------------------------------------------------------------------------

  subroutine make_photo_tables(i_subband)

    integer,intent(in) :: i_subband

    real(kind=dp), dimension(0:NumTau) :: answer

    ! Make photo tables
    call vector_romberg (bb_photo_thick_integrand,vector_weight,NumFreq,NumFreq,NumTau,answer)
    bb_photo_thick_table(:,i_subband) = answer
    call vector_romberg (bb_photo_thin_integrand,vector_weight,NumFreq,NumFreq,NumTau,answer)
    bb_photo_thin_table(:,i_subband) = answer
    call vector_romberg (pl_photo_thick_integrand,vector_weight,NumFreq,NumFreq,NumTau,answer)
    pl_photo_thick_table(:,i_subband) = answer
    call vector_romberg (pl_photo_thin_integrand,vector_weight,NumFreq,NumFreq,NumTau,answer)
    pl_photo_thin_table(:,i_subband) = answer  

  end subroutine make_photo_tables

!---------------------------------------------------------------------------

  subroutine make_heat_tables_HI(table_position)
    
    integer,intent(in) :: table_position

    real(kind=dp), dimension(0:NumTau) :: answer

    call vector_romberg (bb_heat_thick_integrand_HI,vector_weight,NumFreq,NumFreq,NumTau,answer)
    bb_heat_thick_table(:,table_position) = answer
    call vector_romberg (bb_heat_thin_integrand_HI,vector_weight,NumFreq,NumFreq,NumTau,answer)
    bb_heat_thin_table(:,table_position) = answer
    call vector_romberg (pl_heat_thick_integrand_HI,vector_weight,NumFreq,NumFreq,NumTau,answer)
    pl_heat_thick_table(:,table_position) = answer
    call vector_romberg (pl_heat_thin_integrand_HI,vector_weight,NumFreq,NumFreq,NumTau,answer)
    pl_heat_thin_table(:,table_position) = answer
    
  end subroutine make_heat_tables_HI
  
!---------------------------------------------------------------------------

  subroutine make_heat_tables_HeI(table_position)

    integer,intent(in) :: table_position

    real(kind=dp), dimension(0:NumTau) :: answer

    call vector_romberg (bb_heat_thick_integrand_HeI,vector_weight,NumFreq,NumFreq,NumTau,answer)
    bb_heat_thick_table(:,table_position) = answer
    call vector_romberg (bb_heat_thin_integrand_HeI,vector_weight,NumFreq,NumFreq,NumTau,answer)
    bb_heat_thin_table(:,table_position) = answer
    call vector_romberg (pl_heat_thick_integrand_HeI,vector_weight,NumFreq,NumFreq,NumTau,answer)
    pl_heat_thick_table(:,table_position) = answer
    call vector_romberg (pl_heat_thin_integrand_HeI,vector_weight,NumFreq,NumFreq,NumTau,answer)
    pl_heat_thin_table(:,table_position) = answer
    
  end subroutine make_heat_tables_HeI

!---------------------------------------------------------------------------

  subroutine make_heat_tables_HeII(table_position)
    
    integer,intent(in) :: table_position
    
    real(kind=dp), dimension(0:NumTau) :: answer

    call vector_romberg (bb_heat_thick_integrand_HeII,vector_weight,NumFreq,NumFreq,NumTau,answer)
    bb_heat_thick_table(:,table_position) = answer
    call vector_romberg (bb_heat_thin_integrand_HeII,vector_weight,NumFreq,NumFreq,NumTau,answer)
    bb_heat_thin_table(:,table_position) = answer
    call vector_romberg (pl_heat_thick_integrand_HeII,vector_weight,NumFreq,NumFreq,NumTau,answer)
    pl_heat_thick_table(:,table_position) = answer
    call vector_romberg (pl_heat_thin_integrand_HeII,vector_weight,NumFreq,NumFreq,NumTau,answer)
    pl_heat_thin_table(:,table_position) = answer
    
  end subroutine make_heat_tables_HeII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! this subroutine calculates photo-ionization rates at a particular sets of column density
  function photoion_rates (colum_in_HI,colum_out_HI, &
       colum_in_HeI,colum_out_HeI, &
       colum_in_HeII,colum_out_HeII, &
       vol,nsrc,i_state)

    use sourceprops, only: NormFlux,NormFluxPL
    !use cgsphotoconstants

    ! Function type
    type(photrates) :: photoion_rates

    ! Incoming and outgoing HI column density
    real(kind=dp), intent(in) :: colum_in_HI, colum_out_HI

    ! Incoming and outgoing HeI column density
    real(kind=dp), intent(in) :: colum_in_HeI, colum_out_HeI

    ! Incoming and outgoing HeII column density
    real(kind=dp), intent(in) :: colum_in_HeII, colum_out_HeII

    ! Volume of shell cell
    real(kind=dp), intent(in) :: vol

    ! Ionization state of cell
    real(kind=dp), intent(in) :: i_state

    ! Number of the source
    integer, intent(in) :: nsrc 

    integer :: i_subband
    real(kind=dp) :: colum_cell_HI
    real(kind=dp) :: colum_cell_HeI
    real(kind=dp) :: colum_cell_HeII
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_in_all
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_out_all
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_cell_HI
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_cell_HeI
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_cell_HeII
    real(kind=dp), dimension(NumBndin1+1:NumBndin1+NumBndin2+NumBndin3) :: &
         scaling_HI
    real(kind=dp), dimension(NumBndin1+1:NumBndin1+NumBndin2+NumBndin3) :: &
         scaling_HeI
    real(kind=dp), &
         dimension(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+NumBndin3) :: &
         scaling_HeII
    type(tablepos) :: tau_pos_in, tau_pos_out
    type(photrates) :: phi

    ! New source position, set local photo-ionization and heating rates 
    ! to zero. The structure phi is ultimately copied to the result of this function
    call set_photrates_to_zero (phi)

    ! Set the column densities (HI, HeI, HeII) of the current cell
    colum_cell_HI = colum_out_HI-colum_in_HI
    colum_cell_HeI = colum_out_HeI-colum_in_HeI
    colum_cell_HeII = colum_out_HeII-colum_in_HeII 

    ! Calculate the optical depth (incoming, HI, HeI, HeII)
    do i_subband=1,NumFreqBnd
       tau_in_all(i_subband) = colum_in_HI*sigma_HI(i_subband)+ &
            colum_in_HeI*sigma_HeI(i_subband)+ &
            colum_in_HeII*sigma_HeII(i_subband)
    enddo

    ! total tau_out (including HI, HeI, HeII)
    do i_subband=1,NumFreqBnd
       tau_out_all(i_subband) = colum_out_HI*sigma_HI(i_subband)+ &
            colum_out_HeI*sigma_HeI(i_subband)+ &
            colum_out_HeII*sigma_HeII(i_subband)
    enddo

    ! find the table positions for the optical depth (ingoing and outgoing)
    tau_pos_in = set_tau_table_positions(tau_in_all)
    tau_pos_out = set_tau_table_positions(tau_out_all)

    ! Find the scaling factors to be used in distributing rates over
    ! species
    do i_subband=NumBndin1+1,NumBndin1+NumBndin2 ! Band 2
       ! Set the scaling factors to distribute the photo-ionization
       ! over HI and HeI
       call scale_int2(scaling_HI(i_subband),scaling_HeI(i_subband), &
              colum_cell_HI,colum_cell_HeI, i_subband)
    enddo
    do i_subband=NumBndin1+NumBndin2+1,NumBndin1+NumBndin2+NumBndin3 ! Band 3
       ! Set the scaling factors to distribute the photo-ionization
       ! over HI, HeI and HeII
       call scale_int3(scaling_HI(i_subband),scaling_HeI(i_subband), &
            scaling_HeII(i_subband), &
            colum_cell_HI,colum_cell_HeI,colum_cell_HeII,i_subband)
    enddo
    ! Find the photo-ionization rates by looking up the values in
    ! the (appropriate) photo-ionization tables and add to the
    ! rates
    if (NormFlux(nsrc) > 0.0) &  
         phi = phi + photo_lookuptable(tau_pos_in,tau_pos_out, &
         tau_in_all,tau_out_all, &
         NormFlux(nsrc),"B",vol, &
         scaling_HI,scaling_HeI,scaling_HeII)
    !if (colum_in_HI == 0.0) write(logf,*) "After photolookup: ", &
    !     phi%photo_cell_HI, phi%photo_cell_HeI, &
    !              phi%photo_cell_HeII, phi%heat
    if (NormFluxPL(nsrc) > 0.0) &  
         phi = phi + photo_lookuptable(tau_pos_in,tau_pos_out, &
         tau_in_all,tau_out_all, &
         NormFluxPL(nsrc),"P",vol, &
         scaling_HI,scaling_HeI,scaling_HeII)
    
    ! Find the heating rates rates by looking up the values in
    ! the (appropriate) photo-ionization tables and using the
    ! secondary ionization. Add them to the rates.
    if (.not.isothermal) then
       
       ! The optical depths (HI, HeI, HeII) at current cell
       ! These are only needed in heat_lookuptable
       do i_subband=1,NumFreqBnd
          tau_cell_HI(i_subband) = colum_cell_HI*sigma_HI(i_subband)
          tau_cell_HeI(i_subband) = colum_cell_HeI*sigma_HeI(i_subband)
          tau_cell_HeII(i_subband) = colum_cell_HeII*sigma_HeII(i_subband)
       enddo
              
       !if (colum_in_HI == 0.0) then
       !   write(logf,*) "Before heatlookup: ", &
       !        phi%photo_cell_HI, phi%photo_cell_HeI, &
       !        phi%photo_cell_HeII, phi%heat
       !endif
       if (NormFlux(nsrc) > 0.0) & 
            phi = phi + heat_lookuptable(tau_pos_in,tau_pos_out, &
            tau_in_all,tau_out_all, &
            tau_cell_HI,tau_cell_HeI,tau_cell_HeII,NormFlux(nsrc),"B", &
            vol,i_state, &
            scaling_HI,scaling_HeI,scaling_HeII)
       !if (colum_in_HI == 0.0) write(logf,*) "After heatlookup: ", &
       !     phi%photo_cell_HI, phi%photo_cell_HeI, &
       !     phi%photo_cell_HeII, phi%heat

       if (NormFluxPL(nsrc) > 0.0) &  
            phi = phi + heat_lookuptable(tau_pos_in,tau_pos_out, &
            tau_in_all,tau_out_all, &
            tau_cell_HI,tau_cell_HeI,tau_cell_HeII,NormFluxPL(nsrc),"P", &
            vol,i_state, &
            scaling_HI,scaling_HeI,scaling_HeII)
       
    endif

    ! Assign result of function
    photoion_rates = phi

  end function photoion_rates
 
!---------------------------------------------------------------------------

  ! Calculates the table position data for an optical depth tau
  function set_tau_table_positions (tau)

    real(kind=dp), dimension(1:NumFreqBnd),intent(in) :: tau
    type(tablepos) :: set_tau_table_positions
    type(tablepos) :: tau_position

    integer :: i_subband
    
    
    ! fill the table positions structure for the optical depth tau
    do i_subband=1,NumFreqBnd  
       tau_position%tau(i_subband) = log10(max(1.0e-20_dp,tau(i_subband)))
       tau_position%odpos(i_subband) = min(real(NumTau,dp),max(0.0_dp,1.0+ &
            (tau_position%tau(i_subband)-minlogtau)/dlogtau))
       tau_position%ipos(i_subband) = int(tau_position%odpos(i_subband))
       tau_position%residual(i_subband) = tau_position%odpos(i_subband)- &
            real(tau_position%ipos(i_subband),dp)
       tau_position%ipos_p1(i_subband) = min(NumTau, &
            tau_position%ipos(i_subband)+1)
    enddo
    
    ! Set the return value
    set_tau_table_positions=tau_position

  end function set_tau_table_positions

!---------------------------------------------------------------------------

  function read_table(table,tablesize,table_position,i_subband,i_subband2)
    
    integer,intent(in) :: tablesize
    real(kind=dp),dimension(0:NumTau, 1:tablesize),intent(in) :: table
    type(tablepos),intent(in) :: table_position
    integer,intent(in) :: i_subband
    integer,intent(in) :: i_subband2
    
    real(kind=dp) :: read_table
    
    read_table = table(table_position%ipos(i_subband),i_subband2)+ &
         ( table(table_position%ipos_p1(i_subband),i_subband2)- &
         table(table_position%ipos(i_subband),i_subband2) ) * &
         table_position%residual(i_subband)
    
  end function read_table

  !---------------------------------------------------------------------------

  ! find out the correct position in the photo and heating tables
  function photo_lookuptable(tau_pos_in,tau_pos_out, &
       tau_in_all,tau_out_all, &
       NFlux,table_type, &
       vol,scaling_HI,scaling_HeI,scaling_HeII)
    
    !use cgsphotoconstants
    
    ! Function type
    type(photrates) :: photo_lookuptable
    
    ! Optical depth below which we should use the optically thin tables
    real(kind=dp),parameter :: tau_photo_limit = 1.0e-7 

    type(tablepos), intent(in) :: tau_pos_in, tau_pos_out
    real(kind=dp), intent(in) :: NFlux, vol
    real(kind=dp), dimension(NumBndin1+1:NumBndin1+NumBndin2+NumBndin3),intent(in) :: &
         scaling_HI
    real(kind=dp), dimension(NumBndin1+1:NumBndin1+NumBndin2+NumBndin3),intent(in) :: &
         scaling_HeI
    real(kind=dp), &
         dimension(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+NumBndin3),intent(in) :: &
         scaling_HeII
    real(kind=dp), dimension(NumFreqBnd), intent(in) :: tau_in_all, tau_out_all
    character,intent(in) :: table_type
    
    integer ::  i_subband, i
    real(kind=dp) :: phi_photo_in_all, phi_photo_out_all, phi_photo_all
    real(kind=dp), pointer, dimension(:,:) :: photo_thick_table, photo_thin_table
    integer :: Minimum_FreqBnd
    integer :: Maximum_FreqBnd
        
    ! New source. Set all the rates to zero to initialize them.
    call set_photrates_to_zero (photo_lookuptable)

    ! pointers point to the correct tables to use, BB or PL source
    ! Set the maximum frequency band to consider (and limit the
    ! loop over the subbands below)
    if (table_type == "B") then 
       photo_thick_table => bb_photo_thick_table
       photo_thin_table => bb_photo_thin_table
       Minimum_FreqBnd=1
       Maximum_FreqBnd=bb_FreqBnd_UpperLimit
    elseif (table_type == "P") then
       photo_thick_table => pl_photo_thick_table
       photo_thin_table => pl_photo_thin_table
       Minimum_FreqBnd=pl_FreqBnd_LowerLimit
       Maximum_FreqBnd=pl_FreqBnd_UpperLimit
    endif
    
    ! loop through the relevant frequency bands
    do i_subband=Minimum_FreqBnd, Maximum_FreqBnd
       
       ! Incoming total photoionization rate
       phi_photo_in_all = NFlux* &
            read_table(photo_thick_table,NumFreqBnd,tau_pos_in, &
            i_subband,i_subband)
       photo_lookuptable%photo_in = &
            photo_lookuptable%photo_in + phi_photo_in_all

       ! Total cell photo-ionization rate, calculated differently
       ! for optically thick and thin cells
       if (abs(tau_out_all(i_subband)-tau_in_all(i_subband)) > &
            tau_photo_limit) then
          
          ! When current cell is optically thick
          phi_photo_out_all = NFlux* &
               read_table(photo_thick_table,NumFreqBnd, &
               tau_pos_out,i_subband,i_subband)
          phi_photo_all = phi_photo_in_all-phi_photo_out_all
          
       else
          
          ! When current cell is optically thin
          phi_photo_all = NFlux* &
               (tau_out_all(i_subband)-tau_in_all(i_subband))* &
               read_table(photo_thin_table,NumFreqBnd, &
               tau_pos_in,i_subband,i_subband)
          phi_photo_out_all = phi_photo_in_all-phi_photo_all

       endif

      ! Current cell individual photoionization rate of HI, HeI, HeII
      select case (i_subband) 

      ! band 1
      case (NumBndin1) 
     
         ! Assign to the HI photo-ionization rate
         photo_lookuptable%photo_cell_HI = photo_lookuptable%photo_cell_HI + &
              phi_photo_all/vol
         
         ! band 2
      case (NumBndin1+1:NumBndin1+NumBndin2)
                  
         ! Assign to the HI photo-ionization rate
         photo_lookuptable%photo_cell_HI = photo_lookuptable%photo_cell_HI + &
              scaling_HI(i_subband)*phi_photo_all/vol 
         ! Assign to the HeI photo-ionization rate
         photo_lookuptable%photo_cell_HeI = photo_lookuptable%photo_cell_HeI + &
              scaling_HeI(i_subband)*phi_photo_all/vol
         
         ! band 3
      case (NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+NumBndin3)
         
         ! Assign to the HI photo-ionization rate
         photo_lookuptable%photo_cell_HI = photo_lookuptable%photo_cell_HI + &
              scaling_HI(i_subband)*phi_photo_all/vol

         ! Assign the HeI photo-ionization rate
         photo_lookuptable%photo_cell_HeI = photo_lookuptable%photo_cell_HeI +&
              scaling_HeI(i_subband)*phi_photo_all/vol

         ! Assign the HeII photo-ionization rate
         photo_lookuptable%photo_cell_HeII = photo_lookuptable%photo_cell_HeII +&
              scaling_HeII(i_subband)*phi_photo_all/vol

      end select
      
   enddo

 end function photo_lookuptable
 
 !---------------------------------------------------------------------------
 
 ! find out the correct position in the photo and heating tables.
 ! it updates phi
  function heat_lookuptable (tau_pos_in,tau_pos_out, &
       tau_in_all,tau_out_all, &
       tau_cell_HI,tau_cell_HeI,tau_cell_HeII, &
       NFlux,table_type, &
       vol,i_state,scaling_HI,scaling_HeI,scaling_HeII)

    !use cgsphotoconstants

    ! Function type
    type(photrates) :: heat_lookuptable

    ! Optical depth below which we should use the optically thin tables
    real(kind=dp),parameter :: tau_heat_limit = 1.0e-4

    type(tablepos), intent(in) :: tau_pos_in, tau_pos_out
    real(kind=dp), intent(in) :: NFlux, vol, i_state
    real(kind=dp), dimension(NumBndin1+1:NumBndin1+NumBndin2+NumBndin3),intent(in) :: &
         scaling_HI
    real(kind=dp), dimension(NumBndin1+1:NumBndin1+NumBndin2+NumBndin3),intent(in) :: &
         scaling_HeI
    real(kind=dp), &
         dimension(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+NumBndin3),intent(in) :: &
         scaling_HeII
    real(kind=dp), dimension(1:NumFreqBnd), intent(in) :: tau_in_all, &
         tau_out_all
    character,intent(in) :: table_type
    real(kind=dp), dimension(1:NumFreqBnd),intent(in) :: tau_cell_HI, &
         tau_cell_HeI, tau_cell_HeII

    integer ::  i_subband, i
    real(kind=dp) :: phi_heat_HI, phi_heat_HeI, phi_heat_HeII
    real(kind=dp) :: phi_heat_in_HI, phi_heat_in_HeI, phi_heat_in_HeII
    real(kind=dp) :: phi_heat_out_HI, phi_heat_out_HeI, phi_heat_out_HeII
    real(kind=dp) :: f_heat, f_ion_HI, f_ion_HeI
    real(kind=dp) :: fra_sum1, fra_sum2, fra_sum3, fra_sum4
    real(kind=dp), pointer, dimension(:,:) :: heat_thick_table, heat_thin_table
    integer :: Minimum_FreqBnd
    integer :: Maximum_FreqBnd
    ! Related to secondary ionizations
    real(kind=dp), dimension(1:3) :: y1R, y2R
    real(kind=dp) :: xeb

    ! New source. Set all the rates to zero to initialize them.
    call set_photrates_to_zero (heat_lookuptable)

    ! pointers point to the correct tables to use, BB or PL source
    ! Set the maximum frequency band to consider (and limit the
    ! loop over the subbands below)
    if (table_type == "B") then 
       heat_thick_table => bb_heat_thick_table
       heat_thin_table => bb_heat_thin_table
       Minimum_FreqBnd=1
       Maximum_FreqBnd=bb_FreqBnd_UpperLimit
    elseif (table_type == "P") then
       heat_thick_table => pl_heat_thick_table
       heat_thin_table => pl_heat_thin_table
       Minimum_FreqBnd=pl_FreqBnd_LowerLimit
       Maximum_FreqBnd=pl_FreqBnd_UpperLimit
    endif

    ! initialization to zero of local cumulative variables
    ! These variables collect the heating rate (f_heat)
    ! and the additional photo-ionization rates due
    ! to secondary ionization (f_ion_HI and f_ion_HeI).
    ! At the end the values for this source are assigned
    ! to the function result, heat_lookuptable%
    f_heat = 0.0_dp
    f_ion_HI = 0.0_dp
    f_ion_HeI = 0.0_dp
    fra_sum1 = 0.0
    fra_sum2 = 0.0
    fra_sum3 = 0.0
    fra_sum4 = 0.0

    ! Set parameters for secondary ionizations (following Ricotti et al. 2002)
    do i=1,3
       y1R(i)= CR1(i)*(1.0_dp-i_state**bR1(i))**dR1(i)
       xeb=1.0_dp-i_state**bR2(i) 
       y2R(i)= CR2(i)*i_state**aR2(i)*xeb*xeb
    enddo
    !if (tau_in_all(1) == 0.0) write(logf,*) "sec ion parms: ",y1R,xeb,y2R

    ! Current cell individual heating rates of HI, HeI, HeII
    ! loop through the frequency bands
    do i_subband=Minimum_FreqBnd,Maximum_FreqBnd
       
       ! For every subband these will contain the heating due to
       ! the different species by photons in that subband.
       ! The sum over all subbands is collected in f_heat.
       ! All variables starting will phi_ are rates for one subband
       ! and local variables. The heat_lookuptable% is final answer for
       ! all subbands and sources.
       phi_heat_HI = 0.0_dp
       phi_heat_HeI = 0.0_dp
       phi_heat_HeII = 0.0_dp

       select case (i_subband)
          
          ! band 1
       case (NumBndin1)
          
          ! Incoming current cell HI heating rate at band 1
          phi_heat_in_HI = NFlux * read_table(heat_thick_table,NumheatBin, &
               tau_pos_in,i_subband,i_subband)
          
          ! When current cell is HI optically thick     
          if (tau_cell_HI(i_subband) > tau_heat_limit) then
             phi_heat_out_HI = NFlux * read_table(heat_thick_table,NumheatBin, &
                  tau_pos_out,i_subband,i_subband)
             phi_heat_HI = (phi_heat_in_HI-phi_heat_out_HI)/vol
             
             ! When current cell is HI optically thin
          else
             phi_heat_HI = NFlux * tau_cell_HI(i_subband) * &
                  read_table(heat_thin_table,NumheatBin, &
                  tau_pos_in,i_subband,i_subband)
             phi_heat_out_HI = phi_heat_in_HI+phi_heat_HI
             phi_heat_HI = phi_heat_HI/vol
          endif
          
          ! Save the heating in f_heat variable
          f_heat = phi_heat_HI

          ! band 2
       case (NumBndin1+1:NumBndin1+NumBndin2)
          
          ! Band 2, HI
          phi_heat_in_HI = NFlux * read_table(heat_thick_table,NumheatBin, &
               tau_pos_in,i_subband,2*i_subband-2)
          
          ! When current cell is HI optically thick  
          if (tau_cell_HI(i_subband) > tau_heat_limit) then
             phi_heat_out_HI = NFlux * &
                  read_table(heat_thick_table,NumheatBin, &
                  tau_pos_out,i_subband,2*i_subband-2)
             phi_heat_HI = scaling_HI(i_subband)*(phi_heat_in_HI-phi_heat_out_HI)/vol 
             
             ! When current cell is HI optically thin
          else
             phi_heat_HI = NFlux * tau_cell_HI(i_subband) * &
                  read_table(heat_thin_table,NumheatBin, &
                  tau_pos_in,i_subband,2*i_subband-2)
             phi_heat_out_HI = phi_heat_in_HI+phi_heat_HI
             phi_heat_HI = phi_heat_HI/vol
          endif
          
          !if (i_subband == 2) then
          !   write(logf,*) tau_cell_HI(i_subband)
          !   write(logf,*) phi_heat_in_HI,phi_heat_out_HI,phi_heat_HI
          !endif

          ! Band 2, HeI
          phi_heat_in_HeI = NFlux * read_table(heat_thick_table,NumheatBin, &
               tau_pos_in,i_subband,2*i_subband-1)
          
          ! When current cell is HeI optically thick      
          if (tau_cell_HeI(i_subband) > tau_heat_limit) then
             phi_heat_out_HeI = NFlux * &
                  read_table(heat_thick_table,NumheatBin, &
                  tau_pos_out,i_subband,2*i_subband-1)
             phi_heat_HeI = scaling_HeI(i_subband)*(phi_heat_in_HeI-phi_heat_out_HeI)/vol
             
             ! When current cell is HeI optically thin
          else 
             phi_heat_HeI = NFlux * tau_cell_HeI(i_subband) * &
                  read_table(heat_thin_table,NumheatBin, &
                  tau_pos_in,i_subband,2*i_subband-1)
             phi_heat_out_HeI=phi_heat_in_HeI+phi_heat_HeI
             phi_heat_HeI = phi_heat_HeI/vol
          endif
          
          ! Band 2, secondary ionizations
          fra_sum1 = f1ion_HI(i_subband)*phi_heat_HI+ &
               f1ion_HeI(i_subband)*phi_heat_HeI
          fra_sum2 = f2ion_HI(i_subband)*phi_heat_HI+ &
               f2ion_HeI(i_subband)*phi_heat_HeI
          fra_sum3 = f1heat_HI(i_subband)*phi_heat_HI+ &
               f1heat_HeI(i_subband)*phi_heat_HeI
          fra_sum4 = f2heat_HI(i_subband)*phi_heat_HI+ &
               f2heat_HeI(i_subband)*phi_heat_HeI
          
          ! These are all cumulative
          f_ion_HeI = f_ion_HeI+y1R(2)*fra_sum1-y2R(2)*fra_sum2  
          f_ion_HI = f_ion_HI+y1R(1)*fra_sum1-y2R(1)*fra_sum2
          f_heat = f_heat+phi_heat_HI+phi_heat_HeI-y1R(3)*fra_sum3+ &
               y2R(3)*fra_sum4 
          
          ! band 3   
       case (NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+NumBndin3)

          ! Band 3, HI
          phi_heat_in_HI = NFlux * read_table(heat_thick_table,NumheatBin, &
               tau_pos_in,i_subband,3*i_subband-NumBndin2-4)

          ! When current cell is HI optically thick 
          if (tau_cell_HI(i_subband) > tau_heat_limit) then

             phi_heat_out_HI = NFlux * read_table(heat_thick_table,NumheatBin, &
               tau_pos_out,i_subband,3*i_subband-NumBndin2-4)
             phi_heat_HI = scaling_HI(i_subband)* &
                  (phi_heat_in_HI-phi_heat_out_HI)/vol
             
             ! When current cell is HI optically thin
          else
             phi_heat_HI = NFlux * tau_cell_HI(i_subband) * &
                  read_table(heat_thin_table,NumheatBin, &
                  tau_pos_in,i_subband,3*i_subband-NumBndin2-4)
             phi_heat_out_HI = phi_heat_in_HI+phi_heat_HI
             phi_heat_HI = phi_heat_HI/vol
          endif
          
          ! Band 3, HeI
          phi_heat_in_HeI = NFlux * read_table(heat_thick_table,NumheatBin, &
               tau_pos_in,i_subband,3*i_subband-NumBndin2-3)

          ! When current cell is HeI optically thick   
          if (tau_cell_HeI(i_subband) > tau_heat_limit) then 

             phi_heat_out_HeI = NFlux * &
                  read_table(heat_thick_table,NumheatBin, &
                  tau_pos_out,i_subband,3*i_subband-NumBndin2-3)
             phi_heat_HeI = scaling_HeI(i_subband)* &
                  (phi_heat_in_HeI-phi_heat_out_HeI)/vol
             
             ! When current cell is HeI optically thin
          else
             phi_heat_HeI = NFlux * tau_cell_HeI(i_subband) * &
                  read_table(heat_thin_table,NumheatBin, &
                  tau_pos_in,i_subband,3*i_subband-NumBndin2-3)
             phi_heat_out_HeI = phi_heat_in_HeI+phi_heat_HeI
             phi_heat_HeI = phi_heat_HeI/vol
          endif
          
          ! Band 3, HeII
          phi_heat_in_HeII = NFlux * read_table(heat_thick_table,NumheatBin, &
               tau_pos_in,i_subband,3*i_subband-NumBndin2-2)
          
          ! When current cell is HeII optically thick      
          if (tau_cell_HeII(i_subband) > tau_heat_limit) then
             phi_heat_out_HeII = NFlux * &
                  read_table(heat_thick_table,NumheatBin, &
                  tau_pos_out,i_subband,3*i_subband-NumBndin2-2)
             phi_heat_HeII = scaling_HeII(i_subband)* &
                  (phi_heat_in_HeII-phi_heat_out_HeII)/vol
             
             ! When current cell is HeII optically thin
          else 
             phi_heat_HeII = NFlux * tau_cell_HeII(i_subband) * &
                  read_table(heat_thin_table,NumheatBin, &
                  tau_pos_in,i_subband,3*i_subband-NumBndin2-2)
             phi_heat_out_HeII = phi_heat_in_HeII+phi_heat_HeII
             phi_heat_HeII = phi_heat_HeII/vol
          endif
          
          ! Band 3, secondary ionizations
          fra_sum1 = f1ion_HI(i_subband)*phi_heat_HI+ &
               f1ion_HeI(i_subband)*phi_heat_HeI+ &
               f1ion_HeII(i_subband)*phi_heat_HeII
          fra_sum2 = f2ion_HI(i_subband)*phi_heat_HI+ &
               f2ion_HeI(i_subband)*phi_heat_HeI+ &
               f2ion_HeII(i_subband)*phi_heat_HeII
          fra_sum3 = f1heat_HI(i_subband)*phi_heat_HI+ &
               f1heat_HeI(i_subband)*phi_heat_HeI+ &
               f1heat_HeII(i_subband)*phi_heat_HeII
          fra_sum4 = f2heat_HI(i_subband)*phi_heat_HI+ &
               f2heat_HeI(i_subband)*phi_heat_HeI+ &
               f2heat_HeII(i_subband)*phi_heat_HeII
          
          ! These are all cumulative
          f_ion_HeI = f_ion_HeI+y1R(2)*fra_sum1-y2R(2)*fra_sum2
          f_ion_HI = f_ion_HI+y1R(1)*fra_sum1-y2R(1)*fra_sum2
          f_heat = f_heat+phi_heat_HI+phi_heat_HeI+phi_heat_HeII- &
               y1R(3)*fra_sum3+y2R(3)*fra_sum4   
          
       end select
       
    enddo
       
    ! Total heating rate on current cell
    ! Needs to be cumulative because one cell may contain different
    ! types of sources.
    heat_lookuptable%heat = f_heat 
    ! Final HI photoionization rate modified by secondary ionization
    heat_lookuptable%photo_cell_HI = f_ion_HI/(ion_freq_HI*hplanck) 
    ! Final HeI photoionization rate modified by secondary ionization
    heat_lookuptable%photo_cell_HeI = f_ion_HeI/(ion_freq_HeI*hplanck)  
    
  end function heat_lookuptable

  !----------------------------------------------------------------------------

  ! give scalings of species for division of photoionization and 
  ! heating to species for the second frequency bin (He0 ionizing photons).
  ! The scaling is the optical depth of a species over the total optical depth
  ! by all species. This factor is frequency dependent.
  subroutine scale_int2(scaling_HI,scaling_HeI, &
       colum_cell_HI,colum_cell_HeI,i_subband)

    real(kind=dp),intent(in) :: colum_cell_HI, colum_cell_HeI
    integer,intent(in) :: i_subband
    real(kind=dp),intent(out):: scaling_HI, scaling_HeI
    real(kind=dp) :: forscaleing

    forscaleing = 1.0_dp/(sigma_HI(i_subband)*colum_cell_HI+ &
         sigma_HeI(i_subband)*colum_cell_HeI)
    scaling_HI = sigma_HI(i_subband)*colum_cell_HI*forscaleing
    scaling_HeI = sigma_HeI(i_subband)*colum_cell_HeI*forscaleing

  end subroutine scale_int2  

  !----------------------------------------------------------------------------

  ! give scalings of species for division of photoionization and 
  ! heating to species for the third frequency bin (He+ ionizing photons)
  ! The scaling is the optical depth of a species over the total optical depth
  ! by all species. This factor is frequency dependent.
  subroutine scale_int3(scaling_HI, scaling_HeI, scaling_HeII, &
       colum_cell_HI, colum_cell_HeI, colum_cell_HeII, i_subband)

    real(kind=dp),intent(in) :: colum_cell_HI ,colum_cell_HeI, colum_cell_HeII
    integer,intent(in) :: i_subband
    real(kind=dp),intent(out):: scaling_HI, scaling_HeI, scaling_HeII
    real(kind=dp) :: forscaleing

    forscaleing = 1.0_dp/(sigma_HI(i_subband)*colum_cell_HI+ &
         sigma_HeI(i_subband)*colum_cell_HeI+ &
         sigma_HeII(i_subband)*colum_cell_HeII) 
    scaling_HI = colum_cell_HI*sigma_HI(i_subband) *forscaleing
    scaling_HeI = colum_cell_HeI*sigma_HeI(i_subband)*forscaleing
    scaling_HeII = colum_cell_HeII*sigma_HeII(i_subband)*forscaleing 

  end subroutine scale_int3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! set up some scaling factors arrays
  subroutine setup_scalingfactors

    integer :: i_subband

    ! Allocate size of arrays
    allocate(delta_freq(1:NumFreqBnd))  
    allocate(freq_max(1:NumFreqBnd))  
    allocate(freq_min(1:NumFreqBnd))  
    allocate(pl_index_cross_section_HI(1:NumFreqBnd))
    allocate(pl_index_cross_section_HeI(1:NumFreqBnd))
    allocate(pl_index_cross_section_HeII(1:NumFreqBnd))
    allocate(sigma_HI(1:NumFreqBnd))
    allocate(sigma_HeI(1:NumFreqBnd))
    allocate(sigma_HeII(1:NumFreqBnd))

    ! Allocate size of arrays of heating parameters
    if (.not.isothermal) then
      allocate(f1ion_HI(NumBndin1+1:NumFreqBnd))  
      allocate(f1ion_HeI(NumBndin1+1:NumFreqBnd))  
      allocate(f1ion_HeII(NumBndin1+1:NumFreqBnd))  
      allocate(f2ion_HI(NumBndin1+1:NumFreqBnd)) 
      allocate(f2ion_HeI(NumBndin1+1:NumFreqBnd)) 
      allocate(f2ion_HeII(NumBndin1+1:NumFreqBnd)) 
      allocate(f2heat_HI(NumBndin1+1:NumFreqBnd)) 
      allocate(f2heat_HeI(NumBndin1+1:NumFreqBnd)) 
      allocate(f2heat_HeII(NumBndin1+1:NumFreqBnd)) 
      allocate(f1heat_HI(NumBndin1+1:NumFreqBnd)) 
      allocate(f1heat_HeI(NumBndin1+1:NumFreqBnd))
      allocate(f1heat_HeII(NumBndin1+1:NumFreqBnd))
    endif

    ! Assignment of maximum frequency in the sub-bin partition.
    select case (NumBndin1)

      case (1)

        freq_max(NumBndin1) = ion_freq_HeI

    end select

    select case (NumBndin2)

    case (26)
    
      freq_max(NumBndin1+1:26) = ion_freq_HeI* &
                                 (/1.02_dp, 1.05_dp, 1.07_dp, 1.10_dp, 1.15_dp, 1.20_dp, &
		  	           1.25_dp, 1.30_dp, 1.35_dp, 1.40_dp, 1.45_dp, 1.50_dp, &
		                   1.55_dp, 1.60_dp, 1.65_dp, 1.70_dp, 1.75_dp, 1.80_dp, &
		                   1.85_dp, 1.90_dp, 1.95_dp, 2.00_dp, 2.05_dp, 2.10_dp, &
                                   2.15_dp/)
      freq_max(NumBndin1+NumBndin2) = ion_freq_HeII

    case (10)
	
      freq_max(NumBndin1+1:10) = ion_freq_HeI* &
                                 (/1.10_dp, 1.20_dp, 1.30_dp, 1.40_dp, 1.50_dp, &
                                   1.60_dp, 1.70_dp, 1.80_dp, 1.90_dp/)
      freq_max(NumBndin1+NumBndin2) = ion_freq_HeII

    case (6)

      freq_max(NumBndin1+1:6) = ion_freq_HeI* &
                                (/1.15_dp, 1.30_dp, 1.50_dp, 1.70_dp, 1.9557_dp/)
      freq_max(NumBndin1+NumBndin2) = ion_freq_HeII

    case (3)

      freq_max(NumBndin1+1:3) = ion_freq_HeI* &                           
                                (/1.3_dp,1.7_dp/)
      freq_max(NumBndin1+NumBndin2) = ion_freq_HeII

    case (2)   
  
      freq_max(NumBndin1+1) = ion_freq_HeI*1.5_dp
      freq_max(NumBndin1+NumBndin2) = ion_freq_HeII

    case (1)  
      freq_max(NumBndin1+NumBndin2) = ion_freq_HeII

    end select

    select case (NumBndin3)

      case (20)
        freq_max(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = ion_freq_HeII* &
                                (/1.05_dp, 1.10_dp, 1.20_dp, 1.40_dp, 1.70_dp, 2.00_dp, &
			          2.50_dp, 3.00_dp, 4.00_dp, 5.00_dp, 7.00_dp, 10.00_dp, &
			          15.00_dp, 20.00_dp, 30.00_dp, 40.00_dp, 50.00_dp, 70.00_dp, &
                                  90.00_dp, 100.00_dp/)  

      case (16)
        freq_max(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+16) = ion_freq_HeII* &
			        (/1.05_dp, 1.10_dp, 1.20_dp, 1.40_dp, 1.70_dp, 2.00_dp, &
                                  3.00_dp, 5.00_dp, 7.00_dp, 10.00_dp, 15.00_dp, 20.00_dp, &
                                  30.00_dp, 50.00_dp, 70.00_dp, 100.00_dp/)

      case (11)
        freq_max(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+11) = ion_freq_HeII* &
                                (/1.10_dp, 1.20_dp, 1.50_dp, 2.00_dp, 3.00_dp, 4.00_dp, &
                                  7.00_dp, 10.00_dp, 20.00_dp, 50.00_dp, 100.0_dp/)

      case (9)
        freq_max(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+9) = ion_freq_HeII* & 
	                        (/1.50_dp, 2.00_dp, 3.00_dp, 4.00_dp, 7.00_dp, 10.00_dp, &
                                  20.00_dp, 50.00_dp, 100.00_dp/)

      case (4)
        freq_max(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+4) = ion_freq_HeII* &
                                (/2.00_dp, 4.00_dp, 10.00_dp, 100.0_dp/)

      case (1)
        freq_max(NumBndin1+NumBndin2+1) = ion_freq_HeII * 100.00_dp

    end select

    ! Assignment of minimum frequency in the sub-bin partition.
       
    freq_min(NumBndin1) = ion_freq_HI
    do i_subband=2,NumFreqBnd
      freq_min(i_subband) = freq_max(i_subband-1)
    enddo

    ! calculate the width of frequency sub-bin
    do i_subband=1,NumFreqBnd 
       delta_freq(i_subband) = (freq_max(i_subband)-freq_min(i_subband))/real(NumFreq)
    enddo

    ! Assign f_ion and f_heat for secondary ionization
    if (.not.isothermal) then

      select case (NumBndin2)

        case (26)

          f1ion_HI(NumBndin1+1:NumBndin1+26) = (/ 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                  1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                  1.0000_dp, 1.0000_dp/) 

          f1ion_HeI(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 1.0000_dp/) 

          f1ion_HeII(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp/) 

          f2ion_HI(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                 0.9971_dp, 0.9802_dp, 0.9643_dp, 0.9493_dp, &
                                                 0.9350_dp, 0.9215_dp, 0.9086_dp, 0.8964_dp, &
                                                 0.8847_dp, 0.8735_dp/)

          f2ion_HeI(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                  0.0000_dp, 0.9960_dp/) 

          f2ion_HeII(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp/) 

          f1heat_HI(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                    1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                    1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                    1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                    1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                    1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                    1.0000_dp, 1.0000_dp/)

          f1heat_HeI(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 1.0000_dp, &
                                                   1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                   1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                   1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                                                   1.0000_dp, 1.0000_dp /)

          f1heat_HeII(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp /)

          f2heat_HI(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.9704_dp, 0.9290_dp, 0.9037_dp, &
                                                  0.8687_dp, 0.8171_dp, 0.7724_dp, 0.7332_dp, &
                                                  0.6985_dp, 0.6675_dp, 0.6397_dp, 0.6145_dp, &
                                                  0.5916_dp, 0.5707_dp, 0.5514_dp, 0.5337_dp, &
                                                  0.5173_dp, 0.5021_dp, 0.4879_dp, 0.4747_dp, &
                                                  0.4623_dp, 0.4506_dp, 0.4397_dp, 0.4293_dp, &
                                                  0.4196_dp, 0.4103_dp/) 

          f2heat_HeI(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                   0.0000_dp, 0.0000_dp, 0.0000_dp, 0.9959_dp, &
                                                   0.9250_dp, 0.8653_dp, 0.8142_dp, 0.7698_dp, &
                                                   0.7309_dp, 0.6965_dp, 0.6657_dp, 0.6380_dp, &
                                                   0.6130_dp, 0.5903_dp, 0.5694_dp, 0.5503_dp, &
                                                   0.5327_dp, 0.5164_dp/)

          f2heat_HeII(NumBndin1+1:NumBndin1+26) = (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, &
                                                    0.0000_dp, 0.0000_dp/)  

      end select

      select case (NumBndin3)

        case(20)

          f1ion_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp /)

           f1ion_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp /)

           f1ion_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp /)

           f2ion_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) =&
                        (/0.8600_dp, 0.8381_dp, 0.8180_dp, 0.7824_dp, 0.7249_dp, 0.6607_dp, &
                          0.6128_dp, 0.5542_dp, 0.5115_dp, 0.4518_dp, 0.4110_dp, 0.3571_dp, &
                          0.3083_dp, 0.2612_dp, 0.2325_dp, 0.1973_dp, 0.1757_dp, 0.1606_dp, &
                          0.1403_dp, 0.1269_dp /)

           f2ion_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/0.9750_dp, 0.9415_dp, 0.9118_dp, 0.8609_dp, 0.7831_dp, 0.7015_dp, &
                          0.6436_dp, 0.5755_dp, 0.5273_dp, 0.4619_dp, 0.4182_dp, 0.3615_dp, &
                          0.3109_dp, 0.2627_dp, 0.2334_dp, 0.1979_dp, 0.1761_dp, 0.1609_dp, &
                          0.1405_dp, 0.1270_dp /)

           f2ion_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.8841_dp, &
                          0.7666_dp, 0.6518_dp, 0.5810_dp, 0.4940_dp, 0.4403_dp, 0.3744_dp, &
                          0.3183_dp, 0.2668_dp, 0.2361_dp, 0.1993_dp, 0.1771_dp, 0.1616_dp, &
                          0.1409_dp, 0.1273_dp /)

           f1heat_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp /)

           f1heat_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp /)

           f1heat_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, 1.0000_dp, &
                          1.0000_dp, 1.0000_dp /)

           f2heat_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/0.3994_dp, 0.3817_dp, 0.3659_dp, 0.3385_dp, 0.2961_dp, 0.2517_dp, &
                          0.2207_dp, 0.1851_dp, 0.1608_dp, 0.1295_dp, 0.1097_dp, 0.0858_dp, &
                          0.0663_dp, 0.0496_dp, 0.0405_dp, 0.0304_dp, 0.0248_dp, 0.0212_dp, &
                          0.0167_dp, 0.0140_dp /)

           f2heat_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/0.4974_dp, 0.4679_dp, 0.4424_dp, 0.4001_dp, 0.3389_dp, 0.2796_dp, &
                          0.2405_dp, 0.1977_dp, 0.1697_dp, 0.1346_dp, 0.1131_dp, 0.0876_dp, &
                          0.0673_dp, 0.0501_dp, 0.0408_dp, 0.0305_dp, 0.0249_dp, 0.0213_dp, &
                          0.0168_dp, 0.0140_dp /)

           f2heat_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                        (/0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.6202_dp, 0.4192_dp, &
                          0.3265_dp, 0.2459_dp, 0.2010_dp, 0.1513_dp, 0.1237_dp, 0.0932_dp, &
                          0.0701_dp, 0.0515_dp, 0.0416_dp, 0.0309_dp, 0.0251_dp, 0.0214_dp, &
                          0.0169_dp, 0.0141_dp /)

      end select 

    endif

    ! Assign value sigma of HI, HeI and HeII at different frequencies 
    select case (NumBndin1)

      case (1)

        sigma_HI(1) = sigma_HI_at_ion_freq
        sigma_HeI(1) = 0.0_dp
        sigma_HeII(1) = 0.0_dp  

      end select

    select case (NumBndin2)

      case (26) 

        sigma_HI(NumBndin1+1:NumBndin1+26) = (/1.239152e-18_dp, 1.171908e-18_dp, &
             1.079235e-18_dp, 1.023159e-18_dp, 9.455687e-19_dp, 8.329840e-19_dp, &
             7.374876e-19_dp, 6.559608e-19_dp, 5.859440e-19_dp, 5.254793e-19_dp, &
             4.729953e-19_dp, 4.272207e-19_dp, 3.874251e-19_dp, 3.521112e-19_dp, &
             3.209244e-19_dp, 2.932810e-19_dp, 2.686933e-19_dp, 2.467523e-19_dp, &
             2.271125e-19_dp, 2.094813e-19_dp, 1.936094e-19_dp, 1.792838e-19_dp, &
             1.663215e-19_dp, 1.545649e-19_dp, 1.438778e-19_dp, 1.341418e-19_dp/)
        sigma_HeI(NumBndin1+1:NumBndin1+26) = (/7.434699e-18_dp, 7.210641e-18_dp, &
              6.887151e-18_dp, 6.682491e-18_dp, 6.387263e-18_dp, 5.931487e-18_dp, &
              5.516179e-18_dp, 5.137743e-18_dp, 4.792724e-18_dp, 4.477877e-18_dp, &
              4.190200e-18_dp, 3.926951e-18_dp, 3.687526e-18_dp, 3.465785e-18_dp, &
              3.261781e-18_dp, 3.073737e-18_dp, 2.900074e-18_dp, 2.739394e-18_dp, &
              2.590455e-18_dp, 2.452158e-18_dp, 2.323526e-18_dp, 2.203694e-18_dp, &
              2.091889e-18_dp, 1.987425e-18_dp, 1.889687e-18_dp, 1.798126e-18_dp/)
        sigma_HeII(NumBndin1+1:NumBndin1+26) = 0.0_dp 

      case (10) 

        sigma_HI(NumBndin1+1:NumBndin1+10) = (/1.239152e-18_dp, 9.455687e-19_dp, &
             7.374876e-19_dp, 5.859440e-19_dp, 4.729953e-19_dp, 3.874251e-19_dp, &
             3.209244e-19_dp, 2.686933e-19_dp, 2.271125e-19_dp, 1.936094e-19_dp/)
        sigma_HeI(NumBndin1+1:NumBndin1+10) = (/7.434699e-18_dp, 6.387263e-18_dp, &
              5.516179e-18_dp, 4.792724e-18_dp, 4.190200e-18_dp, 3.687526e-18_dp, &
              3.261781e-18_dp, 2.900074e-18_dp, 2.590455e-18_dp, 2.323526e-18_dp/)
        sigma_HeII(NumBndin1+1:NumBndin1+10) = 0.0_dp 

      case (6)
        sigma_HI(NumBndin1+1:NumBndin1+6) = (/1.164e-18_dp, 8.33e-19_dp,5.859e-19_dp, &
                                              3.874e-19_dp,2.687e-19_dp,1.777e-19_dp/)
        sigma_HeI(NumBndin1+1:NumBndin1+6) = (/sigma_HeI_at_ion_freq, 5.9315e-18_dp, &
                                                       4.7927e-18_dp, 3.6875e-18_dp, &
                                                       2.9001e-18_dp, 2.1906e-18_dp/)
        sigma_HeII(NumBndin1+1:NumBndin1+6) = 0.0_dp 

      case (3) 
        sigma_HI(NumBndin1+1:NumBndin1+3) = (/1.239e-18_dp, 5.86e-19_dp, 2.69e-19_dp/)
        sigma_HeI(NumBndin1+1:NumBndin1+3) = (/sigma_HeI_at_ion_freq, 4.793e-18_dp,2.90e-18_dp/)  
        sigma_HeII(NumBndin1+1:NumBndin1+3) = 0.0_dp 

      case (2)
        sigma_HI(NumBndin1+1:NumBndin1+2) = (/1.239e-18_dp, 3.87e-19_dp/)
        sigma_HeI(NumBndin1+1:NumBndin1+2) = (/sigma_HeI_at_ion_freq, 3.688e-18_dp/)  
        sigma_HeII(NumBndin1+1:NumBndin1+2) = 0.0_dp 

      case (1) 
        sigma_HI(NumBndin1+1) = 1.239e-18_dp
        sigma_HeI(NumBndin1+1) = sigma_HeI_at_ion_freq
        sigma_HeII(NumBndin1+1) = 0.0_dp 

    end select

    select case (NumBndin3)

      case (20)

        sigma_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                 (/1.230696e-19_dp, 1.063780e-19_dp, 9.253883e-20_dp, &
                   7.123014e-20_dp, 4.464019e-20_dp, 2.465533e-20_dp, &
                   1.492667e-20_dp, 7.446712e-21_dp, 4.196728e-21_dp, &
                   1.682670e-21_dp, 8.223247e-22_dp, 2.763830e-22_dp, &
                   8.591126e-23_dp, 2.244684e-23_dp, 8.593853e-24_dp, &
                   2.199718e-24_dp, 8.315674e-25_dp, 3.898672e-25_dp, &
                   1.238718e-25_dp, 5.244957e-26_dp/)
        sigma_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                 (/1.690781e-18_dp, 1.521636e-18_dp, 1.373651e-18_dp, &
                   1.128867e-18_dp, 7.845096e-19_dp, 4.825331e-19_dp, &
                   3.142134e-19_dp, 1.696228e-19_dp, 1.005051e-19_dp, &
                   4.278712e-20_dp, 2.165403e-20_dp, 7.574790e-21_dp, &
                   2.429426e-21_dp, 6.519748e-22_dp, 2.534069e-22_dp, &
                   6.599821e-23_dp, 2.520412e-23_dp, 1.189810e-23_dp, &
                   3.814490e-24_dp, 1.624492e-24_dp/)
        sigma_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                 (/1.587280e-18_dp, 1.391911e-18_dp, 1.227391e-18_dp, &
                   9.686899e-19_dp, 6.338284e-19_dp, 3.687895e-19_dp, &
                   2.328072e-19_dp, 1.226873e-19_dp, 7.214988e-20_dp, &
                   3.081577e-20_dp, 1.576429e-20_dp, 5.646276e-21_dp, &
                   1.864734e-21_dp, 5.177347e-22_dp, 2.059271e-22_dp, &
                   5.526508e-23_dp, 2.151467e-23_dp, 1.029637e-23_dp, &
                   3.363164e-24_dp, 1.450239e-24_dp/)

      case (16) 
	
        sigma_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+16) = &
                (/ 1.230696e-19_dp, 1.063780e-19_dp, 9.253883e-20_dp, &
                   7.123014e-20_dp, 4.464019e-20_dp, 2.465533e-20_dp, &	
                   1.492667e-20_dp, 4.196728e-21_dp, 8.223247e-22_dp, &
                   2.763830e-22_dp, 8.591126e-23_dp, 2.244684e-23_dp, &	
                   8.593853e-24_dp, 2.199718e-24_dp, 3.898672e-25_dp, &
                   1.238718e-25_dp/)

        sigma_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+16) = &
                 (/1.690781e-18_dp, 1.521636e-18_dp, 1.373651e-18_dp, &
                   1.128867e-18_dp, 7.845096e-19_dp, 4.825331e-19_dp, &
                   3.142134e-19_dp, 1.005051e-19_dp, 2.165403e-20_dp, &
                   7.574790e-21_dp, 2.429426e-21_dp, 6.519748e-22_dp, &
                   2.534069e-22_dp, 6.599821e-23_dp, 1.189810e-23_dp, &
                   3.814490e-24_dp/)

        sigma_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+16) = &
                 (/1.587280e-18_dp, 1.391911e-18_dp, 1.227391e-18_dp, &
                   9.686899e-19_dp, 6.338284e-19_dp, 3.687895e-19_dp, &
                   2.328072e-19_dp, 7.214988e-20_dp, 1.576429e-20_dp, &
                   5.646276e-21_dp, 1.864734e-21_dp, 5.177347e-22_dp, &
                   2.059271e-22_dp, 5.526508e-23_dp, 1.029637e-23_dp, &
                   3.363164e-24_dp/)

      case (11) 

        sigma_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+11) = &
                       (/1.2307e-19_dp, 9.2539e-20_dp, 7.1230e-20_dp, &
                         3.6176e-20_dp, 1.4927e-20_dp, 4.1967e-21_dp, &
                         1.6827e-21_dp, 2.7638e-22_dp, 8.5911e-23_dp, &
                         8.5939e-24_dp,3.8987e-25_dp/)
        sigma_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+11) = &
                       (/1.6908e-18_dp, 1.3737e-18_dp, 1.1289e-18_dp, &
                         6.6238e-19_dp, 3.1421e-19_dp, 1.0051e-19_dp, &
                         4.2787e-20_dp, 7.5748e-21_dp, 2.4294e-21_dp, &
                         2.5341e-22_dp, 1.1898e-23_dp/)
        sigma_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+11) = &
                       (/1.5873e-18_dp, 1.2274e-18_dp, 9.6869e-19_dp, &
                         5.2339e-19_dp, 2.3281e-19_dp, 7.2150e-20_dp, &
                         3.0816e-20_dp, 5.6463e-21_dp, 1.8647e-21_dp, &
                         2.0593e-22_dp, 1.0296e-23_dp/)

      case (9) 

        sigma_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+9) = &
                 (/1.230696e-19_dp, 3.617600e-20_dp, 1.492667e-20_dp, &
                   4.196728e-21_dp, 1.682670e-21_dp, 2.763830e-22_dp, &
                   8.591126e-23_dp, 8.593853e-24_dp, 3.898672e-25_dp/) 
        sigma_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+9) = &
                 (/1.690781e-18_dp, 6.623773e-19_dp, 3.142134e-19_dp, &
                   1.005051e-19_dp, 4.278712e-20_dp, 7.574790e-21_dp, &
                   2.429426e-21_dp, 2.534069e-22_dp, 1.189810e-23_dp/)
        sigma_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+9) = &
           (/sigma_HeII_at_ion_freq,5.233870e-19_dp, 2.328072e-19_dp, &
                   7.214988e-20_dp, 3.081577e-20_dp, 5.646276e-21_dp, &
                   1.864734e-21_dp, 2.059271e-22_dp, 1.029637e-23_dp/)

      case (4)

        sigma_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+4) = &
          (/1.2307e-19_dp, 1.4927e-20_dp, 1.6827e-21_dp, 8.5900e-23_dp/)
        sigma_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+4) = &
          (/1.6908e-18_dp, 3.1421e-19_dp, 4.2787e-20_dp, 2.4294e-21_dp/)
        sigma_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+4) = &
          (/1.5873e-18_dp, 2.3280e-19_dp, 3.0816e-20_dp, 1.1865e-21_dp/)

      case (1)

        sigma_HI(NumBndin1+NumBndin2+1) = 1.2300e-19_dp
        sigma_HeI(NumBndin1+NumBndin2+1) = 1.691e-18_dp
        sigma_HeII(NumBndin1+NumBndin2+1) = sigma_HeII_at_ion_freq

    end select


    ! Assign power-law index of HI, HeI and HeII at different frequencies (about absorption)
    select case (NumBndin1)

      case (1)

        pl_index_cross_section_HI(1) = 2.761_dp

    end select

    select case (NumBndin2)

      case (26) 

        pl_index_cross_section_HI(NumBndin1+1:NumBndin1+26) = (/2.8277_dp, 2.8330_dp, 2.8382_dp, &
                            2.8432_dp, 2.8509_dp, 2.8601_dp, 2.8688_dp, 2.8771_dp, &
                            2.8850_dp, 2.8925_dp, 2.8997_dp, 2.9066_dp, 2.9132_dp, &
                            2.9196_dp, 2.9257_dp, 2.9316_dp, 2.9373_dp, 2.9428_dp, &
                            2.9481_dp, 2.9532_dp, 2.9582_dp, 2.9630_dp, 2.9677_dp, &
                            2.9722_dp, 2.9766_dp, 2.9813_dp/)
        pl_index_cross_section_HeI(NumBndin1+1:NumBndin1+26) = (/1.5509_dp, 1.5785_dp, 1.6047_dp, &
                             1.6290_dp, 1.6649_dp, 1.7051_dp, 1.7405_dp, 1.7719_dp, &
                             1.8000_dp, 1.8253_dp, 1.8486_dp, 1.8701_dp, 1.8904_dp, &
                             1.9098_dp, 1.9287_dp, 1.9472_dp, 1.9654_dp, 1.9835_dp, &
                             2.0016_dp, 2.0196_dp, 2.0376_dp, 2.0557_dp, 2.0738_dp, &
                             2.0919_dp, 2.1099_dp, 2.1302_dp/)

      case (10) 

        pl_index_cross_section_HI(NumBndin1+1:NumBndin1+10) = (/2.8360_dp, 2.8554_dp, 2.8729_dp, 2.8887_dp, &
                                2.9031_dp,2.9164_dp, 2.9287_dp,2.9400_dp,2.9507_dp,2.9701_dp/)
        pl_index_cross_section_HeI(NumBndin1+1:NumBndin1+10) = (/1.5932_dp, 1.6849_dp, 1.7561_dp, 1.8126_dp, &
                             1.8592_dp, 1.9000_dp, 1.9379_dp, 1.9744_dp, 2.0105_dp, 2.0840_dp/)

      case (6)

        pl_index_cross_section_HI(NumBndin1+1:NumBndin1+6) = (/2.8408_dp, 2.8685_dp, 2.8958_dp, &
                                                 2.9224_dp, 2.9481_dp, 2.9727_dp/)
        pl_index_cross_section_HeI(NumBndin1+1:NumBndin1+6) = (/1.6168_dp, 1.7390_dp, 1.8355_dp, &
                                                  1.9186_dp, 2.0018_dp, 2.0945_dp/)

      case (3) 

        pl_index_cross_section_HI(NumBndin1+1:NumBndin1+3) = (/2.8542_dp, 2.9086_dp, 2.9600_dp/)
        pl_index_cross_section_HeI(NumBndin1+1:NumBndin1+3) = (/1.6770_dp, 1.8758_dp, 2.0458_dp/)

      case (2)

        pl_index_cross_section_HI(NumBndin1+1:NumBndin1+2) = (/2.8697_dp, 2.9486_dp/)
        pl_index_cross_section_HeI(NumBndin1+1:NumBndin1+2) = (/1.7385_dp, 2.0061_dp/)

      case (1) 

        pl_index_cross_section_HI(NumBndin1+1) = 2.9118_dp
        pl_index_cross_section_HeI(NumBndin1+1) = 1.8832_dp

    end select

    select case (NumBndin3)

      case (20)

        pl_index_cross_section_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                            (/2.9884_dp, 2.9970_dp, 3.0088_dp, 3.0298_dp, 3.0589_dp, &
                              3.0872_dp, 3.1166_dp, 3.1455_dp, 3.1773_dp, 3.2089_dp, &
                              3.2410_dp, 3.2765_dp, 3.3107_dp, 3.3376_dp, 3.3613_dp, &
                              3.3816_dp, 3.3948_dp, 3.4078_dp, 3.4197_dp, 3.4379_dp/)
        pl_index_cross_section_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                            (/2.1612_dp, 2.2001_dp, 2.2564_dp, 2.3601_dp, 2.5054_dp, &
                              2.6397_dp, 2.7642_dp, 2.8714_dp, 2.9700_dp, 3.0528_dp, &
                              3.1229_dp, 3.1892_dp, 3.2451_dp, 3.2853_dp, 3.3187_dp, &
                              3.3464_dp, 3.3640_dp, 3.3811_dp, 3.3967_dp, 3.4203_dp/)
        pl_index_cross_section_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+20) = &
                            (/2.6930_dp, 2.7049_dp, 2.7213_dp, 2.7503_dp, 2.7906_dp, &
                              2.8300_dp, 2.8711_dp, 2.9121_dp, 2.9577_dp, 3.0041_dp, &
                              3.0522_dp, 3.1069_dp, 3.1612_dp, 3.2051_dp, 3.2448_dp, &
                              3.2796_dp, 3.3027_dp, 3.3258_dp, 3.3472_dp, 3.3805_dp/)

      case (16) 
	
        pl_index_cross_section_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+16) = &
                            (/2.9884_dp, 2.9970_dp, 3.0088_dp, 3.0298_dp, 3.0589_dp, &
                              3.0872_dp, 3.1303_dp, 3.1920_dp, 3.2410_dp, 3.2765_dp, &
                              3.3107_dp, 3.3376_dp, 3.3613_dp, 3.3878_dp, 3.4078_dp, &
                              3.4343_dp/)
        pl_index_cross_section_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+16) = &
                            (/2.1612_dp, 2.2001_dp, 2.2564_dp, 2.3601_dp, 2.5054_dp, &
                              2.6397_dp, 2.8157_dp, 3.0093_dp, 3.1229_dp, 3.1892_dp, &
                              3.2451_dp, 3.2853_dp, 3.3187_dp, 3.3546_dp, 3.3811_dp, &
                              3.4157_dp/)
        pl_index_cross_section_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+16) = &
                            (/2.6930_dp, 2.7049_dp, 2.7213_dp, 2.7503_dp, 2.7906_dp, &
                              2.8300_dp, 2.8904_dp, 2.9793_dp, 3.0522_dp, 3.1069_dp, &
                              3.1612_dp, 3.2051_dp, 3.2448_dp, 3.2904_dp, 3.3258_dp, &
                              3.3740_dp/)

      case (11) 

        pl_index_cross_section_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+11) = &
                            (/2.9926_dp, 3.0088_dp, 3.0357_dp, 3.0777_dp, 3.1303_dp, &
                              3.1773_dp, 3.2292_dp, 3.2765_dp, 3.3230_dp, 3.3775_dp, &
                              3.4155_dp/)
        pl_index_cross_section_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+11) = &
                            (/2.1803_dp, 2.2564_dp, 2.3901_dp, 2.5951_dp, 2.8157_dp, &
                              2.9700_dp, 3.0976_dp, 3.1892_dp, 3.2636_dp, 3.3407_dp, &
                              3.3913_dp/)
        pl_index_cross_section_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+11) = &
                            (/2.6989_dp, 2.7213_dp, 2.7585_dp, 2.8167_dp, 2.8904_dp, &
                              2.9577_dp, 3.0345_dp, 3.1069_dp, 3.1811_dp, 3.2727_dp, &
                              3.3397_dp/)

      case (9) 

        pl_index_cross_section_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+9) = &
                            (/3.0207_dp, 3.0777_dp, 3.1303_dp, 3.1773_dp, 3.2292_dp, &
                              3.2765_dp, 3.3230_dp, 3.3775_dp, 3.4155_dp/)
        pl_index_cross_section_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+9) = &
                            (/2.3157_dp, 2.5951_dp, 2.8157_dp, 2.9700_dp, 3.0976_dp,&
                              3.1892_dp, 3.2636_dp, 3.3407_dp, 3.3913_dp/)
        pl_index_cross_section_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+9) = &
                            (/2.7377_dp, 2.8167_dp, 2.8904_dp, 2.9577_dp, 3.0345_dp,&
                              3.1069_dp, 3.1811_dp, 3.2727_dp, 3.3397_dp/)

      case (4)

        pl_index_cross_section_HI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+4) = &
                                       (/3.0465_dp, 3.1516_dp, 3.2501_dp, 3.3833_dp/)
        pl_index_cross_section_HeI(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+4) = &
                                       (/2.4431_dp, 2.8878_dp, 3.1390_dp, 3.3479_dp/)
        pl_index_cross_section_HeII(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+4) = &
                                       (/2.7735_dp, 2.9209_dp, 3.0663_dp, 3.2833_dp/)

      case (1)

        pl_index_cross_section_HI(NumBndin1+NumBndin2+1) = 3.3369_dp
        pl_index_cross_section_HeI(NumBndin1+NumBndin2+1) = 3.2681_dp
        pl_index_cross_section_HeII(NumBndin1+NumBndin2+1) = 3.2082_dp

    end select

  end subroutine setup_scalingfactors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function photrates_add (rate1, rate2)
    
    type(photrates) :: photrates_add
    type(photrates),intent(in) :: rate1, rate2
    
    photrates_add%photo_cell_HI = rate1%photo_cell_HI + rate2%photo_cell_HI
    photrates_add%photo_cell_HeI = rate1%photo_cell_HeI + rate2%photo_cell_HeI
    photrates_add%photo_cell_HeII = rate1%photo_cell_HeII + rate2%photo_cell_HeII
    photrates_add%heat_cell_HI = rate1%heat_cell_HI + rate2%heat_cell_HI
    photrates_add%heat_cell_HeI = rate1%heat_cell_HeI + rate2%heat_cell_HeI
    photrates_add%heat_cell_HeII = rate1%heat_cell_HeII + rate2%heat_cell_HeII
    photrates_add%photo_in_HI = rate1%photo_in_HI + rate2%photo_in_HI
    photrates_add%photo_in_HeI = rate1%photo_in_HeI + rate2%photo_in_HeI
    photrates_add%photo_in_HeII = rate1%photo_in_HeII + rate2%photo_in_HeII
    photrates_add%heat_in_HI = rate1%heat_in_HI + rate2%heat_in_HI
    photrates_add%heat_in_HeI = rate1%heat_in_HeI + rate2%heat_in_HeI
    photrates_add%heat_in_HeII = rate1%heat_in_HeII + rate2%heat_in_HeII
    photrates_add%photo_out_HI = rate1%photo_out_HI + rate2%photo_out_HI
    photrates_add%photo_out_HeI = rate1%photo_out_HeI + rate2%photo_out_HeI
    photrates_add%photo_out_HeII = rate1%photo_out_HeII + rate2%photo_out_HeII
    photrates_add%heat_out_HI = rate1%heat_out_HI + rate2%heat_out_HI
    photrates_add%heat_out_HeI = rate1%heat_out_HeI + rate2%heat_out_HeI
    photrates_add%heat_out_HeII = rate1%heat_out_HeII + rate2%heat_out_HeII
    photrates_add%heat = rate1%heat + rate2%heat
    photrates_add%photo_in = rate1%photo_in + rate2%photo_in
    
  end function photrates_add
  
  subroutine set_photrates_to_zero (rate1)
    
    type(photrates),intent(out) :: rate1
    
    rate1%photo_cell_HI = 0.0
    rate1%photo_cell_HeI = 0.0
    rate1%photo_cell_HeII = 0.0
    rate1%heat_cell_HI = 0.0
    rate1%heat_cell_HeI = 0.0
    rate1%heat_cell_HeII = 0.0
    rate1%photo_in_HI = 0.0
    rate1%photo_in_HeI = 0.0
    rate1%photo_in_HeII = 0.0
    rate1%heat_in_HI = 0.0
    rate1%heat_in_HeI = 0.0
    rate1%heat_in_HeII = 0.0
    rate1%photo_out_HI = 0.0
    rate1%photo_out_HeI = 0.0
    rate1%photo_out_HeII = 0.0
    rate1%heat_out_HI = 0.0
    rate1%heat_out_HeI = 0.0
    rate1%heat_out_HeII = 0.0
    rate1%heat = 0.0
    rate1%photo_in = 0.0
    
  end subroutine set_photrates_to_zero
  
end module radiation
