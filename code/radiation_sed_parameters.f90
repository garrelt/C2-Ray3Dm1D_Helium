!     This module contains data and routines which deal with radiative
!     effects. Its main part deal with photo-ionizing radiation, but it
!     also initializes other radiative properties, such as cooling (which
!     are contained in different modules).
!     It can be used in hydrodynamic or stand-alone radiative transfer 
!     calculations.

module radiation_sed_parameters
  
  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use file_admin, only: stdinput, file_input
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
  use c2ray_parameters, only: T_eff_nominal,&            ! Black body effective temperature for for nominal BB SED
                              bb_S_star_nominal, &       ! Ionizing photon rate for nominal BB SED
                              EddLum                     ! Eddington luminosity (for mass_nominal)
#ifdef PL
  use c2ray_parameters, only: pl_index_nominal,&         ! Power law index for for nominal PL SED
                              pl_EddLeff_nominal,&          ! Eddington efficiency for for nominal PL SED
                              pl_S_star_nominal, &       ! Ionizing photon rate for nominal PL SED
                              pl_MinFreq_nominal, &      ! Lowest frequency for nominal PL SED
                              pl_MaxFreq_nominal        ! Highest frequency for nominal PL SED
#endif
#ifdef QUASARS
 use c2ray_parameters, only: qpl_index_nominal,&         ! Power law index for for nominal QU SED
                              qpl_EddLeff_nominal,&          ! Eddington efficiency for for nominal QU SED
                              qpl_S_star_nominal, &       ! Ionizing photon rate for nominal QU SED
                              qpl_MinFreq_nominal, &      ! Lowest frequency for nominal QU SED
                              qpl_MaxFreq_nominal         ! Highest frequency for nominal QU SED

#endif
  use material, only: isothermal

  use radiation_sizes, only: NumFreq,NumFreqBnd,NumBndin1,NumBndin2,NumBndin3
  use radiation_sizes, only: freq_min,freq_max

  implicit none

! Hannah Ross: commenting this out as it doesn't seem to be used
!#ifdef PL
  ! Number of photon of the power law source
!  real(kind=dp) :: pl_input_flux = 0.0
!#endif

  ! Type of source, B=black body, P=power law source
  character :: sourcetype = " "

  ! Stellar properties
  real(kind=dp) :: T_eff  = 0.0      ! Black body effective temperature
  real(kind=dp) :: R_star = 0.0      ! Black body radius
  real(kind=dp) :: L_star = 0.0      ! Black body luminosity
  real(kind=dp) :: L_star_ion = 0.0 ! Black body ionizing luminosity
  real(kind=dp) :: bb_S_star = 0.0      ! Black body ionizing photons rate
  real(kind=dp) :: R_star2 = 0.0     ! Square of R_star
  real(kind=dp) :: h_over_kT    ! Planck constant over k_B * T_eff

  ! The lowest and highest frequency subbands used for the bb and pl source
  integer :: bb_FreqBnd_UpperLimit=NumFreqBnd
  integer :: pl_FreqBnd_LowerLimit
  integer :: pl_FreqBnd_UpperLimit

#ifdef PL
  ! Power law source properties
  real(kind=dp) :: pl_index = 1.0            ! Power law index
  real(kind=dp) :: pl_minfreq           ! Minimum frequency for integration of total power
  real(kind=dp) :: pl_maxfreq           ! Maximum frequency for integration of total power
  real(kind=dp) :: pl_scaling = 1.0     ! The scaling of the flux (needs to be initialized)
  real(kind=dp) :: pl_Edd_Efficiency = 0.0 ! Eddington efficieny
  real(kind=dp) :: pl_S_star = 1.0
#endif

#ifdef QUASARS
  real(kind=dp) :: qpl_index = 1.0            ! Power law index
  real(kind=dp) :: qpl_minfreq           ! Minimum frequency for integration of total power
  real(kind=dp) :: qpl_maxfreq           ! Maximum frequency for integration of total power
  real(kind=dp) :: qpl_scaling = 1.0     ! The scaling of the flux (needs to be initialized)
  real(kind=dp) :: qpl_Edd_Efficiency = 0.0 ! Eddington efficieny
  real(kind=dp) :: qpl_S_star = 1.0
#endif

#ifdef MPI       
    integer,private :: mympierror
#endif

contains

  ! Ask for the parameters of the spectrum
  subroutine spectrum_parms
    
    ! Ask for the input if you are processor 0 and the
    ! spectral parameters are not set in the c2ray_parameters
    ! Note that it is assumed that if teff_nominal is set, 
    ! bb_S_star_nominal is ALSO set.
    if (T_eff_nominal == 0.0) then
       
       if (rank == 0) then
          do while ((sourcetype  /= 'B') .and. (sourcetype /= 'P') .and. &
               (sourcetype /= 'Q') .and. (sourcetype /= 'A'))
             
#if defined(QUASARS) && defined(PL)
             if (.not.file_input) &
                  write(*,"(A,A,A,$)") "Specify source type; ", &
                  "Choices are blackbody source (B), power law source (P), ", &
                  "quasar (Q) or all three (A): "
             ! Read source type, either blackbody, power law or quasar source
             read(stdinput,*) sourcetype 
#elif defined(QUASARS)
             if (.not.file_input) &
                  write(*,"(A,$)") "Specify source type;", &
                  "Choices are blackbody source (B), ", &
                  "quasar (Q) or both (A): "
             ! Read source type, either blackbody or quasar
             read(stdinput,*) sourcetype 
#elif defined(PL)
             if (.not.file_input) &
                  write(*,"(A,$)") "Specify source type;", &
                  "Choices are blackbody source (B), power law source (P), ", &
                  "or both (A): "
             ! Read source type, either blackbody or power law source
             read(stdinput,*) sourcetype 
#else
             ! If not compiled for power law or quasar, source type has to be
             ! black body
             sourcetype = "B"
#endif
          enddo
          
          ! Read in the black body parameters from the command line/input file
          if (sourcetype == 'B' .or. sourcetype == "A") then
             call input_black_body_spectral_parameters()
          endif

#ifdef PL
          if (sourcetype == 'P' .or. sourcetype == 'A') then
             call input_power_law_spectral_parameters()
          endif
#endif

#ifdef QUASARS
          if (sourcetype == 'Q' .or. sourcetype == 'A') then
             call input_quasar_spectral_parameters()
          endif
#endif

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
       call MPI_BCAST(bb_S_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
#ifdef PL
       call MPI_BCAST(pl_Edd_Efficiency,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
       call MPI_BCAST(pl_S_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
       call MPI_BCAST(pl_MinFreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
       call MPI_BCAST(pl_MaxFreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
#endif
#ifdef QUASARS
       call MPI_BCAST(qpl_Edd_Efficiency,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
       call MPI_BCAST(qpl_S_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
       call MPI_BCAST(qpl_MinFreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
       call MPI_BCAST(qpl_MaxFreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
#endif
#endif
       
       ! In case source properties are taken from c2ray_parameters
    else
       ! T_eff and S_star are assumed to have been set in the c2ray_parameter module
       T_eff=T_eff_nominal
       bb_S_star=bb_S_star_nominal

       ! Assign some fiducial values, these are scaled to correspond 
       ! to S_star in routine spec_diag
       R_star=r_solar
       L_star=R_star*R_star*4.0d0*pi*sigma_SB*T_eff**4

#ifdef PL       
       ! Power law source properties set to nominal values from c2ray_parameter module
       pl_index=pl_index_nominal
       pl_MinFreq=pl_MinFreq_nominal
       pl_MaxFreq=pl_MaxFreq_nominal
       pl_S_star=pl_S_star_nominal
       ! Assign some fiducial values, these are scaled to correspond 
       ! to S_star in routine spec_diag
       pl_scaling=1.0
       pl_Edd_Efficiency=0.0
#endif

#ifdef QUASARS
       ! Power law source properties set to nominal values from c2ray_parameter
       ! module
       qpl_index=qpl_index_nominal
       qpl_MinFreq=qpl_MinFreq_nominal
       qpl_MaxFreq=qpl_MaxFreq_nominal
       qpl_S_star=qpl_S_star_nominal
       ! Assign some fiducial values, these are scaled to correspond 
       ! to S_star in routine spec_diag
       qpl_scaling=1.0
       qpl_Edd_Efficiency=0.0
#endif

    endif
    
    ! This is the often needed combination h/kT
    h_over_kT = hplanck/(k_B*T_eff)
    
  end subroutine spectrum_parms

  !-----------------------------------------------------------------------------

  subroutine input_black_body_spectral_parameters()

    integer :: i_choice = 0                 ! option number
    real(kind=dp) :: bb_luminosity_unscaled ! black body total flux
    
    ! In blackbody case, ask for effective temperature of the source. 
    ! The temperature should be bounded below by 2000 and above by 1000000.
    ! and then ask for some parameters of the blackbody.

    write(logf,*) 'Black body source'
             
    do while (T_eff < 2000.0 .or. T_eff > 1000000.) 
       if (.not.file_input) write(*,'(A,$)') 'Give black body effective temperature (2000 <= T <= 1000000): '
       read(stdinput,*) T_eff      ! Read temperature of black body
       write(logf,'(a,es10.3)') 'Temperature of the black body : ', T_eff
       if (T_eff < 2000.0 .or. T_eff > 1000000.) then
          write(*,*) 'Error: Effective temperature out of range. Try again'
       endif
    enddo
    
    ! Find total flux of blackbody (Stefan-Boltzmann law)
    bb_luminosity_unscaled = sigma_SB*T_eff*T_eff*T_eff*T_eff
    
    ! Ask for radius, luminosity, ionizing luminosity or ionizing photon rate?
    if (.not.file_input) then
       write(*,'(A)') 'For a Black Body You can specify' 
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
       write(logf,'(a,es10.3,a)') 'The radius is ', R_star, ' solar radius'
       R_star=R_star*r_solar
       L_star=4.0d0*pi*R_star*R_star*bb_luminosity_unscaled
       bb_S_star=0.0  ! Number of photo-ionizing photons is set to zero here but will be reset in spec_diag routine
       
    case (2)
       write(logf,*) 'A total luminosity is specified'
       if (.not.file_input) write(*,'(A,$)') 'Give total luminosity in solar luminosity: '
       read(stdinput,*) L_star      ! Read luminosity of the black body
       write(logf,'(a,es10.3,a)') 'The luminosity is ', L_star, ' solar luminosity'
       L_star=L_star*l_solar
       R_star=dsqrt(L_star/(4.0d0*pi*bb_luminosity_unscaled))
       bb_S_star=0.0   ! Number of photo-ionizing photons is set to zero here but will be reset in spec_diag routine
       
    case (3)
       write(logf,*) 'Total ionizing luminosity is specified'
       if (.not.file_input) write(*,'(A,$)') 'Give total ionizing luminosity in solar luminosity: '
       read(stdinput,*) L_star_ion   ! Read ionizing luminosity of the black body
       write(logf,'(A,es10.3,A)') 'Total ionizing luminosity is ', L_star_ion, ' solar luminosty'
       L_star_ion=L_star_ion*l_solar
       ! Assign some fiducial values, these are overwritten in routine spec_diag
       R_star=r_solar
       L_star=4.0d0*pi*R_star*R_star*bb_luminosity_unscaled
       bb_S_star=0.0   ! Number of photo-ionizing photons is set to zero here but will be reset in spec_diag routine
       
    case (4)
       write(logf,*) 'Rate of ionzing photons (bb_S_star) is specified'
       if (.not.file_input) write(*,'(A,$)') 'Give the number of ionizing photons per second: '
       read(stdinput,*) bb_S_star
       write(logf,*) 'The number of photons per second is ', bb_S_star
       ! Assign some fiducial values for R_star and L_star, 
       ! these are scaled to correspond to bb_S_star in routine spec_diag
       R_star=r_solar
       L_star=4.0d0*pi*R_star*R_star*bb_luminosity_unscaled
       
    end select
    
  end subroutine input_black_body_spectral_parameters

  !-----------------------------------------------------------------------------

#ifdef PL
  ! This routines handles the reading in of parameters for the power law
  ! source if they are specified interactively or in the input file.
  ! It is identical to the routine input_quasar_spectral_parameters
  ! except for the names of the variables.
  subroutine input_power_law_spectral_parameters

    integer :: i_choice

    write(logf,*) 'Power law source'

    ! Ask for power law index
    if (.not.file_input) write(*,'(A,$)') 'Specify power law index for power law source (for number of photons, not energy) '
    read(stdinput,*) pl_index      ! Read power law index, this number equal to one plus that of energy 
    write(logf,'(A,F10.3)') 'Power law index for power law source is ', pl_index

    ! In power-law case, we ask for number of ionizing photons per second
    ! or Eddington luminosity efficiency 
    if (.not.file_input) then
       write(*,'(A)') 'For a power law you can specify'
       write(*,'(A)') ' 1)  Number of ionizing photons per second '
       write(*,'(A)') ' 2)  Efficiency parameter assuming a 1e6 solar mass BH' 
    endif
    
    i_choice=0
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
       if (.not.file_input) write(*,'(A,$)') 'give number of ionizing photons per second: '
       read(stdinput,*) pl_S_star          ! Read ionizing photons per second
       write(logf,'(A,es10.3)') 'The rate is ', pl_S_star
       
       ! Set the Eddinton luminosity efficiency to the nominal value 
       pl_Edd_Efficiency=pl_EddLeff_nominal
       pl_scaling=1.0 ! fiducial value, to be updated in spec_diag
       
    case (2)
       write(logf,*) 'Efficiency parameter is specified'
       if (.not.file_input) write(*,'(A,$)') 'give efficiency parameter: '
       read(stdinput,*) pl_Edd_Efficiency         ! Read Eddington efficiency
       write(logf,'(A,es10.3)') 'The efficiency parameter is ', pl_Edd_Efficiency
       ! Set some fiducial value, to be updated in spec_diag
       pl_S_star=0.0
       pl_scaling=1.0
    end select
    
    if (.not.file_input) write(*,'(A,$)') 'For power law source specify lower and upper frequency limits in eV: '
    read(stdinput,*) pl_MinFreq,pl_MaxFreq     ! Read lower and upper frequency limits in eV	
    write(logf,'(A,F10.3,A)') 'The lower energy limit is ', pl_MinFreq, ' eV'
    write(logf,'(A,F10.3,A)') 'The upper energy limit is ', pl_MaxFreq, ' eV'
    
    ! Convert eVs to Hz for the frequency limits
    pl_MinFreq = pl_MinFreq * ev2fr
    pl_MaxFreq = pl_MaxFreq * ev2fr
    
  end subroutine input_power_law_spectral_parameters
#endif

  !-----------------------------------------------------------------------------

#ifdef QUASARS
  ! This routines handles the reading in of parameters for the quasar
  ! source if they are specified interactively or in the input file.
  ! It is identical to the routine input_power_law_spectral_parameters
  ! except for the names of the variables.
  subroutine input_quasar_spectral_parameters

    integer :: i_choice

    write(logf,*) 'Quasar source'

    ! Ask for power law index for quasar source
    if (.not.file_input) write(*,'(A,$)') 'Specify quasar spectrum power law index (for number of photons, not energy) '
    read(stdinput,*) qpl_index      ! Read power law index, this number equal to one plus that of energy 
    write(logf,'(A,F10.3)') 'Power law index for Quasar source is ', qpl_index

    ! In quasar case, we ask for number of ionizing photons per second
    ! or Eddinton luminosity efficiency 
    if (.not.file_input) then
       write(*,'(A)') 'For the Quasar source you can specify'
       write(*,'(A)') ' (1)  Number of ionizing photons per second '
       write(*,'(A)') ' (2)  Efficiency parameter assuming a 1e6 solar mass BH'
    endif
    
    i_choice=0
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
       if (.not.file_input) write(*,'(A,$)') 'give number of ionizing photons per second: '
       read(stdinput,*) qpl_S_star          ! Read ionizing photons per second
       write(logf,'(A,es10.3)') 'The rate is ', qpl_S_star
       
       ! Set the Eddinton luminosity efficiency to the nominal value 
       qpl_Edd_Efficiency=qpl_EddLeff_nominal
       qpl_scaling=1.0 ! fiducial value, to be updated in spec_diag
       
    case (2)
       write(logf,*) 'Efficiency parameter is specified'
       if (.not.file_input) write(*,'(A,$)') 'give efficiency parameter: '
       read(stdinput,*) qpl_Edd_Efficiency         ! Read Eddington efficiency
       write(logf,'(A,es10.3)') 'The efficiency parameter is ', qpl_Edd_Efficiency
       ! Set some fiducial value, to be updated in spec_diag
       qpl_S_star=0.0
       qpl_scaling=1.0
    end select

    if (.not.file_input) write(*,'(A,$)') 'For the quasar source specify lower and upper frequency limits in eV: '
    read(stdinput,*) qpl_MinFreq,qpl_MaxFreq     ! Read lower and upper frequency limits in eV   
    write(logf,'(A,F10.3,A)') 'The lower energy limit is ', qpl_MinFreq, ' eV'
    write(logf,'(A,F10.3,A)') 'The upper energy limit is ', qpl_MaxFreq, ' eV'
    
    ! Convert eVs to Hz for the frequency limits
    qpl_MinFreq = qpl_MinFreq * ev2fr
    qpl_MaxFreq = qpl_MaxFreq * ev2fr
    
  end subroutine input_quasar_spectral_parameters
#endif

  !-----------------------------------------------------------------------------
  
  ! Determine spectrum diagnostics
  ! This routine integrates the SEDs in order to achieve proper scaling
  ! of them. The spectrum_parms subroutine has set some of the SED properties
  ! but they are not necessarily consistent. Consistency is set here
  subroutine spec_diag ()
    
    real(kind=dp) :: bb_ionizing_luminosity_unscaled
    real(kind=dp) :: S_star_unscaled, S_scaling

#ifdef PL 
    real(kind=dp) :: pl_S_star_unscaled
    real(kind=dp) :: pl_S_star_wanted
    real(kind=dp) :: pl_ionizing_luminosity_unscaled
    real(kind=dp) :: pl_ionizing_luminosity
    real(kind=dp) :: pl_ionizing_luminosity_wanted
#endif  
#ifdef QUASARS
    real(kind=dp) :: qpl_S_star_unscaled
    real(kind=dp) :: qpl_S_star_wanted
    real(kind=dp) :: qpl_ionizing_luminosity_unscaled
    real(kind=dp) :: qpl_ionizing_luminosity
    real(kind=dp) :: qpl_ionizing_luminosity_wanted
#endif  

    ! Set scaling for black body source if this sourcetype is active
    if (sourcetype == "B" .or. sourcetype == " " .or. sourcetype == "A") then
       ! Determine the scaling factor for the black body source
       call set_scaling_factor("B")
       
       ! This is R_star^2 for the black body source
       R_star2=R_star*R_star
       
       ! Report information about the three frequency bands
       call report_source_band_information("B")
    endif

#ifdef PL
    ! Set scaling for power law source if this sourcetype is active
    if (sourcetype == "P" .or. sourcetype == " " .or. sourcetype == "A") then
       ! Determine the scaling factor pl_scaling
       call set_scaling_factor("P")

       ! Report information about the three frequency bands
       call report_source_band_information("P")
    endif
#endif


#ifdef QUASARS
    ! Set scaling for quasar source if this sourcetype is active
    if (sourcetype == "Q" .or. sourcetype == " " .or. sourcetype == "A") then
       ! Determine the scaling factor qpl_scaling
       call set_scaling_factor("Q")
       
       ! Report information about the three frequency bands
       call report_source_band_information("Q")
    endif
#endif

  end subroutine spec_diag
  
  !-----------------------------------------------------------------------------

  ! Determine the scaling factor needed to achieve either the specified 
  ! ionizing photon rate (S_star) or ionizing luminosity (L_ionizing)

  subroutine set_scaling_factor(sourcetype)

    character,intent(in) :: sourcetype
    
    real(kind=dp) :: S_star
    real(kind=dp) :: S_star_unscaled
    real(kind=dp) :: S_star_wanted
    real(kind=dp) :: L_ionizing
    real(kind=dp) :: L_ionizing_unscaled
    real(kind=dp) :: L_ionizing_wanted
    real(kind=dp) :: Edd_Efficiency
    character(len=10) :: source_string
    real(kind=dp) :: index
    real(kind=dp) :: MinFreq
    real(kind=dp) :: MaxFreq
    real(kind=dp) :: scaling

    scaling=0.0
    select case (sourcetype)
    case("B") 
       source_string="Black Body"
       S_star = bb_S_star
       L_ionizing = L_star_ion
       MinFreq = freq_min(NumBndin1)
       MaxFreq = freq_max(NumBndin1+NumBndin2+NumBndin3)
#ifdef PL
    case("P") 
       source_string="power law"
       index = pl_index
       Edd_Efficiency = pl_Edd_Efficiency
       S_star = pl_S_star
       MinFreq = pl_MinFreq
       MaxFreq = pl_MaxFreq
#endif
#ifdef QUASARS
    case("Q") 
       source_string="quasar"
       index = qpl_index
       Edd_Efficiency = qpl_Edd_Efficiency
       S_star = qpl_S_star
       MinFreq = qpl_MinFreq
       MaxFreq = qpl_MaxFreq
#endif
    end select

    if (S_star > 0.0) then
       ! Total power-law ionizing photon rate is specified (photon sense)
       S_star_unscaled = integrate_sed(MinFreq,MaxFreq,sourcetype,"S")
       S_star_wanted = S_star
       scaling = S_star_wanted/S_star_unscaled
    else
       ! The power-law ionizing luminosity is specified (energy sense). 
       if (sourcetype == "P" .or. sourcetype == "Q") &
            L_ionizing = EddLum*Edd_Efficiency
       if (L_ionizing > 0.0) then
          L_ionizing_unscaled = integrate_sed(MinFreq,MaxFreq,sourcetype,"L")
          L_ionizing_wanted = L_ionizing
          scaling = L_ionizing_wanted/L_ionizing_unscaled
       else
          ! Only for black body: R_star has already been set, so scaling is 1
          scaling = 1.0
       endif
    endif
    
    ! Put the scaling factor into the right variable for different source types
    select case (sourcetype)
    case("B")
       R_star = sqrt(scaling)*R_star
       L_star = R_star*R_star*4.0_dp*pi*sigma_SB*T_eff**4
#ifdef PL       
    case("P") 
       pl_scaling = scaling
#endif
#ifdef QUASARS
    case("Q") 
       qpl_scaling = scaling
#endif
    end select

    ! Now that scaling is set, recalculate the various quantities for
    ! diagnostic reasons. These should reproduce whatever was needed.
    S_star = integrate_sed(MinFreq,MaxFreq,sourcetype,"S")
    L_ionizing = integrate_sed(MinFreq,MaxFreq,sourcetype,"L")
    Edd_Efficiency = L_ionizing / EddLum

    ! Save the ionizing photon rate (S_star) for each of the source types
    select case (sourcetype)
    case("B")
       bb_S_star = S_star
#ifdef PL
    case("P") 
       pl_S_star = S_star
       pl_Edd_Efficiency=Edd_Efficiency
#endif
#ifdef QUASARS
    case("Q") 
       qpl_S_star = S_star
       qpl_Edd_Efficiency=Edd_Efficiency
#endif
    end select

    ! Report back to the log file
    if (rank == 0) then
       write(logf,'(/a,a,a)')     "Using a ",source_string," source with"
       select case (sourcetype)
       case("P","Q")
          write(logf,'(a,es10.3)')   ' Power law index = ', index
          write(logf,'(a,es10.3)')   ' Efficiency parameter = ', Edd_Efficiency
       case("B")
          write(logf,'(a,es10.3,a)')   ' Teff =       ', T_eff, ' K'
          write(logf,'(a,es10.3,a)')   ' Radius =     ', R_star/r_solar, ' R_solar'
          write(logf,'(a,es10.3,a)')   ' Luminosity = ', L_star/l_solar, ' L_solar'
       end select
       write(logf,'(a,es10.3,a)')   ' Ionizing photon rate = ', S_star, ' s^-1'
       write(logf,'(a,es10.3,a)')   ' Ionizing luminosity = ', L_ionizing/l_solar, ' L_solar'
       write(logf,'(a,2f8.3,a)')   ' between energies ', &
            MinFreq/(1e3*ev2fr),MaxFreq/(1e3*ev2fr),' kEv'
    endif

  end subroutine set_scaling_factor

  !-----------------------------------------------------------------------------

  subroutine report_source_band_information(sourcetype)
    
    character,intent(in) :: sourcetype

    real(kind=dp) :: S_star_band1, S_star_band2, S_star_band3
    real(kind=dp) :: L_star_band1, L_star_band2, L_star_band3
    real(kind=dp) :: MinFreq, MaxFreq
    real(kind=dp) :: freq_min_Bnd1, freq_min_Bnd2, freq_min_Bnd3
    real(kind=dp) :: freq_max_Bnd1, freq_max_Bnd2, freq_max_Bnd3
    character(len=3) :: source_string

    select case (sourcetype)
    case("B") 
       source_string="BB"
       MinFreq = freq_min(NumBndin1)
#ifdef PL
    case("P") 
       source_string="PL"
       MinFreq = pl_MinFreq
       MaxFreq = pl_MaxFreq
#endif
#ifdef QUASARS
    case("Q") 
       source_string="QSO"
       MinFreq = qpl_MinFreq
       MaxFreq = qpl_MaxFreq
#endif
    end select

    ! Take care of frequency limits
    freq_min_Bnd1=max(freq_min(NumBndin1),MinFreq)
    freq_max_Bnd1=min(freq_max(NumBndin1),MaxFreq)
    freq_min_Bnd2=max(freq_min(NumBndin1+1),MinFreq)
    freq_max_Bnd2=min(freq_max(NumBndin1+NumBndin2),MaxFreq)
    freq_min_Bnd3=max(freq_min(NumBndin1+NumBndin2+1),MinFreq)
    freq_max_Bnd3=min(freq_max(NumBndin1+NumBndin2+NumBndin3),MaxFreq)

    ! Find out the number of photons in each band
    ! Band 1
    if (freq_min_Bnd1 < freq_max(NumBndin1)) then
       S_star_band1 = integrate_sed(freq_min_Bnd1, &
            freq_max_Bnd1,sourcetype,"S")
       L_star_band1 = integrate_sed(freq_min_Bnd1, &
            freq_max_Bnd1,sourcetype,"L")
    else
       S_star_band1 = 0.0
       L_star_band1 = 0.0
    endif
    
    ! Band 2
    if (freq_min_Bnd2 < freq_max(NumBndin1+NumBndin2)) then
       S_star_band2 = integrate_sed(freq_min_Bnd2, &
            freq_max_Bnd2,sourcetype,"S")
       L_star_band2 = integrate_sed(freq_min_Bnd2, &
            freq_max_Bnd2,sourcetype,"L")
    else
       S_star_band2 = 0.0
       L_star_band2 = 0.0
    endif

    ! Band 3
    if (freq_min_Bnd3 < freq_max(NumBndin1+NumBndin2+NumBndin3)) then
       S_star_band3 = integrate_sed(freq_min_Bnd3, &
            freq_max_Bnd3,sourcetype,"S")
       L_star_band3 = integrate_sed(freq_min_Bnd3, &
            freq_max_Bnd3,sourcetype,"L")
    else
       S_star_band3 = 0.0
       L_star_band3 = 0.0
    endif
    
    ! Report back to the log file
    if (rank == 0) then
       write(logf,'(A,A,A,(ES12.5),A//)') " Number of ",source_string," photons in band 1: ", &
            S_star_band1, " s^-1"
       write(logf,'(A,A,(ES12.5),A//)') source_string," luminosity in band 1: ", &
            L_star_band1, " erg s^-1"
       write(logf,'(A,A,A,(ES12.5),A//)') " Number of ",source_string," photons in band 2: ", &
            S_star_band2, " s^-1"
       write(logf,'(A,A,(ES12.5),A//)') source_string," luminosity in band 2: ", &
            L_star_band2, " erg s^-1"
       write(logf,'(A,A,A(ES12.5),A//)') " Number of ",source_string," photons in band 3: ", &
            S_star_band3, " s^-1"
       write(logf,'(A,A,(ES12.5),A//)') source_string," luminosity in band 3: ", &
            S_star_band3, " erg s^-1"
       write(logf,'(A,A,A,(ES12.5),A//)') &
            " Total number of ionizing photons in band 1,2,3 (",source_string,"): ", &
            S_star_band1+S_star_band2+S_star_band3, "s^-1"
       write(logf,'(A,A,A,(ES12.5),A,(ES12.5),A,//)') &
            " Total luminosity in band 1,2,3 (",source_string,"): ", &
            L_star_band1+L_star_band2+L_star_band3, "erg s^-1 = ", &
            (L_star_band1+L_star_band2+L_star_band3)/L_solar, " L0"
       
    endif
  end subroutine report_source_band_information

  !-----------------------------------------------------------------------------

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

    ! Set default answer
    integrate_sed=0.0

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
#ifdef PL
    case("P")
       do i_freq=0,NumFreq
          ! this power-law is in number of photon sense
          integrand(i_freq)=frequency(i_freq)**(-pl_index)
          if (sedtype == "L") integrand(i_freq) = hplanck*frequency(i_freq)*integrand(i_freq)
       enddo
       integrate_sed = pl_scaling*scalar_romberg(integrand,weight,NumFreq,NumFreq,0)
#endif
#ifdef QUASARS
    case("Q")
       do i_freq=0,NumFreq
          ! this power-law is in number of photon sense
          integrand(i_freq)=frequency(i_freq)**(-qpl_index)
          if (sedtype == "L") integrand(i_freq) = hplanck*frequency(i_freq)*integrand(i_freq)
       enddo
       integrate_sed = qpl_scaling*scalar_romberg(integrand,weight,NumFreq,NumFreq,0)
#endif
    end select

  end function integrate_sed

end module radiation_sed_parameters
