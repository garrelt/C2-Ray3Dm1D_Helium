!     This module contains data and routines which deal with radiative
!     effects. Its main part deal with photo-ionizing radiation, but it
!     also initializes other radiative properties, such as cooling (which
!     are contained in different modules).
!     It can be used in hydrodynamic or stand-alone radiative transfer 
!     calculations.

module radiation_sed_parameters
  
  use precision, only: dp
  use my_mpi
  use file_admin, only: logf, stdinput, file_input
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
  use sed_parameters, only: ask_for_sed, &               ! Ask for SED parameters if true
                            T_eff_nominal,&            ! Black body  effective temperature for for nominal BB SED
#ifdef PL
                              pl_index_nominal,&         ! Power law index for for nominal PL SED
                              EddLeff_nominal,&          ! Eddington efficiency for for nominal PL SED
                              EddLum, &                  ! Eddington luminosity for for nominal PL SED
                              pl_S_star_nominal, &       ! Ionizing photon rate for nominal PL SED
                              pl_MinFreq_nominal, &      ! Lowest frequency for nominal PL SED
                              pl_MaxFreq_nominal, &         ! Highest frequency for nominal PL SED
#endif
#ifdef QUASARS
                              qpl_index_nominal,&         ! Power law index for for nominal QU SED
                              qEddLeff_nominal,&          ! Eddington efficiency for for nominal QU SED
                              qEddLum, &                  ! Eddington luminosity for for nominal QU SED
                              qpl_S_star_nominal, &       ! Ionizing photon rate for nominal QU SED
                              qpl_MinFreq_nominal, &      ! Lowest frequency for nominal QU SED
                              qpl_MaxFreq_nominal, &         ! Highest frequency for nominal QU SED

#endif
                              S_star_nominal
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

#ifdef PL
  ! Power law source properties
  real(kind=dp) :: pl_index = 1.0            ! Power law index
  real(kind=dp) :: pl_minfreq           ! Minimum frequency for integration of total power
  real(kind=dp) :: pl_maxfreq           ! Maximum frequency for integration of total power
  real(kind=dp) :: pl_scaling = 1.0     ! The scaling of the flux (needs to be initialized)
  real(kind=dp) :: Edd_Efficiency = 0.0 ! Eddinton efficieny
  real(kind=dp) :: pl_S_star = 1.0
  real(kind=dp) :: pl_ionizing_luminosity
#endif

#ifdef QUASARS
  real(kind=dp) :: qpl_index = 1.0            ! Power law index
  real(kind=dp) :: qpl_minfreq           ! Minimum frequency for integration of total power
  real(kind=dp) :: qpl_maxfreq           ! Maximum frequency for integration of total power
  real(kind=dp) :: qpl_scaling = 1.0     ! The scaling of the flux (needs to be initialized)
  real(kind=dp) :: qEdd_Efficiency = 0.0 ! Eddinton efficieny
  real(kind=dp) :: qpl_S_star = 1.0
  real(kind=dp) :: qpl_ionizing_luminosity
#endif

#ifdef MPI       
    integer,private :: mympierror
#endif

contains

  ! Ask for the parameters of the spectrum
  subroutine spectrum_parms
    
    ! Ask for the input if you are processor 0 and the
    ! the ask_for_sed parameter is set in sed_parameters.
    if (ask_for_sed) then
       
       if (rank == 0) then
          do while ((sourcetype  /= 'B') .and. (sourcetype /= 'P'))
             
             ! Read source type
#if defined(QUASARS) && defined(PL)
             if (.not.file_input) &
                  write(*,"(A,$)") "Specify source type;", &
                  "Choices are blackbody source (B), power law source (P), ", &
                  "quasar power law source (Q) or all (A)"
             read(stdinput,*) sourcetype 
#elif defined(QUASARS)
             if (.not.file_input) &
                  write(*,"(A,$)") "Specify source type;", &
                  "Choices are blackbody source (B), ", &
                  "quasar power law source (Q) or both (A)"
             read(stdinput,*) sourcetype 
#elif defined(PL)
             if (.not.file_input) &
                  write(*,"(A,$)") "Specify source type;", &
                  "Choices are blackbody source (B), power law source (P), ", &
                  "or both (A)"
             read(stdinput,*) sourcetype 
#else
             sourcetype="B"  
#endif
          enddo
          
          ! Ask for SED parameters of the blackbody.
          if (sourcetype == 'B' .or. sourcetype == "A") then
             call ask_for_blackbody_parameters()
          endif
#ifdef PL
          ! Ask for SED parameters of the power law.
          if (sourcetype == 'P' .or. sourcetype == 'A') then
             call ask_for_powerlaw_parameters()
          endif
#endif
          
#ifdef QUASARS
          ! Ask for SED parameters of the quasar power law.
          if (sourcetype == 'Q' .or. sourcetype == 'A') then
             call ask_for_quasar_powerlaw_parameters()
          endif
#endif
          
       endif ! rank 0 test
       
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
#ifdef PL
       call MPI_BCAST(Edd_Efficiency,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
       call MPI_BCAST(pl_S_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
       call MPI_BCAST(pl_MinFreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
       call MPI_BCAST(pl_MaxFreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
#endif
#ifdef QUASAR
       call MPI_BCAST(qEdd_Efficiency,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
       call MPI_BCAST(qpl_S_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
       call MPI_BCAST(qpl_MinFreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
       call MPI_BCAST(qpl_MaxFreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
            mympierror)
#endif
#endif
       
    else

       ! Source properties are taken from sed_parameters
       call set_sed_parameters_from_nominal_values()

    endif

  end subroutine spectrum_parms

!-----------------------------------------------------------------------------

  subroutine set_sed_parameters_from_nominal_values()

    ! Black body source
    T_eff=T_eff_nominal
    S_star=S_star_nominal
    ! Assign some fiducial values for radius and luminosity
    ! these are scaled to correspond to S_star in routine spec_diag
    R_star=r_solar
    L_star=R_star*R_star*4.0d0*pi*sigma_SB*T_eff**4
    ! This is h/kT
    h_over_kT = hplanck/(k_B*T_eff)
    
#ifdef PL       
    ! Power law source
    pl_index=pl_index_nominal
    pl_MinFreq=pl_MinFreq_nominal
    pl_MaxFreq=pl_MaxFreq_nominal
    pl_S_star=pl_S_star_nominal
    ! Assign some fiducial values to pl_scaling and Edd_Efficiency.
    ! These are scaled to correspond to pl_S_star in routine spec_diag
    pl_scaling=1.0
    Edd_Efficiency=0.0
#endif
#ifdef QUASARS
    ! Quasar power law source
    qpl_index=qpl_index_nominal
    qpl_MinFreq=qpl_MinFreq_nominal
    qpl_MaxFreq=qpl_MaxFreq_nominal
    qpl_S_star=qpl_S_star_nominal
    ! Assign some fiducial values to qpl_scaling and qEdd_Efficiency.
    ! These are scaled to correspond to qpl_S_star in routine spec_diag
    qpl_scaling=1.0
    qEdd_Efficiency=0.0
#endif
    
  end subroutine set_sed_parameters_from_nominal_values

!-----------------------------------------------------------------------------

  subroutine ask_for_black_body_parameters ()
    
    integer :: i_choice = 0                 ! option number
    real(kind=dp) :: bb_luminosity_unscaled ! black body total flux

    write(logf,*) 'Black body source'
    
    ! In blackbody case, ask for effective temperature of the source. 
    ! The temperature should be bounded below by 2000 and above by 1000000.
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
    
    ! set some fiducial values for the BB source here, though they are
    ! not useful
    R_star=r_solar
    S_star=0.0
    L_star=0.0
    T_eff=1.0e5
  
  end subroutine ask_for_black_body_parameters

!-----------------------------------------------------------------------------

#ifdef PL
  subroutine ask_for_powerlaw_parameters ()
    
    integer :: i_choice = 0                 ! option number

    ! In power-law case, we ask for number of ionizing photons per second 
    ! or Eddinton luminosity efficiency 
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
    
    ! Convert eVs to Hz for the frequency limits
    pl_MinFreq = pl_MinFreq * ev2fr
    pl_MaxFreq = pl_MaxFreq * ev2fr
    
  end subroutine ask_for_powerlaw_parameters
#endif

!-----------------------------------------------------------------------------

#ifdef QUASARS
  subroutine ask_for_quasar_powerlaw_parameters

    integer :: i_choice = 0                 ! option number

    ! In power-law case, we ask for number of ionizing photons per second
    ! or Eddinton luminosity efficiency 

    write(logf,*) 'Quasar source'
    if (.not.file_input) then
       write(*,'(A)') 'You can specify'
       write(*,'(A)') ' (1)  Number of ionizing photons per second '
       write(*,'(A)') ' (2)  Efficiency parameter assuming a 1e6 solar mass BH'
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
       read(stdinput,*) qpl_S_star          ! Read ionizing photons per second
       write(logf,*) 'The rate is ', qpl_S_star
       
       ! Set the Eddinton luminosity efficiency to the nominal value 
       qEdd_Efficiency=qEddLeff_nominal
       qpl_scaling=1.0 ! fiducial value, to be updated in spec_diag
       
    case (2)
       write(logf,*) 'Efficiency parameter is specified'
       if (.not.file_input) write(*,'(A,$)') 'give efficiency parameter'
       read(stdinput,*) qEdd_Efficiency         ! Read Eddington efficiency
       write(logf,*) 'The efficiency parameter is ', qEdd_Efficiency
       ! Set some fiducial value, to be updated in spec_diag
       qpl_S_star=0.0
       qpl_scaling=1.0
    end select
    if (.not.file_input) write(*,'(A,$)') 'Specify power law index (for number of photons, not energy) '
    read(stdinput,*) qpl_index      ! Read power law index, this number equal to one plus that of energy 
    write(logf,*) 'Power law index is ', qpl_index
    if (.not.file_input) write(*,'(A,$)') 'give lower and upper frequency limits in eV '
    read(stdinput,*) qpl_MinFreq,qpl_MaxFreq     ! Read lower and upper frequency limits in eV   
    write(logf,*) 'The lower energy limit is ', qpl_MinFreq, ' eV'
    write(logf,*) 'The upper energy limit is ', qpl_MaxFreq, ' eV'
    
    
    ! Convert eVs to Hz for the frequency limits
    qpl_MinFreq = qpl_MinFreq * ev2fr
    qpl_MaxFreq = qpl_MaxFreq * ev2fr
    
  end subroutine ask_for_quasar_powerlaw_parameters
#endif

!-----------------------------------------------------------------------------
  
  ! This routine integrates the SEDs in order to achieve proper scaling
  ! of them. The spectrum_parms subroutine has set some of the SED properties
  ! but they are not necessarily consistent. Consistency is set here
  subroutine normalize_seds ()
    
    call normalize_blackbody()

#ifdef PL
    call normalize_powerlaw()
#endif

#ifdef QUASARS
    call normalize_quasars()
#endif

  end subroutine normalize_seds

!-----------------------------------------------------------------------------

  subroutine report_on_seds ()

    real(kind=dp) :: bb_S_star_band1, bb_S_star_band2, bb_S_star_band3
    real(kind=dp) :: pl_S_star_band1, pl_S_star_band2, pl_S_star_band3
    real(kind=dp) :: qpl_S_star_band1, qpl_S_star_band2, qpl_S_star_band3
    
    if (rank == 0) then
       if (sourcetype == 'B' .or. sourcetype == " " .or. sourcetype == "A") then
          write(logf,'(/a)')           'Using a black body with'
          write(logf,'(a,es10.3,a)')   ' Teff =       ', T_eff, ' K'
          write(logf,'(a,es10.3,a)')   ' Radius =     ', R_star/r_solar, ' R_solar'
          write(logf,'(a,es10.3,a)')   ' Luminosity = ', L_star/l_solar, ' L_solar'
          write(logf,'(a,es10.3,a)')   ' Ionzing photon rate = ', S_star, ' s^-1'
          ! Find out the number of photons in each band
          bb_S_star_band1 = integrate_sed(freq_min(NumBndin1), &
               freq_max(NumBndin1),"B","S")
          bb_S_star_band2 = integrate_sed(freq_min(NumBndin1+1), &
               freq_max(NumBndin1+NumBndin2),"B","S")
          bb_S_star_band3 = integrate_sed(freq_min(NumBndin1+NumBndin2+1), &
               freq_max(NumBndin1+NumBndin2+NumBndin3),"B","S")
          
          ! Report back to the log file
          write(logf,'(/A,(ES12.5),A//)') &
               ' Number of BB photons in band 1: ', &
               bb_S_star_band1, ' s^-1'
          write(logf,'(A,(ES12.5),A//)') &
               ' Number of BB photons in band 2: ', &
               bb_S_star_band2, ' s^-1'
          write(logf,'(A,(ES12.5),A//)') &
               ' Number of BB photons in band 3: ', &
               bb_S_star_band3, ' s^-1'
          write(logf,'(A,(ES12.5),A//)') &
               ' Total number of ionizing photons (BB): ', &
               bb_S_star_band1+bb_S_star_band2+bb_S_star_band3, 's^-1'
       endif

#ifdef PL
       if (sourcetype == 'P' .or. sourcetype == " " .or. sourcetype == "A") then
          write(logf,'(/a)')           'Using a power law source with'
          write(logf,'(a,es10.3)')   ' Power law index = ', pl_index
          write(logf,'(a,es10.3)')   ' Efficiency parameter = ', Edd_Efficiency
          write(logf,'(a,es10.3)')   ' Ionizing photon rate = ', pl_S_star
          write(logf,'(a,es10.3)')   ' Ionizing luminosity = ', pl_ionizing_luminosity
          write(logf,'(a,2f8.3,a)')   ' between energies ', &
               pl_MinFreq/(1e3*ev2fr),pl_MaxFreq/(1e3*ev2fr),' kEv'
          
          ! Find out the number of photons in each band
          pl_S_star_band1 = integrate_sed(freq_min(NumBndin1), &
               freq_max(NumBndin1),"P","S")
          pl_S_star_band2 = integrate_sed(freq_min(NumBndin1+1), &
               freq_max(NumBndin1+NumBndin2),"P","S")
          pl_S_star_band3 = integrate_sed(freq_min(NumBndin1+NumBndin2+1), &
               freq_max(NumBndin1+NumBndin2+NumBndin3),"P","S")
          
          ! Report back to the log file
          write(logf,'(A,(ES12.5),A//)') ' Number of PL photons in band 1: ', &
               pl_S_star_band1, ' s^-1'
          write(logf,'(A,(ES12.5),A//)') ' Number of PL photons in band 2: ', &
               pl_S_star_band2, ' s^-1'
          write(logf,'(A,(ES12.5),A//)') ' Number of PL photons in band 3: ', &
               pl_S_star_band3, ' s^-1'
          write(logf,'(A,(ES12.5),A//)') &
               ' Total number of ionizing photons (PL): ', &
               pl_S_star_band1+pl_S_star_band2+pl_S_star_band3, 's^-1'
       endif
#endif

#ifdef QUASARS
       if (sourcetype == 'Q' .or. sourcetype == " " .or. sourcetype == "A") then
          write(logf,'(/a)')           'Using a quasar source with'
          write(logf,'(a,es10.3)')   ' Power law index = ', qpl_index
          write(logf,'(a,es10.3)')   ' Efficiency parameter = ', qEdd_Efficiency
          write(logf,'(a,es10.3)')   ' Ionizing photon rate = ', qpl_S_star
          write(logf,'(a,es10.3)')   ' Ionizing luminosity = ',qpl_ionizing_luminosity
          write(logf,'(a,2f8.3,a)')   ' between energies ', &
               qpl_MinFreq/(1e3*ev2fr),qpl_MaxFreq/(1e3*ev2fr),' kEv'
          
          ! Find out the number of photons in each band
          qpl_S_star_band1 = integrate_sed(freq_min(NumBndin1), &
               freq_max(NumBndin1),"Q","S")
          qpl_S_star_band2 = integrate_sed(freq_min(NumBndin1+1), &
               freq_max(NumBndin1+NumBndin2),"Q","S")
          qpl_S_star_band3 = integrate_sed(freq_min(NumBndin1+NumBndin2+1), &
               freq_max(NumBndin1+NumBndin2+NumBndin3),"Q","S")
          
          ! Report back to the log file
          write(logf,'(A,(ES12.5),A//)') ' Number of Q photons in band 1: ', &
               qpl_S_star_band1, ' s^-1'
          write(logf,'(A,(ES12.5),A//)') ' Number of Q photons in band 2: ', &
               qpl_S_star_band2, ' s^-1'
          write(logf,'(A,(ES12.5),A//)') ' Number of Q photons in band 3: ', &
               qpl_S_star_band3, ' s^-1'
          write(logf,'(A,(ES12.5),A//)') &
               ' Total number of ionizing photons (Q): ', &
               qpl_S_star_band1+qpl_S_star_band2+qpl_S_star_band3, 's^-1'
       endif
#endif
    endif

end subroutine report_on_seds

!-----------------------------------------------------------------------------

  subroutine normalize_blackbody ()

    real(kind=dp) :: bb_ionizing_luminosity_unscaled
    real(kind=dp) :: S_star_unscaled, S_scaling

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
    
    ! This is R_star^2
    R_star2=R_star*R_star
    
    ! Report back to the log file
    
  end subroutine normalize_blackbody

!-----------------------------------------------------------------------------

#ifdef PL
  subroutine normalize_powerlaw ()

    real(kind=dp) :: pl_S_star_unscaled
    real(kind=dp) :: pl_S_star_wanted
    real(kind=dp) :: pl_ionizing_luminosity_unscaled
    real(kind=dp) :: pl_ionizing_luminosity_wanted

    ! Determine the scaling factor pl_scaling
    ! needed to achieve either the specified ionizing photon rate or ionizing luminosity
    if (pl_S_star > 0.0) then
       ! Total power-law ionizing photon rate is specified (photon sense)
       pl_S_star_unscaled = integrate_sed(pl_MinFreq,pl_MaxFreq,"P","S")
       pl_S_star_wanted = pl_S_star
       pl_scaling = pl_S_star_wanted/pl_S_star_unscaled
       pl_ionizing_luminosity = integrate_sed(pl_MinFreq,pl_MaxFreq,"P","L")
       Edd_Efficiency = pl_ionizing_luminosity / EddLum
    else
       ! The power-law ionizing luminosity is specified (energy sense). 
       pl_ionizing_luminosity = EddLum*Edd_Efficiency
       pl_ionizing_luminosity_unscaled = integrate_sed(pl_MinFreq,pl_MaxFreq,"P","L")
       pl_ionizing_luminosity_wanted = pl_ionizing_luminosity
       pl_scaling = pl_ionizing_luminosity_wanted/pl_ionizing_luminosity_unscaled
       pl_S_star = integrate_sed(pl_MinFreq,pl_MaxFreq,"P","S")
    endif
    
  end subroutine normalize_powerlaw
#endif

!-----------------------------------------------------------------------------

#ifdef QUASARS
  subroutine normalize_quasars ()

    real(kind=dp) :: qpl_S_star_unscaled
    real(kind=dp) :: qpl_S_star_wanted
    real(kind=dp) :: qpl_ionizing_luminosity_unscaled
    real(kind=dp) :: qpl_ionizing_luminosity_wanted

    ! Determine the scaling factor pl_scaling
    ! needed to achieve either the specified ionizing photon rate or ionizing
    ! luminosity
    if (qpl_S_star > 0.0) then
       ! Total power-law ionizing photon rate is specified (photon sense)
       qpl_S_star_unscaled = integrate_sed(qpl_MinFreq,qpl_MaxFreq,"Q","S")
       qpl_S_star_wanted = qpl_S_star
       qpl_scaling = qpl_S_star_wanted/qpl_S_star_unscaled
       qpl_ionizing_luminosity = integrate_sed(qpl_MinFreq,qpl_MaxFreq,"Q","L")
       qEdd_Efficiency = qpl_ionizing_luminosity / qEddLum
    else
       ! The power-law ionizing luminosity is specified (energy sense). 
       qpl_ionizing_luminosity = qEddLum*qEdd_Efficiency
       qpl_ionizing_luminosity_unscaled = integrate_sed(qpl_MinFreq,qpl_MaxFreq,"Q","L")
       qpl_ionizing_luminosity_wanted = qpl_ionizing_luminosity
       qpl_scaling =qpl_ionizing_luminosity_wanted/qpl_ionizing_luminosity_unscaled
       qpl_S_star = integrate_sed(qpl_MinFreq,qpl_MaxFreq,"Q","S")
    endif

  end subroutine normalize_quasars
#endif

!-----------------------------------------------------------------------------

  function integrate_sed (freq_min, freq_max, sourcetype, sedtype)
    
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
          integrand(i_freq) = blackbody_sed(frequency(i_freq),sedtype)
       enddo
       integrate_sed = 4.0*pi*R_star*R_star*scalar_romberg(integrand,weight,NumFreq,NumFreq,0)
#ifdef PL
    case("P")
       do i_freq=0,NumFreq
          ! this power-law is in number of photon sense
          integrand(i_freq)=powerlaw_sed (frequency(i_freq),pl_index,sedtype)
       enddo
       integrate_sed = pl_scaling*scalar_romberg(integrand,weight,NumFreq,NumFreq,0)
#endif
#ifdef QUASARS
    case("Q")
       do i_freq=0,NumFreq
          integrand(i_freq)=powerlaw_sed (frequency(i_freq),qpl_index,sedtype)
       enddo
       integrate_sed = qpl_scaling*scalar_romberg(integrand,weight,NumFreq,NumFreq,0)
#endif
    end select

  end function integrate_sed

!-----------------------------------------------------------------------------

  function blackbody_sed (frequency,sedtype)
    
    ! function type
    real(kind=dp) :: blackbody_sed

    real(kind=dp),intent(in) :: frequency
    character(len=1),intent(in) :: sedtype ! L or S
   
    if (frequency*h_over_kT <= 709.0_dp) then
       ! this blackbody is in number of photon sense
       blackbody_sed = two_pi_over_c_square*frequency*frequency/ &
            (exp(frequency*h_over_kT)-1.0_dp)  
       ! when the argument of the exponential function gets too high
    else
       blackbody_sed = two_pi_over_c_square*frequency*frequency/ &
         (exp((frequency*h_over_kT)/2.0_dp))/ &
         (exp((frequency*h_over_kT)/2.0_dp))
    endif
    if (sedtype == "L") blackbody_sed = hplanck*frequency* blackbody_sed

  end function blackbody_sed

!-----------------------------------------------------------------------------

  function powerlaw_sed (frequency,index,sedtype)
    
    ! function type
    real(kind=dp) :: powerlaw_sed

    real(kind=dp),intent(in) :: frequency
    real(kind=dp),intent(in) :: index
    character(len=1),intent(in) :: sedtype ! L or S
   
    ! this powerlaw is in number of photon sense
    powerlaw_sed = frequency**(-index)

    if (sedtype == "L") powerlaw_sed = hplanck * frequency * powerlaw_sed

  end function powerlaw_sed

!-----------------------------------------------------------------------------

  function double_powerlaw_sed (frequency,index1,index2,transition_frequency,sedtype)
    
    ! function type
    real(kind=dp) :: double_powerlaw_sed

    real(kind=dp),intent(in) :: frequency
    real(kind=dp),intent(in) :: index1
    real(kind=dp),intent(in) :: index2
    real(kind=dp),intent(in) :: transition_frequency
    character(len=1),intent(in) :: sedtype ! L or S
   
    ! this powerlaw is in number of photon sense
    if (frequency <= transition_frequency) then
       double_powerlaw_sed = (frequency/transition_frequency)**(-index1)
    else
       double_powerlaw_sed = (frequency/transition_frequency)**(-index2)
    endif

    if (sedtype == "L") double_powerlaw_sed = hplanck * frequency * double_powerlaw_sed

  end function double_powerlaw_sed

end module radiation_sed_parameters
