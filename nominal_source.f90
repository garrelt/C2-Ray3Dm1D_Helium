module nominal_source
  
  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use mathconstant, only: pi
  use cgsconstants, only: sigma_SB, &                    ! Stefan-Boltzmann constant
                          hplanck, &                     ! Planck constant
                          k_B, &                         ! Boltzmann constant
                          two_pi_over_c_square           ! two times pi over c aquare
  use cgsphotoconstants, only: ion_freq_HI,&             ! HI ionization energy in frequency
                               ion_freq_HeI,&            ! HeI ionization energy in frequency
                               ion_freq_HeII
  use astroconstants, only: r_solar, &                   ! Solar radius
                            l_solar                      ! Solar luminosity
  use romberg, only: scalar_romberg                      ! 1D integration subroutine
  use c2ray_parameters, only: bb_S_star_nominal, &          ! Ionizing photon rate for for nominal SED
                              pl_index_nominal,&         ! Power law index for for nominal SED
                              pl_S_star_nominal       ! Ionizing photon rate for nominal PL SED
  use parameter, only: NumFreq, NumBndin1, NumBndin2, NumBndin3, NumFreqBnd, &
                       R_star_nominal, L_star_nominal, pl_scaling, input_source_temperature
  use array, only: freq_max, freq_min

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Determine spectrum diagnostics   
  subroutine nominal_source_initialization()

    integer :: i_freq
    real(kind=dp) :: h_over_kT, freq_step 
    real(kind=dp) :: bb_photo_flux_unscaled, pl_photo_flux_unscaled, pl_photo_flux_wanted
    real(kind=dp) :: bb_S_star_unscaled, bb_scaling
    real(kind=dp) :: bb_S_star_band1, bb_S_star_band2, bb_S_star_band3
    real(kind=dp) :: pl_S_star_band1, pl_S_star_band2, pl_S_star_band3
    real(kind=dp) :: bb_S_star_band123
    real(kind=dp) :: pl_S_star_band123
    real(kind=dp), dimension(0:NumFreq) :: frequency, weight
    real(kind=dp), dimension(0:NumFreq) :: bb_photon, bb_energy
    real(kind=dp), dimension(0:NumFreq) :: pl_photon, pl_energy
    real(kind=dp) :: sigma_T_4

    if (rank .eq. 0) then
      write(logf,*) "Beginning of nominal source initialization" 
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Black-Body !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    sigma_T_4=sigma_SB*input_source_temperature*input_source_temperature*input_source_temperature*input_source_temperature
   
    ! This is h/kT
    h_over_kT = hplanck/(k_B*input_source_temperature)

    ! Frequency step width
    freq_step = (freq_max(NumFreqBnd)-freq_min(1))/real(NumFreq)

    ! Fill the arrays for integration(frequency, weight, spectrum)
    do i_freq=0,NumFreq
      frequency(i_freq) = freq_min(1)+freq_step*real(i_freq)
      weight(i_freq) = freq_step
    enddo

    do i_freq=0,NumFreq
      if (frequency(i_freq)*h_over_kT .le. 709.0_dp)then
        ! this blackbody is in number of photon sense
        bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                            (exp(frequency(i_freq)*h_over_kT)-1.0_dp)  
        ! when the argument of the exponential function gets too high
      else
        bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                            (exp((frequency(i_freq)*h_over_kT)/2.0_dp))/(exp((frequency(i_freq)*h_over_kT)/2.0_dp))
      endif     
    enddo
    ! Black-body flux (photon sense)
    bb_photo_flux_unscaled = scalar_romberg(bb_photon,weight,NumFreq,NumFreq,0) 
    ! Find out S_star with the given radius.
    bb_S_star_unscaled = 4.0*pi*bb_photo_flux_unscaled 
    bb_scaling = bb_S_star_nominal/bb_S_star_unscaled
    R_star_nominal = sqrt(bb_scaling)
    L_star_nominal = 4.0d0*pi*sigma_T_4*bb_scaling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Power-Law !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Find out the power law scaling factor for the case Eddinton luminosity efficiency is provided.
    do i_freq=0,NumFreq
      ! this power-law is in number of photon sense
      pl_photon(i_freq) = frequency(i_freq)**(-pl_index_nominal)       
    enddo
    ! power law flux (photon sense) 
    pl_photo_flux_unscaled = scalar_romberg(pl_photon,weight,NumFreq,NumFreq,0)
    pl_scaling = pl_S_star_nominal/pl_photo_flux_unscaled

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Find out the number of photons in band 1
    freq_step=(freq_max(NumBndin1)-freq_min(NumBndin1))/real(NumFreq)

    do i_freq=0,NumFreq
      frequency(i_freq)=ion_freq_HI+freq_step*real(i_freq)
      weight(i_freq)=freq_step
      if (frequency(i_freq)*h_over_kT .le. 709.0_dp) then
        ! this blackbody is in number of photon sense
        bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                            (exp(frequency(i_freq)*h_over_kT)-1.0_dp)  
      ! when the argument of the exponential function gets too high
      else
        bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                            (exp((frequency(i_freq)*h_over_kT)/2.0_dp))/ &
                            (exp((frequency(i_freq)*h_over_kT)/2.0_dp))
      endif
    enddo
    bb_photo_flux_unscaled = scalar_romberg(bb_photon,weight,NumFreq,NumFreq,0)
    bb_S_star_band1 = 4.0*pi*R_star_nominal*R_star_nominal*bb_photo_flux_unscaled

    do i_freq=0,NumFreq
      ! this power-law is in number of photon sense
      pl_photon(i_freq)=frequency(i_freq)**(-pl_index_nominal)
    enddo
    pl_photo_flux_unscaled = scalar_romberg(pl_photon,weight,NumFreq,NumFreq,0)
    pl_S_star_band1 = pl_photo_flux_unscaled*pl_scaling

    ! Find out the number of photons in band 2
    freq_step = (freq_max(NumBndin1+NumBndin2)-freq_min(NumBndin1+1))/real(NumFreq)

    do i_freq=0,NumFreq
      frequency(i_freq) = ion_freq_HeI+freq_step*real(i_freq)
      weight(i_freq) = freq_step
      if (frequency(i_freq)*h_over_kT .le. 709.0_dp) then
        ! this blackbody is in number of photon sense
        bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                            (exp(frequency(i_freq)*h_over_kT)-1.0_dp)  
        ! when the argument of the exponential function gets too high
      else
        bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                            (exp((frequency(i_freq)*h_over_kT)/2.0_dp))/ &
                            (exp((frequency(i_freq)*h_over_kT)/2.0_dp))
      endif
      ! this power-law is in number of photon sense
      pl_photon(i_freq) = frequency(i_freq)**(-pl_index_nominal)
    enddo

    bb_photo_flux_unscaled = scalar_romberg(bb_photon,weight,NumFreq,NumFreq,0)
    pl_photo_flux_unscaled = scalar_romberg(pl_photon,weight,NumFreq,NumFreq,0)
    bb_S_star_band2 = 4.0*pi*R_star_nominal*R_star_nominal*bb_photo_flux_unscaled
    pl_S_star_band2 = pl_photo_flux_unscaled*pl_scaling
 
    ! Find out the number of photons in band 3
    freq_step = (freq_max(NumBndin1+NumBndin2+NumBndin3)-freq_min(NumBndin1+NumBndin2+1))/real(NumFreq)
    do i_freq=0,NumFreq
      frequency(i_freq) = ion_freq_HeII+freq_step*real(i_freq)
      weight(i_freq) = freq_step
      if (frequency(i_freq)*h_over_kT .le. 709.0_dp) then
        ! this blackbody is in number of photon sense
        bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                            (exp(frequency(i_freq)*h_over_kT)-1.0_dp)  
        ! when the argument of the exponential function gets too high
      else
        bb_photon(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/ &
                            (exp((frequency(i_freq)*h_over_kT)/2.0_dp))/ &
                            (exp((frequency(i_freq)*h_over_kT)/2.0_dp))
      endif
       ! this power-law is in number of photon sense
       pl_photon(i_freq) = frequency(i_freq)**(-pl_index_nominal)
    enddo

    bb_photo_flux_unscaled = scalar_romberg(bb_photon,weight,NumFreq,NumFreq,0)
    pl_photo_flux_unscaled = scalar_romberg(pl_photon,weight,NumFreq,NumFreq,0)
    bb_S_star_band3 = 4.0*pi*R_star_nominal*R_star_nominal*bb_photo_flux_unscaled
    pl_S_star_band3 = pl_photo_flux_unscaled*pl_scaling

    ! Find out the number of photons in all bands
    bb_S_star_band123 = bb_S_star_band1+bb_S_star_band2+bb_S_star_band3
    pl_S_star_band123 = pl_S_star_band1+pl_S_star_band2+pl_S_star_band3

    ! Report back to the log file
    if (rank == 0) then

      write(logf,*) "The nominal black body source:"
      write(logf,*) "Teff = ",input_source_temperature,"K"
      write(logf,*) "Radius = ",R_star_nominal/r_solar,"R_solar"
      write(logf,*) "Luminosity = ",L_star_nominal/l_solar,"L_solar"
      write(logf,*) "Ionzing photon production rate = ",bb_S_star_nominal," s^-1"

      write(logf,*) "BB Number of photons in band 1: ",bb_S_star_band1,"s^-1"
      write(logf,*) "BB Number of photons in band 2: ",bb_S_star_band2,"s^-1"
      write(logf,*) "BB Number of photons in band 3: ",bb_S_star_band3,"s^-1"
      write(logf,*) "BB Total number of ionizing photons: ",bb_S_star_band123,"s^-1" 
      write(logf,*) "BB ionizing photons error: ",(bb_S_star_band123-bb_S_star_nominal)*100.0/bb_S_star_nominal,"%"

      write(logf,*) "Using a power law source with"
      write(logf,*) "Power law index = ",pl_index_nominal
      write(logf,*) "Ionzing photon production rate = ",pl_S_star_nominal,"s^-1"

      write(logf,*) "PL Number of photons in band 1: ",pl_S_star_band1,"s^-1"
      write(logf,*) "PL Number of photons in band 2: ",pl_S_star_band2,"s^-1"
      write(logf,*) "PL Number of photons in band 3: ",pl_S_star_band3,"s^-1"
      write(logf,*) "PL Total number of ionizing photons: ",pl_S_star_band123,"s^-1" 
      write(logf,*) "PL ionizing photons error: ",(pl_S_star_band123-pl_S_star_nominal)*100.0/pl_S_star_nominal,"%"

      write(logf,*) "End of nominal source initialization" 
      write(logf,*)      
      flush(logf)
    
    endif

  end subroutine nominal_source_initialization
  
end module nominal_source
