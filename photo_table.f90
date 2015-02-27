module photo_table

  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use mathconstant, only: pi
  use cgsconstants, only: hplanck, &                     ! Planck constant
                          k_B, &                         ! Boltzmann constant
                          two_pi_over_c_square           ! two times pi over c aquare
  use cgsphotoconstants, only: ion_freq_HI,&             ! HI ionization energy in frequency
                               ion_freq_HeI,&            ! HeI ionization energy in frequency
                               ion_freq_HeII           ! HeII ionization energy in frequency
  use romberg, only: vector_romberg                      ! 1D integration subroutine
  use c2ray_parameters, only: pl_index_nominal
  use parameter, only: NumFreq, NumTau, NumBndin1, NumBndin2, NumBndin3, NumFreqBnd, NumheatBin,&
                       minlogtau, maxlogtau, dlogtau, R_star_nominal, &
                       pl_scaling, do_thermal, input_source_temperature
  use array, only: delta_freq, freq_min, &
                   pl_index_cross_section_HI, &
                   pl_index_cross_section_HeI, &
                   pl_index_cross_section_HeII, &
                   bb_photo_thick_integrand, bb_photo_thin_integrand, &
                   pl_photo_thick_integrand, pl_photo_thin_integrand, &
                   bb_heat_thick_integrand_HI, bb_heat_thick_integrand_HeI, bb_heat_thick_integrand_HeII, &
                   bb_heat_thin_integrand_HI, bb_heat_thin_integrand_HeI, bb_heat_thin_integrand_HeII, &
                   pl_heat_thick_integrand_HI, pl_heat_thick_integrand_HeI, pl_heat_thick_integrand_HeII, &
                   pl_heat_thin_integrand_HI, pl_heat_thin_integrand_HeI, pl_heat_thin_integrand_HeII, &
                   bb_photo_thick_table, bb_photo_thin_table, &
                   pl_photo_thick_table, pl_photo_thin_table, &
                   bb_heat_thick_table, bb_heat_thin_table, &
                   pl_heat_thick_table, pl_heat_thin_table

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Generate photoionization tables and heating tables   
  subroutine photo_table_initialization()

    integer :: i_freq, i_tau, i_subband
    real(kind=dp) :: R_star2, h_over_kT
    real(kind=dp), dimension(0:NumTau) :: tau
    real(kind=dp), dimension(0:NumTau) :: answer
    real(kind=dp), dimension(0:NumFreq) :: exponent_HI
    real(kind=dp), dimension(0:NumFreq) :: exponent_HeI
    real(kind=dp), dimension(0:NumFreq) :: exponent_HeII
    real(kind=dp), dimension(0:NumFreq) :: frequency
    real(kind=dp), dimension(0:NumFreq, 0:NumTau) :: weight

    if (rank .eq. 0) then
      write(logf,*) "Beginning of photo table initialization"
      write(logf,*) "Frequency resolution = ",NumFreq 
      write(logf,*) "Optical depth resolution = ",NumTau 
    endif

    ! Photoionization integrand as a function of frequency and tau
    allocate(bb_photo_thick_integrand(0:NumFreq, 0:NumTau))    
    allocate(bb_photo_thin_integrand(0:NumFreq, 0:NumTau)) 
    allocate(pl_photo_thick_integrand(0:NumFreq, 0:NumTau))    
    allocate(pl_photo_thin_integrand(0:NumFreq, 0:NumTau)) 

    if (do_thermal) then
    ! Heating integrand as a function of frequency and tau
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

    ! Photoionization table as a function of photo sub-bin and tau
    allocate(bb_photo_thick_table(0:NumTau, 1:NumFreqBnd))
    allocate(bb_photo_thin_table(0:NumTau, 1:NumFreqBnd))
    allocate(pl_photo_thick_table(0:NumTau, 1:NumFreqBnd))
    allocate(pl_photo_thin_table(0:NumTau, 1:NumFreqBnd))

    if (do_thermal) then
    ! Heating table as a function of heating sub-bin and tau
	allocate(bb_heat_thick_table(0:NumTau, 1:NumheatBin))
    	allocate(bb_heat_thin_table(0:NumTau, 1:NumheatBin))
	allocate(pl_heat_thick_table(0:NumTau, 1:NumheatBin))
    	allocate(pl_heat_thin_table(0:NumTau, 1:NumheatBin))
    endif

    ! This is h/kT
    h_over_kT=hplanck/(k_B*input_source_temperature)

    ! This is R_star^2
    R_star2=R_star_nominal*R_star_nominal

    ! fill the optical depth array used to fill the tables 
    ! it is filled in NumTau logarithmic steps 
    ! from minlogtau to maxlogtau
    do i_tau = 1,NumTau
       tau(i_tau) = 10.0**(minlogtau+dlogtau*real(i_tau-1))
    enddo

    ! Position zero corresponds to zero optical depth
    tau(0)=0.0

    ! In frequency band 1, fill in integrands and make tables

    ! Go through all the sub-bin in band 1
    do i_subband=1,NumBndin1  
     
      ! Assign values to exponent array
      do i_freq=0,NumFreq

        frequency(i_freq) = freq_min(i_subband)+delta_freq(i_subband)*real(i_freq)

        exponent_HI(i_freq) = ((frequency(i_freq)/freq_min(i_subband))**(-pl_index_cross_section_HI(i_subband)))        


      enddo

      ! Loop through the tau partition
      do i_tau=0,NumTau 

        ! Loop through the frequency partition
        do i_freq=0,NumFreq
          weight(i_freq,i_tau) = delta_freq(i_subband)  

          ! Assign values to the photo integrands
          if (tau(i_tau)*exponent_HI(i_freq) < 700.0) then   
            bb_photo_thick_integrand(i_freq,i_tau) = 4.0_dp*pi*R_star2*two_pi_over_c_square*frequency(i_freq)* &
                                                     frequency(i_freq)*exp(-tau(i_tau)*exponent_HI(i_freq))/ &
                                                     (exp(frequency(i_freq)*h_over_kT)-1.0)
            bb_photo_thin_integrand(i_freq,i_tau) = 4.0_dp*pi*R_star2*two_pi_over_c_square*frequency(i_freq)* &
                                                    frequency(i_freq)*exponent_HI(i_freq)*exp(-tau(i_tau)* &
                                                    exponent_HI(i_freq))/(exp(frequency(i_freq)*h_over_kT)-1.0)
            pl_photo_thick_integrand(i_freq,i_tau) = pl_scaling*frequency(i_freq)**(-pl_index_nominal)* &
                                                     exp(-tau(i_tau)*exponent_HI(i_freq))
            pl_photo_thin_integrand(i_freq,i_tau) = pl_scaling*frequency(i_freq)**(-pl_index_nominal)*exponent_HI(i_freq) &
                                                      *exp(-tau(i_tau)*exponent_HI(i_freq))
          else
            bb_photo_thick_integrand(i_freq,i_tau) = 0.0
            bb_photo_thin_integrand(i_freq,i_tau) = 0.0
            pl_photo_thick_integrand(i_freq,i_tau) = 0.0
            pl_photo_thin_integrand(i_freq,i_tau) = 0.0
          endif

          ! Assign values to the heating integrands
            bb_heat_thick_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                       bb_photo_thick_integrand(i_freq,i_tau)
            bb_heat_thin_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                      bb_photo_thin_integrand(i_freq,i_tau)
            pl_heat_thick_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                       pl_photo_thick_integrand(i_freq,i_tau)
            pl_heat_thin_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                      pl_photo_thin_integrand(i_freq,i_tau)

        enddo

      enddo

      ! Make photo tables
      call vector_romberg (bb_photo_thick_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      bb_photo_thick_table(:,1) = answer
      call vector_romberg (bb_photo_thin_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      bb_photo_thin_table(:,1) = answer
      call vector_romberg (pl_photo_thick_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      pl_photo_thick_table(:,1) = answer
      call vector_romberg (pl_photo_thin_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      pl_photo_thin_table(:,1) = answer

      ! Make heating tables
        call vector_romberg (bb_heat_thick_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thick_table(:,1) = answer
        call vector_romberg (bb_heat_thin_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thin_table(:,1) = answer
        call vector_romberg (pl_heat_thick_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thick_table(:,1) = answer
        call vector_romberg (pl_heat_thin_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thin_table(:,1) = answer

    enddo

    ! In frequency band 2, fill in integrands and make tables

    ! Go through all the sub-bin in band 2
    do i_subband=NumBndin1+1,NumBndin1+NumBndin2
       
      ! Assign values to exponent array
      do i_freq=0,NumFreq

        frequency(i_freq) = freq_min(i_subband)+delta_freq(i_subband)*real(i_freq)

        exponent_HeI(i_freq) = ((frequency(i_freq)/freq_min(i_subband))**(-pl_index_cross_section_HeI(i_subband)))

      enddo
         
      ! Loop through the tau partition
      do i_tau=0,NumTau 

        ! Loop through the frequency partition
        do i_freq=0,NumFreq
          weight(i_freq,i_tau) = delta_freq(i_subband)  
              
          ! Assign values to the photo integrands
          if (tau(i_tau)*exponent_HeI(i_freq) < 700.0) then 
            bb_photo_thick_integrand(i_freq,i_tau) = 4.0_dp*pi*R_star2*two_pi_over_c_square*frequency(i_freq)* &
                                                     frequency(i_freq)*exp(-(tau(i_tau)*exponent_HeI(i_freq)))/ &
                                                     (exp(frequency(i_freq)*h_over_kT)-1.0)
            bb_photo_thin_integrand(i_freq,i_tau) = 4.0_dp*pi*R_star2*two_pi_over_c_square*frequency(i_freq)* &
                                                    frequency(i_freq)* exponent_HeI(i_freq)*exp(-(tau(i_tau)* &
                                                    exponent_HeI(i_freq)))/(exp(frequency(i_freq)*h_over_kT)-1.0)
            pl_photo_thick_integrand(i_freq,i_tau) = pl_scaling*frequency(i_freq)**(-pl_index_nominal)* &
                                                     exp(-tau(i_tau)*exponent_HeI(i_freq))
            pl_photo_thin_integrand(i_freq,i_tau) = pl_scaling*frequency(i_freq)**(-pl_index_nominal)*exponent_HeI(i_freq)* &
                                                    exp(-tau(i_tau)*exponent_HeI(i_freq))
          else
            bb_photo_thick_integrand(i_freq,i_tau) = 0.0  
            bb_photo_thin_integrand(i_freq,i_tau) = 0.0   
            pl_photo_thick_integrand(i_freq,i_tau) = 0.0  
            pl_photo_thin_integrand(i_freq,i_tau) = 0.0  
          endif            

          ! Assign values to the heating integrands
            bb_heat_thick_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                       bb_photo_thick_integrand(i_freq,i_tau)
            bb_heat_thick_integrand_HeI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                        bb_photo_thick_integrand(i_freq,i_tau)
            bb_heat_thin_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)*  &
                                                      bb_photo_thin_integrand(i_freq,i_tau)
            bb_heat_thin_integrand_HeI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                       bb_photo_thin_integrand(i_freq,i_tau)
            pl_heat_thick_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                       pl_photo_thick_integrand(i_freq,i_tau)
            pl_heat_thick_integrand_HeI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                        pl_photo_thick_integrand(i_freq,i_tau)
            pl_heat_thin_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                      pl_photo_thin_integrand(i_freq,i_tau)
            pl_heat_thin_integrand_HeI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                       pl_photo_thin_integrand(i_freq,i_tau)

        enddo 
        
      enddo   

      ! Make photo tables
      call vector_romberg (bb_photo_thick_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      bb_photo_thick_table(:,i_subband) = answer
      call vector_romberg (bb_photo_thin_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      bb_photo_thin_table(:,i_subband) = answer
      call vector_romberg (pl_photo_thick_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      pl_photo_thick_table(:,i_subband) = answer
      call vector_romberg (pl_photo_thin_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      pl_photo_thin_table(:,i_subband) = answer 
 
      ! Make heating tables
        call vector_romberg (bb_heat_thick_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thick_table(:,i_subband*2-2) = answer
        call vector_romberg (bb_heat_thick_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thick_table(:,i_subband*2-1) = answer
        call vector_romberg (bb_heat_thin_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thin_table(:,i_subband*2-2) = answer
        call vector_romberg (bb_heat_thin_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thin_table(:,i_subband*2-1) = answer
        call vector_romberg (pl_heat_thick_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thick_table(:,i_subband*2-2) = answer
        call vector_romberg (pl_heat_thick_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thick_table(:,i_subband*2-1) = answer
        call vector_romberg (pl_heat_thin_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thin_table(:,i_subband*2-2) = answer
        call vector_romberg (pl_heat_thin_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thin_table(:,i_subband*2-1) = answer

    enddo            

    ! In frequency band 3, fill in integrands and make tables

    ! Go through all the sub-bin in band 3
    do i_subband=NumBndin1+NumBndin2+1,NumBndin1+NumBndin2+NumBndin3

      ! Assign values to exponent array
      do i_freq=0,NumFreq

        frequency(i_freq) = freq_min(i_subband)+delta_freq(i_subband)*real(i_freq)

        exponent_HeII(i_freq) = ((frequency(i_freq)/freq_min(i_subband))**(-pl_index_cross_section_HeII(i_subband)))

      enddo

      ! Loop through the tau partition
      do i_tau=0,NumTau 
      
        ! Loop through the frequency partition 
        do i_freq=0,NumFreq
           weight(i_freq,i_tau) = delta_freq(i_subband)  

          ! Assign values to the photo integrands
          if (tau(i_tau)*exponent_HeII(i_freq) < 700.0) then  
            bb_photo_thick_integrand(i_freq,i_tau) = 4.0_dp*pi*R_star2*two_pi_over_c_square*frequency(i_freq)* &
                                                     frequency(i_freq)*exp(-tau(i_tau)*exponent_HeII(i_freq))/ &
                                                     (exp(frequency(i_freq)*h_over_kT)-1.0)  
            bb_photo_thin_integrand(i_freq,i_tau) = 4.0_dp*pi*R_star2*two_pi_over_c_square*frequency(i_freq)* &
                                                    frequency(i_freq)*exponent_HeII(i_freq)*exp(-tau(i_tau)* &
                                                    exponent_HeII(i_freq))/(exp(frequency(i_freq)*h_over_kT)-1.0)  
            pl_photo_thick_integrand(i_freq,i_tau) = pl_scaling*frequency(i_freq)**(-pl_index_nominal)* &
                                                     exp(-tau(i_tau)*exponent_HeII(i_freq))
            pl_photo_thin_integrand(i_freq,i_tau) = pl_scaling*frequency(i_freq)**(-pl_index_nominal)*exponent_HeII(i_freq)* &
                                                    exp(-tau(i_tau)*exponent_HeII(i_freq))
	  else
            bb_photo_thick_integrand(i_freq,i_tau) = 0.0
            bb_photo_thin_integrand(i_freq,i_tau) = 0.0
            pl_photo_thick_integrand(i_freq,i_tau) = 0.0
            pl_photo_thin_integrand(i_freq,i_tau) = 0.0
	  endif 

          ! Assign values to the heating integrands
            bb_heat_thick_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                       bb_photo_thick_integrand(i_freq,i_tau)
            bb_heat_thick_integrand_HeI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                        bb_photo_thick_integrand(i_freq,i_tau)
            bb_heat_thick_integrand_HeII(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeII)* &
                                                         bb_photo_thick_integrand(i_freq,i_tau)
            bb_heat_thin_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                      bb_photo_thin_integrand(i_freq,i_tau)
            bb_heat_thin_integrand_HeI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                       bb_photo_thin_integrand(i_freq,i_tau)
            bb_heat_thin_integrand_HeII(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeII)* &
                                                        bb_photo_thin_integrand(i_freq,i_tau)
            pl_heat_thick_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)*  &
                                                       pl_photo_thick_integrand(i_freq,i_tau)
            pl_heat_thick_integrand_HeI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                        pl_photo_thick_integrand(i_freq,i_tau)
            pl_heat_thick_integrand_HeII(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeII)* &
                                                         pl_photo_thick_integrand(i_freq,i_tau)
            pl_heat_thin_integrand_HI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
                                                      pl_photo_thin_integrand(i_freq,i_tau)
            pl_heat_thin_integrand_HeI(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeI)* &
                                                       pl_photo_thin_integrand(i_freq,i_tau)
            pl_heat_thin_integrand_HeII(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HeII)* &
                                                        pl_photo_thin_integrand(i_freq,i_tau)

        enddo

      enddo  

      ! Make photo tables
      call vector_romberg (bb_photo_thick_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      bb_photo_thick_table(:,i_subband) = answer
      call vector_romberg (bb_photo_thin_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      bb_photo_thin_table(:,i_subband) = answer
      call vector_romberg (pl_photo_thick_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      pl_photo_thick_table(:,i_subband) = answer
      call vector_romberg (pl_photo_thin_integrand,weight,NumFreq,NumFreq,NumTau,answer)
      pl_photo_thin_table(:,i_subband) = answer  

      ! Make heating tables
        call vector_romberg (bb_heat_thick_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thick_table(:,i_subband*3-NumBndin2-4) = answer
        call vector_romberg (bb_heat_thick_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thick_table(:,i_subband*3-NumBndin2-3) = answer
        call vector_romberg (bb_heat_thick_integrand_HeII,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thick_table(:,i_subband*3-NumBndin2-2) = answer
        call vector_romberg (bb_heat_thin_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thin_table(:,i_subband*3-NumBndin2-4) = answer
        call vector_romberg (bb_heat_thin_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thin_table(:,i_subband*3-NumBndin2-3) = answer
        call vector_romberg (bb_heat_thin_integrand_HeII,weight,NumFreq,NumFreq,NumTau,answer)
        bb_heat_thin_table(:,i_subband*3-NumBndin2-2) = answer
        call vector_romberg (pl_heat_thick_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thick_table(:,i_subband*3-NumBndin2-4) = answer
        call vector_romberg (pl_heat_thick_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thick_table(:,i_subband*3-NumBndin2-3) = answer
        call vector_romberg (pl_heat_thick_integrand_HeII,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thick_table(:,i_subband*3-NumBndin2-2) = answer
        call vector_romberg (pl_heat_thin_integrand_HI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thin_table(:,i_subband*3-NumBndin2-4) = answer
        call vector_romberg (pl_heat_thin_integrand_HeI,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thin_table(:,i_subband*3-NumBndin2-3) = answer
        call vector_romberg (pl_heat_thin_integrand_HeII,weight,NumFreq,NumFreq,NumTau,answer)
        pl_heat_thin_table(:,i_subband*3-NumBndin2-2) = answer

    enddo           

    ! deallocate the useless photo integrand
    deallocate(bb_photo_thick_integrand)
    deallocate(bb_photo_thin_integrand)
    deallocate(pl_photo_thick_integrand)
    deallocate(pl_photo_thin_integrand)

    ! deallocate the useless heating integrand
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

    if (rank .eq. 0) then
      write(logf,*) "End of photo table initialization" 
      write(logf,*)      
      flush(logf) 
    endif

  end subroutine photo_table_initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module photo_table
