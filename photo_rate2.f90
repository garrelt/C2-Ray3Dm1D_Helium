module photo_rate
  
  use precision, only: dp
  use file_admin, only: logf
  use cgsconstants, only: hplanck                     ! Planck constant
  use cgsphotoconstants, only: ion_freq_HI,&             ! HI ionization energy in frequency
                               ion_freq_HeI,&            ! HeI ionization energy in frequency
                               ion_freq_HeII           ! HeII ionization energy in frequency
  use parameter, only: NumTau, NumBndin1, NumBndin2, NumBndin3, NumFreqBnd, NumheatBin,&
                       minlogtau, dlogtau, CR1, CR2, bR1, dR1, aR2, bR2, y1R, y2R, xeb, &
                       do_secondary
  use array, only: sigma_HI, sigma_HeI, sigma_HeII, &
                   f1ion_HI, f1ion_HeI, f1ion_HeII, &
                   f2ion_HI, f2ion_HeI, f2ion_HeII, &
                   f2heat_HI, f2heat_HeI, f2heat_HeII, &
                   f1heat_HI, f1heat_HeI, f1heat_HeII, &
                   bb_photo_thick_table, bb_photo_thin_table, &
                   pl_photo_thick_table, pl_photo_thin_table, &
                   bb_heat_thick_table, bb_heat_thin_table, &
                   pl_heat_thick_table, pl_heat_thin_table
  use type_definition, only: photrates, tablepos
    use sourceprops, only: BB_NormFlux,PL_NormFlux
    use cgsphotoconstants

  implicit none

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! this subroutine calculates photo-ionization rates at a particular sets of column density
  subroutine photoion_rates (phi,colum_in_HI,colum_out_HI,colum_in_HeI,colum_out_HeI, &
	       	       colum_in_HeII,colum_out_HeII,vol,nsrc,i_state)

    implicit none

    type(photrates),intent(out) :: phi

    ! Incoming and outgoing HI column density
    real(kind=dp), intent(in) :: colum_in_HI, colum_out_HI

    ! Incoming and outgoing HeI column density
    real(kind=dp), intent(in) :: colum_in_HeI, colum_out_HeI

    ! Incoming and outgoing HeII column density
    real(kind=dp), intent(in) :: colum_in_HeII, colum_out_HeII

    ! Volume of shell cell
    real(kind=dp), intent(in) :: vol

    real(kind=dp), intent(in) :: i_state

    ! Number of the source
    integer, intent(in) :: nsrc 

    character :: sourcetype

    integer :: i_tau, i_subband
    real(kind=dp) :: colum_cell_HI
    real(kind=dp) :: colum_cell_HeI
    real(kind=dp) :: colum_cell_HeII
    real(kind=dp) :: NFlux
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_in_all
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_out_all
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_cell_HI
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_cell_HeI
    real(kind=dp), dimension(1:NumFreqBnd) :: tau_cell_HeII
    type(tablepos) :: tau_pos_in, tau_pos_out


    ! The column densities (HI, HeI, HeII) at current cell
    colum_cell_HI = colum_out_HI-colum_in_HI
    colum_cell_HeI = colum_out_HeI-colum_in_HeI
    colum_cell_HeII = colum_out_HeII-colum_in_HeII 

    ! The optical depths (HI, HeI, HeII) at current cell
    do i_subband=1,NumFreqBnd
      tau_cell_HI(i_subband) = colum_cell_HI*sigma_HI(i_subband)
      tau_cell_HeI(i_subband) = colum_cell_HeI*sigma_HeI(i_subband)
      tau_cell_HeII(i_subband) = colum_cell_HeII*sigma_HeII(i_subband)
    enddo      
              
    ! total tau_in (including HI, HeI, HeII)
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

    ! find the table positions for the optical depth (ingoing)
    do i_subband=1,NumFreqBnd  
      tau_pos_in%tau(i_subband) = log10(max(1.0e-20_dp,tau_in_all(i_subband)))
      tau_pos_in%odpos(i_subband) = min(real(NumTau,dp),max(0.0_dp,1.0+ &
                                      (tau_pos_in%tau(i_subband)-minlogtau)/dlogtau))
      tau_pos_in%ipos(i_subband) = int(tau_pos_in%odpos(i_subband))
      tau_pos_in%residual(i_subband) = tau_pos_in%odpos(i_subband)-real(tau_pos_in%ipos(i_subband),dp)
      tau_pos_in%ipos_p1(i_subband) = min(NumTau,tau_pos_in%ipos(i_subband)+1)
    enddo

    ! find the table positions for the optical depth (outgoing)
    do i_subband=1,NumFreqBnd  
      tau_pos_out%tau(i_subband) = log10(max(1.0e-20_dp,tau_out_all(i_subband)))
      tau_pos_out%odpos(i_subband) = min(real(NumTau,dp),max(0.0_dp,1.0+ &
                                     (tau_pos_out%tau(i_subband)-minlogtau)/dlogtau))
      tau_pos_out%ipos(i_subband) = int(tau_pos_out%odpos(i_subband))
      tau_pos_out%residual(i_subband) = tau_pos_out%odpos(i_subband)-real(tau_pos_out%ipos(i_subband),dp)
      tau_pos_out%ipos_p1(i_subband) = min(NumTau,tau_pos_out%ipos(i_subband)+1)
    enddo 

    if (nsrc.eq.0) then
      NFlux=BB_NormFlux(0)
      sourcetype='B'
    elseif (PL_NormFlux(nsrc) .le. 1.0e-100_dp) then  
      NFlux=BB_NormFlux(nsrc) 
      sourcetype='B'
    else
      NFlux=PL_NormFlux(nsrc)
      sourcetype='P'
    endif
!$OMP CRITICAL
    call lookuptable(tau_pos_in,tau_pos_out,phi,tau_in_all,tau_out_all, &
                     tau_cell_HI,tau_cell_HeI,tau_cell_HeII,NFlux,sourcetype, &
                     vol,i_state,colum_cell_HI,colum_cell_HeI, colum_cell_HeII )
!$OMP END CRITICAL
  end subroutine photoion_rates
 
  ! find out the correct position in the photo and heating tables
  subroutine lookuptable(tau_pos_in,tau_pos_out,phi,tau_in_all,tau_out_all, &
                         tau_cell_HI,tau_cell_HeI,tau_cell_HeII,NFlux,sourcetype, &
                         vol,i_state,colum_cell_HI,colum_cell_HeI,colum_cell_HeII)

    implicit none

    type(photrates), intent(out) :: phi
    type(tablepos), intent(in) :: tau_pos_in, tau_pos_out
    real(kind=dp), intent(in) :: NFlux, vol, i_state
    real(kind=dp), intent(in) :: colum_cell_HI, colum_cell_HeI, colum_cell_HeII
    real(kind=dp), dimension(NumFreqBnd), intent(in) :: tau_in_all, tau_out_all
    character,intent(in) :: sourcetype
    real(kind=dp), dimension(NumFreqBnd),intent(in) :: tau_cell_HI, tau_cell_HeI, tau_cell_HeII

    integer ::  n, i_subband, i
    real(kind=dp) :: phi_heat_HI, phi_heat_HeI, phi_heat_HeII
    real(kind=dp) :: f_heat, f_ion_HI, f_ion_HeI
    real(kind=dp) :: phi_photo_in_all, phi_photo_out_all, phi_photo_all
    real(kind=dp) :: fra_sum1, fra_sum2, fra_sum3, fra_sum4
    real(kind=dp), pointer, dimension(:,:) :: photo_thick_table, photo_thin_table
    real(kind=dp), pointer, dimension(:,:) :: heat_thick_table, heat_thin_table
    real(kind=dp), dimension(NumBndin1+1:NumBndin1+NumBndin2+NumBndin3) :: scaling_HI
    real(kind=dp), dimension(NumBndin1+1:NumBndin1+NumBndin2+NumBndin3) :: scaling_HeI
    real(kind=dp), dimension(NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+NumBndin3) :: scaling_HeII
    real(kind=dp), dimension(1:NumheatBin) :: scaling
    real(kind=dp) :: test1, test2
    real(kind=dp) :: tau_photo_limit = 1.0e-7 
    real(kind=dp) :: tau_heat_limit = 1.0e-4

    ! pointers point to some variables
    if (sourcetype.eq.'B') then 
      photo_thick_table => bb_photo_thick_table
      photo_thin_table => bb_photo_thin_table
      heat_thick_table => bb_heat_thick_table
      heat_thin_table => bb_heat_thin_table
    elseif (sourcetype.eq.'P') then
      photo_thick_table => pl_photo_thick_table
      photo_thin_table => pl_photo_thin_table
      heat_thick_table => pl_heat_thick_table
      heat_thin_table => pl_heat_thin_table
    endif

    ! initialization
    phi%photo_cell_HI = 0.0_dp
    phi%photo_cell_HeI = 0.0_dp
    phi%photo_cell_HeII = 0.0_dp
    phi_heat_HI = 0.0_dp
    phi_heat_HeI = 0.0_dp
    phi_heat_HeII = 0.0_dp
    phi%heat_cell_HI = 0.0_dp
    phi%heat_cell_HeI = 0.0_dp
    phi%heat_cell_HeII = 0.0_dp
    phi%photo_in = 0.0_dp
    f_heat = 0.0_dp
    f_ion_HI = 0.0_dp
    f_ion_HeI = 0.0_dp

    ! loop through all frequency band
    do i_subband=1,NumFreqBnd

      ! Incoming, outcoming, current cell total photoionization rate
      phi_photo_in_all = (photo_thick_table(tau_pos_in%ipos(i_subband),i_subband)+ &
                         (photo_thick_table(tau_pos_in%ipos_p1(i_subband),i_subband)- &
                         photo_thick_table(tau_pos_in%ipos(i_subband),i_subband))* &
                         tau_pos_in%residual(i_subband))*NFlux
      phi%photo_in = phi%photo_in+phi_photo_in_all

      ! When current cell is optically thick
      if (abs(tau_out_all(i_subband)-tau_in_all(i_subband)) .gt. tau_photo_limit) then
        phi_photo_out_all = (photo_thick_table(tau_pos_out%ipos(i_subband),i_subband)+ &
                            (photo_thick_table(tau_pos_out%ipos_p1(i_subband),i_subband)- &
                            photo_thick_table(tau_pos_out%ipos(i_subband),i_subband))* &
                            tau_pos_out%residual(i_subband))*NFlux 
        phi_photo_all = phi_photo_in_all-phi_photo_out_all

      ! When current cell is optically thin
      else
        phi_photo_all = ((photo_thin_table(tau_pos_in%ipos(i_subband),i_subband)+ &
                        (photo_thin_table(tau_pos_in%ipos_p1(i_subband),i_subband)- &
                        photo_thin_table(tau_pos_in%ipos_p1(i_subband),i_subband))* &
                        tau_pos_in%residual(i_subband))* &
                        (tau_out_all(i_subband)-tau_in_all(i_subband)))*NFlux
        phi_photo_out_all = phi_photo_in_all-phi_photo_all
      endif

      ! Current cell individual photoionization rate of HI, HeI, HeII
      select case (i_subband) 

      ! band 1
      case (NumBndin1) 
    
        phi%photo_cell_HI = phi_photo_all/vol

      ! band 2
      case (NumBndin1+1:NumBndin1+NumBndin2)
        
        call scale_int2(scaling_HI(i_subband),scaling_HeI(i_subband),colum_cell_HI,colum_cell_HeI, i_subband)

        phi%photo_cell_HI = phi%photo_cell_HI+scaling_HI(i_subband)*phi_photo_all/vol 
        phi%photo_cell_HeI = phi%photo_cell_HeI+scaling_HeI(i_subband)*phi_photo_all/vol

      ! band 3
      case (NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+NumBndin3)

        call scale_int3(scaling_HI(i_subband),scaling_HeI(i_subband),scaling_HeII(i_subband), &
                        colum_cell_HI,colum_cell_HeI,colum_cell_HeII,i_subband)
        phi%photo_cell_HI = phi%photo_cell_HI+scaling_HI(i_subband)*phi_photo_all/vol
        phi%photo_cell_HeI = phi%photo_cell_HeI+scaling_HeI(i_subband)*phi_photo_all/vol
        phi%photo_cell_HeII = phi%photo_cell_HeII+scaling_HeII(i_subband)*phi_photo_all/vol

      end select
     
    enddo

      ! in general, I'm following Ricotti et al 2002
      CR1 = (/0.3908_dp, 0.0554_dp, 1.0_dp/)
      bR1 = (/0.4092_dp, 0.4614_dp, 0.2663_dp/)
      dR1 = (/1.7592_dp, 1.6660_dp, 1.3163_dp/)
      CR2 = (/0.6941_dp,0.0984_dp,3.9811_dp/)
      aR2 = (/0.2_dp,0.2_dp,0.4_dp/)
      bR2 = (/0.38_dp,0.38_dp,0.34_dp/)
      test1 = 0.0_dp
      test2 = 0.0_dp

      do i=1,3
        y1R(i) = CR1(i)*(1.0_dp-i_state**bR1(i))**dR1(i)
        xeb = 1.0_dp-i_state**bR2(i) 
        y2R(i) = CR2(i)*i_state**aR2(i)*xeb*xeb
      enddo

      ! Current cell individual heating rates of HI, HeI, HeII
      do i_subband=1,NumFreqBnd 
        phi_heat_HI = 0.0_dp
        phi_heat_HeI = 0.0_dp
        phi_heat_HeII = 0.0_dp
      
        select case (i_subband)

        ! Incoming, outcoming, current cell HI heating rate at band 1
        case (NumBndin1)

          phi%heat_in_HI = (heat_thick_table(tau_pos_in%ipos(1),1)+(heat_thick_table(tau_pos_in%ipos_p1(1),1)- &
                           heat_thick_table(tau_pos_in%ipos(1),1))*tau_pos_in%residual(1))*NFlux 

          ! When current cell is HI optically thick     
          if (abs(tau_cell_HI(i_subband)) .gt. tau_heat_limit) then
            phi%heat_out_HI = (heat_thick_table(tau_pos_out%ipos(1),1)+(heat_thick_table(tau_pos_out%ipos_p1(1),1)- &
                              heat_thick_table(tau_pos_out%ipos(1),1))*tau_pos_out%residual(1))*NFlux 
            phi_heat_HI = (phi%heat_in_HI-phi%heat_out_HI)/vol
 
          ! When current cell is HI optically thin
          else
            phi_heat_HI = (heat_thin_table(tau_pos_in%ipos(1),1)+(heat_thin_table(tau_pos_in%ipos_p1(1),1)- &
                          heat_thin_table(tau_pos_in%ipos(1),1))*tau_pos_in%residual(1))* &   
                          (tau_out_all(1)-tau_in_all(1))*NFlux 
            phi%heat_out_HI = phi%heat_in_HI+phi_heat_HI
            phi_heat_HI = phi_heat_HI/vol
          endif

          f_heat = phi_heat_HI

        ! Incoming, outcoming, current cell HI, HeI heating rate at band 2      
        case (NumBndin1+1:NumBndin1+NumBndin2)
       
          phi%heat_in_HI = (heat_thick_table(tau_pos_in%ipos(i_subband),2*i_subband-2)+ &
                           (heat_thick_table(tau_pos_in%ipos_p1(i_subband),2*i_subband-2)- &
                           heat_thick_table(tau_pos_in%ipos(i_subband),2*i_subband-2))* &
                           tau_pos_in%residual(i_subband))*NFlux

          phi%heat_in_HeI = (heat_thick_table(tau_pos_in%ipos(i_subband),2*i_subband-1)+ &
                            (heat_thick_table(tau_pos_in%ipos_p1(i_subband),2*i_subband-1)- &
                            heat_thick_table(tau_pos_in%ipos(i_subband),2*i_subband-1))* &
                            tau_pos_in%residual(i_subband))*NFlux

          ! When current cell is HI optically thick  
          if (abs(tau_cell_HI(i_subband)) .gt. tau_heat_limit) then
            phi%heat_out_HI = (heat_thick_table(tau_pos_out%ipos(i_subband),2*i_subband-2)+ &
                              (heat_thick_table(tau_pos_out%ipos_p1(i_subband),2*i_subband-2)- &
                              heat_thick_table(tau_pos_out%ipos(i_subband),2*i_subband-2))* &
                              tau_pos_out%residual(i_subband))*NFlux
            phi_heat_HI = scaling_HI(i_subband)*(phi%heat_in_HI-phi%heat_out_HI)/vol 

          ! When current cell is HI optically thin
          else
            phi_heat_HI = scaling_HI(i_subband)*(((heat_thin_table(tau_pos_in%ipos(i_subband),2*i_subband-2)+ &
                          (heat_thin_table(tau_pos_in%ipos_p1(i_subband),2*i_subband-2)- &
                          heat_thin_table(tau_pos_in%ipos_p1(i_subband),2*i_subband-2))* &
                          tau_pos_in%residual(i_subband))*(tau_out_all(i_subband)-tau_in_all(i_subband))))*NFlux
            phi%heat_out_HI = phi%heat_in_HI+phi_heat_HI
            phi_heat_HI = phi_heat_HI/vol
          endif

          ! When current cell is HeI optically thick      
          if (abs(tau_cell_HeI(i_subband)) .gt. tau_heat_limit) then
            phi%heat_out_HeI = (heat_thick_table(tau_pos_out%ipos(i_subband),2*i_subband-1)+ &
                               (heat_thick_table(tau_pos_out%ipos_p1(i_subband),2*i_subband-1)- &
                               heat_thick_table(tau_pos_out%ipos(i_subband),2*i_subband-1))* &
                               tau_pos_out%residual(i_subband))*NFlux
            phi_heat_HeI = scaling_HeI(i_subband)*(phi%heat_in_HeI-phi%heat_out_HeI)/vol

          ! When current cell is HeI optically thin
          else 
            phi_heat_HeI = scaling_HeI(i_subband)*((heat_thin_table(tau_pos_in%ipos(i_subband),2*i_subband-1)+&
                           (heat_thin_table(tau_pos_in%ipos_p1(i_subband),2*i_subband-1)- &
                           heat_thin_table(tau_pos_in%ipos_p1(i_subband),2*i_subband-1))* &
                           tau_pos_in%residual(i_subband))*(tau_out_all(i_subband)-tau_in_all(i_subband)))*NFlux
            phi%heat_out_HeI=phi%heat_in_HeI+phi_heat_HeI
            phi_heat_HeI = phi_heat_HeI/vol
          endif

          fra_sum1 = f1ion_HI(i_subband)*phi_heat_HI+f1ion_HeI(i_subband)*phi_heat_HeI
          fra_sum2 = f2ion_HI(i_subband)*phi_heat_HI+f2ion_HeI(i_subband)*phi_heat_HeI
          fra_sum3 = f1heat_HI(i_subband)*phi_heat_HI+f1heat_HeI(i_subband)*phi_heat_HeI
          fra_sum4 = f2heat_HI(i_subband)*phi_heat_HI+f2heat_HeI(i_subband)*phi_heat_HeI

          if (do_secondary) then
            ! These are all cumulative
            f_ion_HeI = f_ion_HeI+y1R(2)*fra_sum1-y2R(2)*fra_sum2  
            f_ion_HI = f_ion_HI+y1R(1)*fra_sum1-y2R(1)*fra_sum2
            f_heat = f_heat+phi_heat_HI+phi_heat_HeI-y1R(3)*fra_sum3+y2R(3)*fra_sum4 
          else
            f_heat = f_heat+phi_heat_HI+phi_heat_HeI
          endif

        ! Incoming, outcoming, current cell HI, HeI, HeII heating rate at band 3   
        case (NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+NumBndin3)
          phi%heat_in_HI = (heat_thick_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-4)+ &
                           (heat_thick_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-4)- &
                           heat_thick_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-4))* &
                           tau_pos_in%residual(i_subband))*NFlux 
          phi%heat_in_HeI = (heat_thick_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-3)+ &
                            (heat_thick_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-3)- &
                            heat_thick_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-3))* &
                            tau_pos_in%residual(i_subband))*NFlux
          phi%heat_in_HeII = (heat_thick_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-2)+ &
                             (heat_thick_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-2)- &
                             heat_thick_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-2))* &
                             tau_pos_in%residual(i_subband))*NFlux

          ! When current cell is HI optically thick 
          if (abs(tau_cell_HI(i_subband)) .gt. tau_heat_limit) then
            phi%heat_out_HI = (heat_thick_table(tau_pos_out%ipos(i_subband),3*i_subband-NumBndin2-4)+ &
                              (heat_thick_table(tau_pos_out%ipos_p1(i_subband),3*i_subband-NumBndin2-4)- &
                              heat_thick_table(tau_pos_out%ipos(i_subband),3*i_subband-NumBndin2-4))* &
                              tau_pos_out%residual(i_subband))*NFlux
            phi_heat_HI = scaling_HI(i_subband)*(phi%heat_in_HI-phi%heat_out_HI)/vol

          ! When current cell is HI optically thin
          else
            phi_heat_HI = scaling_HI(i_subband)*((heat_thin_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-4)+ &
                          (heat_thin_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-4)- &
                          heat_thin_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-4))* &
                          tau_pos_in%residual(i_subband))*(tau_out_all(i_subband)-tau_in_all(i_subband)))*NFlux
            phi%heat_out_HI = phi%heat_in_HI+phi_heat_HI
            phi_heat_HI = phi_heat_HI/vol
          endif

          ! When current cell is HeI optically thick   
          if (abs(tau_cell_HeI(i_subband)) .gt. tau_heat_limit) then 
            phi%heat_out_HeI = (heat_thick_table(tau_pos_out%ipos(i_subband),3*i_subband-NumBndin2-3)+ &
                               (heat_thick_table(tau_pos_out%ipos_p1(i_subband),3*i_subband-NumBndin2-3)- &
                               heat_thick_table(tau_pos_out%ipos(i_subband),3*i_subband-NumBndin2-3))* &
                               tau_pos_out%residual(i_subband))*NFlux
            phi_heat_HeI = scaling_HeI(i_subband)*(phi%heat_in_HeI-phi%heat_out_HeI)/vol

          ! When current cell is HeI optically thin
          else
            phi_heat_HeI = scaling_HeI(i_subband)*((heat_thin_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-3)+ &
                           (heat_thin_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-3)- &
                           heat_thin_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-3))* &
                           tau_pos_in%residual(i_subband))*(tau_out_all(i_subband)-tau_in_all(i_subband)))*NFlux
            phi%heat_out_HeI = phi%heat_in_HeI+phi_heat_HeI
            phi_heat_HeI = phi_heat_HeI/vol
          endif

          ! When current cell is HeII optically thick      
          if (abs(tau_cell_HeII(i_subband)) .gt. tau_heat_limit) then
            phi%heat_out_HeII = (heat_thick_table(tau_pos_out%ipos(i_subband),3*i_subband-NumBndin2-2)+ &
                                (heat_thick_table(tau_pos_out%ipos_p1(i_subband),3*i_subband-NumBndin2-2)- &
                                heat_thick_table(tau_pos_out%ipos(i_subband),3*i_subband-NumBndin2-2))* &
                                tau_pos_out%residual(i_subband))*NFlux
            phi_heat_HeII = scaling_HeII(i_subband)*(phi%heat_in_HeII-phi%heat_out_HeII)/vol

          ! When current cell is HeII optically thin
          else 
            phi_heat_HeII = scaling_HeII(i_subband)*((heat_thin_table(tau_pos_in%ipos(i_subband),3*i_subband-NumBndin2-2)+ &
                            (heat_thin_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-2)- &
                            heat_thin_table(tau_pos_in%ipos_p1(i_subband),3*i_subband-NumBndin2-2))* &
                            tau_pos_in%residual(i_subband))*(tau_out_all(i_subband)- tau_in_all(i_subband)))*NFlux
            phi%heat_out_HeII = phi%heat_in_HeII+phi_heat_HeII
            phi_heat_HeII = phi_heat_HeII/vol
          endif
          
          fra_sum1 = f1ion_HI(i_subband)*phi_heat_HI+f1ion_HeI(i_subband)*phi_heat_HeI+f1ion_HeII(i_subband)*phi_heat_HeII
          fra_sum2 = f2ion_HI(i_subband)*phi_heat_HI+f2ion_HeI(i_subband)*phi_heat_HeI+f2ion_HeII(i_subband)*phi_heat_HeII
          fra_sum3 = f1heat_HI(i_subband)*phi_heat_HI+f1heat_HeI(i_subband)*phi_heat_HeI+f1heat_HeII(i_subband)*phi_heat_HeII
          fra_sum4 = f2heat_HI(i_subband)*phi_heat_HI+f2heat_HeI(i_subband)*phi_heat_HeI+f2heat_HeII(i_subband)*phi_heat_HeII
          
          if (do_secondary) then
            ! These are all cumulative
            f_ion_HeI = f_ion_HeI+y1R(2)*fra_sum1-y2R(2)*fra_sum2
            f_ion_HI = f_ion_HI+y1R(1)*fra_sum1-y2R(1)*fra_sum2
            f_heat = f_heat+phi_heat_HI+phi_heat_HeI+phi_heat_HeII-y1R(3)*fra_sum3+y2R(3)*fra_sum4   
          else
            f_heat = f_heat+phi_heat_HI+phi_heat_HeI+phi_heat_HeII 
          endif

        end select
	
      enddo 

      !Total heating rate on current cell
      phi%heat = f_heat 
      !Final HI photoionization rate modified by secondary ionization
      phi%photo_cell_HI = phi%photo_cell_HI+f_ion_HI/(ion_freq_HI*hplanck) 
      !Final HeI photoionization rate modified by secondary ionization
      phi%photo_cell_HeI = phi%photo_cell_HeI+f_ion_HeI/(ion_freq_HeI*hplanck)  

  end subroutine lookuptable

  ! give scalings of species for division of photoionization and heating to species
  subroutine scale_int2(scaling_HI,scaling_HeI,colum_cell_HI,colum_cell_HeI,i_subband)

    real(kind=dp),intent(in) :: colum_cell_HI, colum_cell_HeI
    integer,intent(in) :: i_subband
    real(kind=dp),intent(out):: scaling_HI, scaling_HeI
    real(kind=dp) :: forscaleing

    forscaleing = 1.0_dp/(sigma_HI(i_subband)*colum_cell_HI+sigma_HeI(i_subband)*colum_cell_HeI)
    scaling_HI = sigma_HI(i_subband)*colum_cell_HI*forscaleing
    scaling_HeI = sigma_HeI(i_subband)*colum_cell_HeI*forscaleing

  end subroutine scale_int2  

  ! give scalings of species for division of photoionization and heating to species
  subroutine scale_int3(scaling_HI, scaling_HeI, scaling_HeII, colum_cell_HI, colum_cell_HeI, colum_cell_HeII, i_subband)

    real(kind=dp),intent(in) :: colum_cell_HI ,colum_cell_HeI, colum_cell_HeII
    integer,intent(in) :: i_subband
    real(kind=dp),intent(out):: scaling_HI, scaling_HeI, scaling_HeII
    real(kind=dp) :: forscaleing

    forscaleing = 1.0_dp/(sigma_HI(i_subband)*colum_cell_HI+sigma_HeI(i_subband)*colum_cell_HeI+ &
                  sigma_HeII(i_subband)*colum_cell_HeII) 
    scaling_HI = colum_cell_HI*sigma_HI(i_subband) *forscaleing
    scaling_HeI = colum_cell_HeI*sigma_HeI(i_subband)*forscaleing
    scaling_HeII = colum_cell_HeII*sigma_HeII(i_subband)*forscaleing 

  end subroutine scale_int3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module photo_rate
