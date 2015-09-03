!     This module contains data and routines which deal with radiative
!     effects. Its main part deal with photo-ionizing radiation, but it
!     also initializes other radiative properties, such as cooling (which
!     are contained in different modules).
!     It can be used in hydrodynamic or stand-alone radiative transfer 
!     calculations.

module radiation_photoionrates
  
  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use cgsconstants, only: hplanck                     ! Planck constant
  use cgsphotoconstants, only: ion_freq_HI,&             ! HI ionization energy in frequency
                               ion_freq_HeI,&            ! HeI ionization energy in frequency
                               ion_freq_HeII             ! HeII ionization energy in frequency

  use material, only: isothermal

  use radiation_sizes, only: NumFreqBnd,NumBndin1,NumBndin2,NumBndin3
  use radiation_sizes, only: sigma_HI, sigma_HeI, sigma_HeII
  use radiation_sizes, only: NumTau, NumheatBin
  use radiation_sizes, only: f1ion_HI, f1ion_HeI, f1ion_HeII
  use radiation_sizes, only: f2ion_HI, f2ion_HeI, f2ion_HeII
  use radiation_sizes, only: f1heat_HI, f1heat_HeI, f1heat_HeII
  use radiation_sizes, only: f2heat_HI, f2heat_HeI, f2heat_HeII
  
  use radiation_tables, only: minlogtau, dlogtau
  use radiation_tables, only: bb_photo_thick_table, bb_photo_thin_table 
  use radiation_tables, only: bb_heat_thick_table, bb_heat_thin_table 
  use radiation_tables, only: bb_FreqBnd_UpperLimit, bb_FreqBnd_LowerLimit

#ifdef PL
  use radiation_tables, only: pl_photo_thick_table, pl_photo_thin_table
  use radiation_tables, only: pl_heat_thick_table, pl_heat_thin_table 
  use radiation_tables, only: pl_FreqBnd_UpperLimit, pl_FreqBnd_LowerLimit
#endif
#ifdef QUASARS
  use radiation_tables, only: qpl_photo_thick_table, qpl_photo_thin_table
  use radiation_tables, only: qpl_heat_thick_table, qpl_heat_thin_table
  use radiation_tables, only: qpl_FreqBnd_UpperLimit, qpl_FreqBnd_LowerLimit
#endif
  implicit none

  ! Variables connected to secondary ionizations.
  logical,parameter :: do_secondary_ionizations=.true.

  ! in general, I'm following Ricotti et al 2002
  real(kind=dp), dimension(1:3) :: CR1=(/0.3908_dp, 0.0554_dp, 1.0_dp/)
  real(kind=dp), dimension(1:3) :: bR1=(/0.4092_dp, 0.4614_dp, 0.2663_dp/)
  real(kind=dp), dimension(1:3) :: dR1=(/1.7592_dp, 1.6660_dp, 1.3163_dp/)
    
  real(kind=dp), dimension(1:3) :: CR2=(/0.6941_dp,0.0984_dp,3.9811_dp/)
  real(kind=dp), dimension(1:3) :: aR2=(/0.2_dp,0.2_dp,0.4_dp/)
  real(kind=dp), dimension(1:3) :: bR2=(/0.38_dp,0.38_dp,0.34_dp/)
  !real(kind=dp), dimension(1:3) :: dR2=(/2.0_dp,2.0_dp,2.0_dp/) write explicitly ^2 -> introduce xeb

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
     real(kind=dp) :: photo_out               ! Total photoionization rate incoming to the cell
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! this subroutine calculates photo-ionization rates at a particular sets of column density
  function photoion_rates (colum_in_HI,colum_out_HI, &
       colum_in_HeI,colum_out_HeI, &
       colum_in_HeII,colum_out_HeII, &
       vol,nsrc,i_state)

    use sourceprops, only: NormFlux
#ifdef PL
    use sourceprops, only: NormFluxPL
#endif
#ifdef QUASARS
    use sourceprops, only: NormFluxQPL
#endif
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
#ifdef PL
    if (NormFluxPL(nsrc) > 0.0) &  
         phi = phi + photo_lookuptable(tau_pos_in,tau_pos_out, &
         tau_in_all,tau_out_all, &
         NormFluxPL(nsrc),"P",vol, &
         scaling_HI,scaling_HeI,scaling_HeII)
#endif    
#ifdef QUASARS
    if (NormFluxQPL(nsrc) > 0.0) &
         phi = phi + photo_lookuptable(tau_pos_in,tau_pos_out, &
         tau_in_all,tau_out_all, &
         NormFluxQPL(nsrc),"Q",vol, &
         scaling_HI,scaling_HeI,scaling_HeII)
#endif
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
#ifdef PL
       if (NormFluxPL(nsrc) > 0.0) &  
            phi = phi + heat_lookuptable(tau_pos_in,tau_pos_out, &
            tau_in_all,tau_out_all, &
            tau_cell_HI,tau_cell_HeI,tau_cell_HeII,NormFluxPL(nsrc),"P", &
            vol,i_state, &
            scaling_HI,scaling_HeI,scaling_HeII)
#endif       
#ifdef QUASARS
       if (NormFluxQPL(nsrc) > 0.0) &
            phi = phi + heat_lookuptable(tau_pos_in,tau_pos_out, &
            tau_in_all,tau_out_all, &
            tau_cell_HI,tau_cell_HeI,tau_cell_HeII,NormFluxQPL(nsrc),"Q", &
            vol,i_state, &
            scaling_HI,scaling_HeI,scaling_HeII)
#endif
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
    real(kind=dp), pointer, dimension(:,:),intent(in) :: table
    !real(kind=dp),dimension(0:NumTau, 1:tablesize),intent(in) :: table
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
    
    integer ::  i_subband
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
#ifdef PL
    elseif (table_type == "P") then
       photo_thick_table => pl_photo_thick_table
       photo_thin_table => pl_photo_thin_table
       Minimum_FreqBnd=pl_FreqBnd_LowerLimit
       Maximum_FreqBnd=pl_FreqBnd_UpperLimit
#endif
#ifdef QUASARS
    elseif (table_type == "Q") then
       photo_thick_table => qpl_photo_thick_table
       photo_thin_table => qpl_photo_thin_table
       Minimum_FreqBnd=qpl_FreqBnd_LowerLimit
       Maximum_FreqBnd=qpl_FreqBnd_UpperLimit
#endif
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

       ! Collect all outgoing photons
       photo_lookuptable%photo_out = &
            photo_lookuptable%photo_out + phi_photo_out_all

       ! Current cell individual photoionization rate of HI, HeI, HeII
       select case (i_subband) 
          
          ! band 1
       case (1:NumBndin1) 
          
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
    real(kind=dp) :: df_heat, df_ion_HI, df_ion_HeI
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
#ifdef PL
    elseif (table_type == "P") then
       heat_thick_table => pl_heat_thick_table
       heat_thin_table => pl_heat_thin_table
       Minimum_FreqBnd=pl_FreqBnd_LowerLimit
       Maximum_FreqBnd=pl_FreqBnd_UpperLimit
#endif
#ifdef QUASARS
    elseif (table_type == "Q") then
       heat_thick_table => qpl_heat_thick_table
       heat_thin_table => qpl_heat_thin_table
       Minimum_FreqBnd=qpl_FreqBnd_LowerLimit
       Maximum_FreqBnd=qpl_FreqBnd_UpperLimit
#endif
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
    df_ion_HI = 0.0
    df_ion_HeI = 0.0
    df_heat=0.0

    ! Set parameters for secondary ionizations (following Ricotti et al. 2002)
    if (do_secondary_ionizations) then
       do i=1,3
          y1R(i)= CR1(i)*(1.0_dp-i_state**bR1(i))**dR1(i)
          xeb=1.0_dp-i_state**bR2(i) 
          y2R(i)= CR2(i)*i_state**aR2(i)*xeb*xeb
       enddo
    endif
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
       case (1:NumBndin1)
          
          ! Incoming current cell HI heating rate at band 1
          phi_heat_in_HI = NFlux * read_table(heat_thick_table,NumheatBin, &
               tau_pos_in,i_subband,i_subband)
          
          ! When current cell is HI optically thick (in total tau)
          if (abs(tau_out_all(i_subband)-tau_in_all(i_subband)) > &
               tau_heat_limit) then
             phi_heat_out_HI = NFlux * read_table(heat_thick_table,NumheatBin, &
                  tau_pos_out,i_subband,i_subband)
             phi_heat_HI = (phi_heat_in_HI-phi_heat_out_HI)/vol
             
             ! When current cell is HI optically thin
          else
             phi_heat_HI = NFlux * tau_cell_HI(i_subband) * &
                  read_table(heat_thin_table,NumheatBin, &
                  tau_pos_in,i_subband,i_subband)
             !phi_heat_out_HI = phi_heat_in_HI-phi_heat_HI
             phi_heat_HI = phi_heat_HI/vol
          endif
          
          ! Save the heating in f_heat variable
          df_heat = phi_heat_HI

          ! band 2
       case (NumBndin1+1:NumBndin1+NumBndin2)
          
          ! Band 2, HI
          phi_heat_in_HI = NFlux * read_table(heat_thick_table,NumheatBin, &
               tau_pos_in,i_subband,2*i_subband-NumBndin1-1)
          phi_heat_in_HeI = NFlux * read_table(heat_thick_table,NumheatBin, &
               tau_pos_in,i_subband,2*i_subband-NumBndin1)
          
          ! When current cell is optically thick  (in total tau)
          if (abs(tau_out_all(i_subband)-tau_in_all(i_subband)) > &
               tau_heat_limit) then
             phi_heat_out_HI = NFlux * &
                  read_table(heat_thick_table,NumheatBin, &
                  tau_pos_out,i_subband,2*i_subband-NumBndin1-1)
             phi_heat_HI = scaling_HI(i_subband)*(phi_heat_in_HI-phi_heat_out_HI)/vol 
             
             phi_heat_out_HeI = NFlux * &
                  read_table(heat_thick_table,NumheatBin, &
                  tau_pos_out,i_subband,2*i_subband-NumBndin1)
             phi_heat_HeI = scaling_HeI(i_subband)*(phi_heat_in_HeI-phi_heat_out_HeI)/vol

             ! When current cell optically thin
          else
             ! By multiplying here with tau_cell_HI we already
             ! apply the scaling_HI factor!
             phi_heat_HI = NFlux * tau_cell_HI(i_subband) * &
                  read_table(heat_thin_table,NumheatBin, &
                  tau_pos_in,i_subband,2*i_subband-NumBndin1-1)
             !phi_heat_out_HI = phi_heat_in_HI-phi_heat_HI
             phi_heat_HI = phi_heat_HI/vol

             ! By multiplying here with tau_cell_HeI we already
             ! apply the scaling_HeI factor!
             phi_heat_HeI = NFlux * tau_cell_HeI(i_subband) * &
                  read_table(heat_thin_table,NumheatBin, &
                  tau_pos_in,i_subband,2*i_subband-NumBndin1)
             !phi_heat_out_HeI=phi_heat_in_HeI-phi_heat_HeI
             phi_heat_HeI = phi_heat_HeI/vol

          endif ! Optical depth test
          
          ! Collect contribution from Band 2
          df_heat = phi_heat_HI+phi_heat_HeI

          if (do_secondary_ionizations) then
             ! Band 2, secondary ionizations
             fra_sum1 = f1ion_HI(i_subband)*phi_heat_HI+ &
                  f1ion_HeI(i_subband)*phi_heat_HeI
             fra_sum2 = f2ion_HI(i_subband)*phi_heat_HI+ &
                  f2ion_HeI(i_subband)*phi_heat_HeI
             fra_sum3 = f1heat_HI(i_subband)*phi_heat_HI+ &
                  f1heat_HeI(i_subband)*phi_heat_HeI
             fra_sum4 = f2heat_HI(i_subband)*phi_heat_HI+ &
                  f2heat_HeI(i_subband)*phi_heat_HeI
          
             ! Collect contribution from this band
             df_ion_HeI = y1R(2)*fra_sum1-y2R(2)*fra_sum2  
             df_ion_HI = y1R(1)*fra_sum1-y2R(1)*fra_sum2
             df_heat = df_heat-y1R(3)*fra_sum3+y2R(3)*fra_sum4 
          endif
          
          ! band 3   
       case (NumBndin1+NumBndin2+1:NumBndin1+NumBndin2+NumBndin3)

          ! Band 3, HI
          phi_heat_in_HI = NFlux * read_table(heat_thick_table,NumheatBin, &
               tau_pos_in,i_subband,3*i_subband-NumBndin2-NumBndin1*2-2)
          ! Band 3, HeI
          phi_heat_in_HeI = NFlux * read_table(heat_thick_table,NumheatBin, &
               tau_pos_in,i_subband,3*i_subband-NumBndin2-NumBndin1*2-1)
          ! Band 3, HeII
          phi_heat_in_HeII = NFlux * read_table(heat_thick_table,NumheatBin, &
               tau_pos_in,i_subband,3*i_subband-NumBndin2-NumBndin1*2)

          ! When current cell is optically thick (in total tau)
          if (abs(tau_out_all(i_subband)-tau_in_all(i_subband)) > &
               tau_heat_limit) then

             phi_heat_out_HI = NFlux * read_table(heat_thick_table,NumheatBin, &
               tau_pos_out,i_subband,3*i_subband-NumBndin2-NumBndin1*2-2)
             phi_heat_HI = scaling_HI(i_subband)* &
                  (phi_heat_in_HI-phi_heat_out_HI)/vol

             phi_heat_out_HeI = NFlux * &
                  read_table(heat_thick_table,NumheatBin, &
                  tau_pos_out,i_subband,3*i_subband-NumBndin2-NumBndin1*2-1)
             phi_heat_HeI = scaling_HeI(i_subband)* &
                  (phi_heat_in_HeI-phi_heat_out_HeI)/vol
             
             phi_heat_out_HeII = NFlux * &
                  read_table(heat_thick_table,NumheatBin, &
                  tau_pos_out,i_subband,3*i_subband-NumBndin2-NumBndin1*2)
             phi_heat_HeII = scaling_HeII(i_subband)* &
                  (phi_heat_in_HeII-phi_heat_out_HeII)/vol
             
             ! When current cell is optically thin
          else

             ! By multiplying here with tau_cell_HI we already
             ! apply the scaling_HI factor!
             phi_heat_HI = NFlux * tau_cell_HI(i_subband) * &
                  read_table(heat_thin_table,NumheatBin, &
                  tau_pos_in,i_subband,3*i_subband-NumBndin2-NumBndin1*2-2)
             !phi_heat_out_HI = phi_heat_in_HI-phi_heat_HI
             phi_heat_HI = phi_heat_HI/vol

             ! By multiplying here with tau_cell_HeI we already
             ! apply the scaling_HeI factor!
             phi_heat_HeI = NFlux * tau_cell_HeI(i_subband) * &
                  read_table(heat_thin_table,NumheatBin, &
                  tau_pos_in,i_subband,3*i_subband-NumBndin2-NumBndin1*2-1)
             !phi_heat_out_HeI = phi_heat_in_HeI-phi_heat_HeI
             phi_heat_HeI = phi_heat_HeI/vol

             ! By multiplying here with tau_cell_HeII we already
             ! apply the scaling_HeII factor!
             phi_heat_HeII = NFlux * tau_cell_HeII(i_subband) * &
                  read_table(heat_thin_table,NumheatBin, &
                  tau_pos_in,i_subband,3*i_subband-NumBndin2-NumBndin1*2)
             !phi_heat_out_HeII = phi_heat_in_HeII-phi_heat_HeII
             phi_heat_HeII = phi_heat_HeII/vol

          endif
          
          ! Collect the total contributions to the ionization
          ! and heating rates for this subband
          df_heat=phi_heat_HI+phi_heat_HeI+phi_heat_HeII

          ! Apply effect of secondary ionizations if parameter set
          if (do_secondary_ionizations) then
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
             
             !          if (i_subband == 40) then
             !          endif
             df_ion_HeI = y1R(2)*fra_sum1-y2R(2)*fra_sum2
             df_ion_HI = y1R(1)*fra_sum1-y2R(1)*fra_sum2
             df_heat=df_heat-y1R(3)*fra_sum3+y2R(3)*fra_sum4
          endif
          
       end select
       
       ! Add the subband contribution to the total rates
       f_heat = f_heat + df_heat
       f_ion_HI = f_ion_HI + df_ion_HI
       f_ion_HeI = f_ion_HeI + df_ion_HeI

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
    photrates_add%photo_out = rate1%photo_out + rate2%photo_out
    
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
    rate1%photo_out = 0.0
    
  end subroutine set_photrates_to_zero
  
end module radiation_photoionrates
