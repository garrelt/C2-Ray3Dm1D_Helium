!     This module sets the sizes of radiation quantities

module radiation_sizes

  use precision, only: dp
  use cgsphotoconstants, only: ion_freq_HI,&             ! HI ionization energy in frequency
                               ion_freq_HeI,&            ! HeI ionization energy in frequency
                               ion_freq_HeII,&           ! HeII ionization energy in frequency
                               sigma_HI_at_ion_freq,&    ! HI cross section at its ionzing frequency
                               sigma_HeI_at_ion_freq,&   ! HeI cross section at its ionzing frequency
                               sigma_HeII_at_ion_freq    ! HeII cross section at its ionzing frequency
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

contains

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

end module radiation_sizes
