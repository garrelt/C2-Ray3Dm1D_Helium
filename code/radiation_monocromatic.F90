module radiation
  
  !     This module contains data and routines which deal with radiative
  !     effects. Its main part deal with photo-ionizing radiation, but it
  !     also initializes other radiative properties, such as cooling (which
  !     are contained in different modules).
  !     It can be used in hydrodynamic or stand-alone radiative transfer 
  !     calculations.
  !
  !     Author: Garrelt Mellema
  ! 
  !     Date: 31-Jan-2008 (02-Jun-2004 (04-Mar-2004)
  
  ! Version
  ! Simplified version
  ! - Only hydrogen
  ! - Option for Grey photo-ionization cross section
  ! - MPI enabled (broadcasts of radiative parameters to all nodes).

  ! Notes:
  ! - the initialization of the radiative cooling does not really belong
  !   here.
  ! - isothermal is sometimes an input parameter, and sometimes a compile
  !   time parameter. This needs to be streamlined. Probably along similar
  !   lines as the stellar parameters are dealt with.

  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use mathconstants, only: pi
  use cgsconstants, only: sigmasb, hplanck, kb, tpic2
  use cgsphotoconstants, only: frth0,frthe0, frthe1, &!, frtop1, frtop2, &
      betah0, betahe, sigh, sighe0, sighe1, ini_He_factors
  use astroconstants, only: R_SOLAR, L_SOLAR
  use romberg, only: scalar_romberg,vector_romberg,romberg_initialisation
  use c2ray_parameters, only: teff_nominal, S_star_nominal, epsilon
  use material, only: isothermal

  implicit none

  integer,parameter :: NumTau=1500  !< Number of table points for the optical depth
  ! This parameter sets the optical depth at the entrance of the grid.
  ! It can be used if radiation enters the simulation volume from the
  ! outside.
  real(kind=dp) :: tauHI=0.0
  real(kind=dp) :: tauHeI=0.0
  real(kind=dp) :: tauHeII=0.0
  real(kind=dp) :: sh0=3.0, she0=2.4,she1=2.8
  real(kind=dp) :: h0ffr!=1.0_dp
  real(kind=dp) :: he0ffr!=1.0_dp
  real(kind=dp) :: he1ffr!=1.0_dp
  character :: sourcetype='B'
  ! Parameters defining the optical depth entries in the table.
  ! minlogtau is log10(lowest optical depth) (table position 1)
  ! maxlogtau is log10(highest optical depth) (table position NumTau)
  ! dlogtau is the step size in log10(tau) between table entries
  real(kind=dp),parameter :: minlogtau=-20.0
  real(kind=dp),parameter :: maxlogtau=4.0
  real(kind=dp),parameter :: dlogtau=(maxlogtau-minlogtau)/real(NumTau)
  logical,parameter :: grey=.false. ! use grey opacities?
  ! stellar properties
  real(kind=dp) :: S_star,rydfactor,teff=0.0,rstar=0.0,lstar=0.0
  
  ! Photo-ionization integral cores
  real(kind=dp):: steph0

  real(kind=dp),dimension(:,:),allocatable   :: hcoreint
  real(kind=dp),dimension(:,:),allocatable   :: hcoreintthin

  real(kind=dp),dimension(:),allocatable :: coreint
  real(kind=dp),dimension(:),allocatable :: coreintthin

  real(kind=dp),dimension(:,:),allocatable   :: hphotint  
  real(kind=dp),dimension(:,:),allocatable   :: hphotintthin

  real(kind=dp),dimension(:),allocatable   :: photint  
  real(kind=dp),dimension(:),allocatable   :: photintthin

  ! This type contains all the photo-ionization rates
  ! The in and out rates are used to ensure photon-conservation.
  ! See the C2-Ray paper. 

  type photrates    
     real(kind=dp) :: h          !< total H ionizing rate           
     real(kind=dp) :: he(0:1)        !< total He0 ionizing rate    
     real(kind=dp),dimension(0:6) :: int_out   !< ionizing rate in Int1:Int3d out
     real(kind=dp) :: hv_h       !< total H heating rate      
     real(kind=dp) :: hv_he(0:1)     !< total He0 heating rate          
     real(kind=dp) :: h_in       !< H ionizing in-rate       
     real(kind=dp) :: he_in(0:1)     !< He0 ionizing in-rate       
     real(kind=dp) :: hv_h_in    !< H in-heating rate     
     real(kind=dp) :: hv_he_in(0:1)  !< He0 in-heating rate        
     real(kind=dp) :: h_out      !< H ionizing out-rate       
     real(kind=dp) :: he_out(0:1)    !< He0 ionizing out-rate       
     real(kind=dp) :: hv_h_out   !< H out-heating rate     
     real(kind=dp) :: hv_he_out(0:1) !< He0 out-heating rate      
  end type photrates

#ifdef MPI       
    integer,private :: ierror
#endif

contains

!=======================================================================

  subroutine rad_ini ()

    ! initializes constants and tables for radiation processes
    ! (heating, cooling and ionization)

    use radiative_cooling, only: setup_cool
    ! Initialize integration routines
!    call romberg_initialisation(NumFreq)

    ! Ask for the parameters of the spectrum
    call spectrum_parms ()

    ! Calculate spectral integral cores
    call spec_integr_cores ()

    ! Find the photo-ionization integrals for this spectrum
    call spec_integr ()

    call ini_He_factors ()

    ! Setup cooling is now done in C2ray directly
    if (.not.isothermal) call setup_cool () ! SHOULD BE CALLED ELSEWHERE

  end subroutine rad_ini

  !=======================================================================

  subroutine spectrum_parms

    ! Input routine: establish the ionizing spectrum     
    ! Author: Garrelt Mellema
    ! Update: 18-Feb-2004

    use file_admin, only: stdinput
    
    integer :: nchoice
    real(kind=dp) :: totflux

    ! Ask for input
    if (rank == 0 .and. teff_nominal == 0.0) then
          write(*,'(A,$)') 'Give S_* (ionizing photons s^-1): '
          read(stdinput,*) S_star
          write(*,'(A,$)') 'Give frequency in rybberg: '
          read(stdinput,*) rydfactor
    endif

#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(rydfactor,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         ierror)
    call MPI_BCAST(S_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         ierror)
#endif

  end subroutine spectrum_parms

  !=======================================================================
  
  !=======================================================================

  subroutine spec_integr_cores ()

    integer :: i,n, n2, nfrq,m
    real(kind=dp) :: tau(0:NumTau)
    real(kind=dp) :: fr, xi, yi, Fi
    real(kind=dp) :: Eth(1:3),Emax(1:3),E0(1:3),sig0(1:3) 
    real(kind=dp),target ::     sig(1:3)
    real(kind=dp) ::   ya(1:3) ,P(1:3), yw(1:3), y0(1:3),y1(1:3),E
    real(kind=dp) :: Mb=1e-18  !Mb in cm^2
    real(kind=dp),pointer :: myffr
!-----------------------------------------------------------
!-----------------------------------------------------------  
    ! Allocate the spectral integral cores
    allocate(coreint(0:NumTau))     ! for 7 intervals
    allocate(coreintthin(0:NumTau)) ! for 7 intervals 
!---THERMAL-------------------------------------------------
    if (.not.isothermal) then
    allocate(hcoreint(0:NumTau,1:3))       !for 1 interval a 1 species, 2 intervals a 2 species
    allocate(Hcoreintthin(0:NumTau,1:3))   !and 4 intervals a 3 species
    endif
!-----------------------------------------------------------
          	Eth=(/13.6_dp, 24.59_dp, 54.42_dp/)
          	Emax=(/5.0e4_dp, 5.0e4_dp, 5.0e4_dp/)
          	E0=(/0.4298_dp, 13.61_dp, 1.720_dp/)
		sig0=(/5.475e4_dp, 9.492e2_dp, 1.369e4_dp/)*Mb
		ya=(/3.288e1_dp, 1.469_dp, 3.288e1_dp/)
		P=(/2.963_dp, 3.188_dp, 2.963_dp/)
		yw=(/0.0_dp, 2.039_dp, 0.0_dp/)
		y0=(/0.0_dp, 4.434e-1_dp, 0.0_dp/)
		y1=(/0.0_dp, 2.136_dp, 0.0_dp/)

    ! fill the optical depth array used to fill the tables 
    ! it is filled in NumTau logarithmic steps 
    ! from minlogtau to maxlogtau
    do n=1,NumTau
       tau(n)=10.0**(minlogtau+dlogtau*real(n-1))
    enddo

    ! Position zero corresponds to zero optical depth
    tau(0)=0.0
    E=13.6_dp*rydfactor
    ! Warn about grey opacities:
    if (grey .and. rank == 0) write(logf,*) 'WARNING: Using grey opacities'

          fr=frth0*rydfactor   ! set here which frequency I want to use
       
       m=1   
       if (E.ge.Eth(2))  m=2 
       if (E.ge.Eth(3))  m=3        		
          		
          		
	    do i=1,m
    	 	xi=E/E0(i)-y0(i) 
    	 	yi=sqrt(xi*xi+y1(i)*y1(i))  
    	 	Fi=((xi-1.0_dp)*(xi-1.0_dp)+yw(i)*yw(i))*yi**(0.5_dp*P(i)-5.5_dp)*(1.0_dp+sqrt(yi/ya(i)))**(-P(i))
    	 	sig(i)=sig0(i)*Fi
          enddo    		
                        
      if (E.ge.Eth(3)) then
         m=3
         myffr => sig(3)
      elseif (E.ge.Eth(2)) then
         m=2
         sig(3)=0.0_dp         
         myffr => sig(2)
      else 
         m=1
         sig(2:3)=0.0_dp         
         myffr => sig(1)
      endif
                   
           h0ffr=sig(1)
           he0ffr=sig(2)
           he1ffr=sig(3)


          do n=0,NumTau
           if (tau(n)*he1ffr < 700.0) then  
                coreint(n)    = S_star*          exp(-tau(n))!*myffr)  
                coreintthin(n)= S_star*exp(-tau(n))!*myffr)  
	     else
                coreint(n)    =0.0
                coreintthin(n)=0.0!
	     endif

            if (.not.isothermal) then
                hcoreint(n,1)    =hplanck*(fr-frth0) *coreint(n)
                hcoreint(n,2)    =hplanck*(fr-frthe0)*coreint(n)
                hcoreint(n,3)    =hplanck*(fr-frthe1)*coreint(n)
                hcoreintthin(n,1)=hplanck*(fr-frth0) *coreintthin(n)
                hcoreintthin(n,2)=hplanck*(fr-frthe0)*coreintthin(n)
                hcoreintthin(n,3)=hplanck*(fr-frthe1)*coreintthin(n)
            endif
          enddo	 

  end subroutine spec_integr_cores
  
  ! =======================================================================
  
  subroutine spec_integr ()

    ! Calculates photo ionization integrals

    ! Author: Garrelt Mellema
    ! Date: 19-Feb-2004

    ! Version: Simplified from Coral

    ! Two types of integrals are evaluated: one for optically thick cells
    ! (hphot, hheat) and one for optically thin cells (hphot1, hheat1).

    integer :: i,n,nfrq, n2,j
    real(kind=dp) :: rstar2,rfr,frmax
    real(kind=dp) :: fr
    real(kind=dp) :: weight(0:NumTau),phot(0:NumTau)

! Allocate photo-ionization tables
    allocate(photint(0:NumTau))
    allocate(photintthin(0:NumTau))

    if (.not.isothermal) then
	allocate(hphotint(0:NumTau,1:3))
    	allocate(hphotintthin(0:NumTau,1:3))
    endif

 ! *** optically thick ***

       fr=frth0*rydfactor    !use the frequency I set somewhere earlier
       do n=0,NumTau
          weight(n)=steph0
          photint(n)=coreint(n)
          if (.not.isothermal) then
          hphotint(n,1)=hcoreint(n,1)
          hphotint(n,2)=hcoreint(n,2)
          hphotint(n,3)=hcoreint(n,3)
          endif
       enddo

!   write(*,*) photint
!-----------------------------------------------------------   
 ! *** optically thin ***
       do n=0,NumTau
          photintthin(n)=coreintthin(n)
          if (.not.isothermal) then
          hphotintthin(n,1)=hcoreintthin(n,1)
          hphotintthin(n,2)=hcoreintthin(n,2)
          hphotintthin(n,3)=hcoreintthin(n,3)
          endif
       enddo
!-----------------------------------------------------------
!-----------------------------------------------------------  
    ! Deallocate the cores
    deallocate(coreint)
    deallocate(coreintthin)
    if (.not.isothermal) then
    	deallocate(hcoreint)
    	deallocate(hcoreintthin)
    endif
  end subroutine spec_integr

  ! =======================================================================
  subroutine photoion (phi,hcolum_in,hcolum_out,hecolum_in,hecolum_out,&
	       	      vol,nsrc,i_state)
    ! Calculates photo-ionization rates
    
    ! Author: Garrelt Mellema
    ! Date: 11-May-2005 (f90) (18 feb 2004
    
    ! Version:
    ! Simplified version derived from Coral version, for testing
    ! photon conservation. Only hydrogen is dealt with, and
    ! one frequency band is used.
!*** I HAVE TO INCLUDE THE MULTIOPLICATION WITH NSRC!
    !use sourceprops, only: NormFlux
    use cgsphotoconstants
    type(photrates),intent(out) :: phi
    real(kind=dp),intent(in) :: hcolum_in,  hcolum_out,  vol,i_state
    integer,intent(in) :: nsrc !< number of the source

    real(kind=dp),dimension(0:1),intent(in) :: hecolum_in,hecolum_out !(0:1)
    real(kind=dp) ::  hcolum_cell
    real(kind=dp), dimension(0:1):: hecolum_cell

    real(kind=dp) :: tauh_in,  tauh_out, phi_hv_h, phi_hv_he0, phi_hv_he1
    real(kind=dp) :: tauhe0_in,tauhe0_out, f_heat, f_ion_h, f_ion_he0
    real(kind=dp) :: tauhe1_in,tauhe1_out

    real(kind=dp) :: tau_i,      tau_o
    real(kind=dp) :: odpos_i,    odpos_o
    real(kind=dp) :: residual_i, residual_o
    integer     :: ipos_i,     ipos_o 
    integer       :: ipos_p1_i,  ipos_p1_o
    real(kind=dp) ::NormFlux(1)=1.0

!-----------------------------------------------------------
! For the positions of n, n2 osv 
    real(kind=dp) :: N_He0_ov_H0,    N_He1_ov_He0,    N_H0_ov_He0,    N_He1_ov_H0    ,N_H0_ov_He1,    N_He0_ov_He1

    integer       ::  n
    real(kind=dp) :: delta_tau_out
    real(kind=dp) :: delta_tau_in
    real(kind=dp), dimension(3):: delt_tau_o


    real(kind=dp) ::  phi_in,   phi_out, phi_tot     

    real(kind=dp) :: NFlux
    real(kind=dp),dimension(3) :: scaling
	hcolum_cell=hcolum_out-hcolum_in
	hecolum_cell(:)=hecolum_out(:)-hecolum_in(:)

        N_He0_ov_H0         = max(epsilon,hecolum_cell(0))/max(epsilon,hcolum_cell)
        N_He0_ov_He1        = max(epsilon,hecolum_cell(0))/max(epsilon,hecolum_cell(1))
        N_H0_ov_He0         = max(epsilon,hcolum_cell)  /max(epsilon,hecolum_cell(0))
        N_H0_ov_He1         = max(epsilon,hcolum_cell)  /max(epsilon,hecolum_cell(1))
        N_He1_ov_He0        = max(epsilon,hecolum_cell(1))/max(epsilon,hecolum_cell(0))
        N_He1_ov_H0         = max(epsilon,hecolum_cell(1))/max(epsilon,hcolum_cell)
        

    ! find the optical depths (in and outgoing) at the corresponding threshold values
   
    tauh_in   =h0ffr *hcolum_in
    tauh_out  =h0ffr *hcolum_out  
    tauhe0_in =he0ffr*hecolum_in(0)
    tauhe0_out=he0ffr*hecolum_out(0)
    tauhe1_in =he1ffr*hecolum_in(1)
    tauhe1_out=he1ffr*hecolum_out(1)   

    delta_tau_in=tauh_in+tauhe0_in+tauhe1_in

    ! find the table positions for the optical depth (ingoing)--------------------!
      tau_i      =log10(max(1.0e-20_dp,delta_tau_in))                           !
    	odpos_i    =min(real(NumTau,dp),max(0.0_dp,1.0+(tau_i-minlogtau)/dlogtau))  !
    	ipos_i     =int(odpos_i)                                                    !
    	residual_i =odpos_i-real(ipos_i,dp)                                         ! 
    	ipos_p1_i  =min(NumTau,ipos_i+1)                                            ! 
    !-----------------------------------------------------------------------------!  
!write(*,*) ipos_i, tauh_in, tauh_out
!pause
! What is the sigma ratio at the chosen frequency?  !No, should be as on paper 0, no? 
! Yes, ratio at of H(nu=nuHe1)/H(nu=nuH) and He0(nu=nuHe1)/He0(nu=nuHe0,fit)


    delta_tau_out=tauh_out+tauhe0_out+tauhe1_out

   ! find the table positions for the optical depth (outgoing)
   !find the table positions-------------------------------------------------------!
      tau_o      =log10(max(1.0e-20_dp,delta_tau_out))                           !
    	odpos_o    =min(real(NumTau,dp),max(0.0_dp,1.0+(tau_o-minlogtau)/dlogtau)) !
    	ipos_o     =int(odpos_o)                                                   !
    	residual_o =odpos_o-real(ipos_o,dp)                                        !
    	ipos_p1_o  =min(NumTau,ipos_o+1)                                           !
   !-------------------------------------------------------------------------------!
    NFlux=NormFlux(nsrc)  

   phi_in= (photint(ipos_i)+(photint(ipos_p1_i)  -photint(ipos_i))*residual_i)*NFlux

  scaling(1)= (tauh_out-tauh_in)/(delta_tau_out-delta_tau_in)
  scaling(2)= (tauhe0_out-tauhe0_in)/(delta_tau_out-delta_tau_in)
  scaling(3)= (tauhe1_out-tauhe1_in)/(delta_tau_out-delta_tau_in)

 !write(*,*) scaling(1),scaling(2),scaling(3)

 !*********total        
       if (abs(delta_tau_out-delta_tau_in).gt.-1.0e-9) then
          phi_out=    (photint(ipos_o) +(photint(ipos_p1_o)  -photint(ipos_o)) *residual_o)*NFlux 
       else
          phi_tot=((photintthin(ipos_i)+(photintthin(ipos_p1_i)-photintthin(ipos_p1_i))*residual_i)&
             *(delta_tau_out-delta_tau_in))*NFlux
          phi_out=phi_in-phi_tot
       endif

    phi_tot=phi_in-phi_out
    phi%h=scaling(1)*phi_tot/vol
    phi%he(0)=scaling(2)*phi_tot/vol
    phi%he(1)=scaling(3)*phi_tot/vol
    phi%int_out=phi_out

!---------------------------------------NOT ISOTHERMAL------
    if(.not.isothermal) then 
       phi%hv_h_in  =(hphotint(ipos_i,1)  +(hphotint(ipos_p1_i,1) -hphotint(ipos_i,1))*residual_i)*NFlux 
       phi%hv_he_in(0)=(hphotint(ipos_i,2)  +(hphotint(ipos_p1_i,2) -hphotint(ipos_i,2))*residual_i)*NFlux 
       phi%hv_he_in(1)=(hphotint(ipos_i,3)  +(hphotint(ipos_p1_i,3) -hphotint(ipos_i,3))*residual_i)*NFlux  
	!********H
          if (abs(delt_tau_o(1)-delta_tau_in).gt.1.0e-2) then
          	phi%hv_h_out= (hphotint(ipos_o,1)+&
                              (hphotint(ipos_p1_o,1)-hphotint(ipos_o,1))*residual_o)*NFlux
	        phi_hv_h=scaling(1)*(phi%hv_h_in-phi%hv_h_out)/vol
          else
    		phi_hv_h =scaling(1)*((hphotintthin(ipos_i,1)+&
                          (hphotintthin(ipos_p1_i,1)-hphotintthin(ipos_p1_i,1))*residual_i)&
                          *(delta_tau_out-delta_tau_in))/vol*NFlux
          endif
! 	!********He0      
          if (abs(delt_tau_o(2)-delta_tau_in).gt.1.0e-2) then 
          	phi%hv_he_out(0)=(hphotint(ipos_o,2)+&
			    (hphotint(ipos_p1_o,2)-hphotint(ipos_o,2))*residual_o)*NFlux
	        phi_hv_he0 =scaling(2)*(phi%hv_he_in(0)-phi%hv_he_out(0))/vol
          else
     		phi_hv_he0=scaling(2)*((hphotintthin(ipos_i,2)+&
			   (hphotintthin(ipos_p1_i,2)-hphotintthin(ipos_p1_i,2))*residual_i) &
                           *(delta_tau_out-delta_tau_in))/vol*NFlux
          endif

! 	!********He1       
          if (abs(delt_tau_o(3)-delta_tau_in).gt.1.0e-2) then
         	phi%hv_he_out(1)=  (hphotint(ipos_o,3)+&
			       (hphotint(ipos_p1_o,3)-hphotint(ipos_o,3))*residual_o)*NFlux
                phi_hv_he1  =  scaling(3)*(phi%hv_he_in(1)-phi%hv_he_out(1))/vol
          else 
          	phi_hv_he1 =scaling(3)*((hphotintthin(ipos_i,3)+&
			    (hphotintthin(ipos_p1_i,3)-hphotintthin(ipos_p1_i,3)) *residual_i) &
                            *(delta_tau_out- delta_tau_in))/vol*NFlux
          endif

    endif
!-----------------------------------------------------------

! probably I should actually distingues between "primary electrons" above and below 100 eV. Below that, I do have a 
! dependence on the energy, above that I do not. 

    if(.not.isothermal) then 
   !fractions taken from Shull and van Steenberg 1985
! Since the primary electrons from H-ionization have energiese above 50 eV, I use the fitting formular
! He1 have energies below 10 eV, so everything goes into heat.
! He0 have energies around 40 eV, so I use Shull and van Steenberg but only for H, not He0.

    phi_hv_he0=phi_hv_h+phi_hv_he0
    f_heat=0.9971_dp*(1.0_dp-(1.0_dp-i_state**0.2663_dp)**1.3163_dp)
    f_ion_H=0.3908_dp*(1.0_dp-i_state**0.4092_dp)**1.7592_dp
    f_ion_He0=0.0554_dp*(1.0_dp-i_state**0.4614_dp)**1.6660_dp*0.1_dp
      
! This is to test the effect of secondary ionizations
!    f_heat=1.0_dp
!    f_ion_H=0.0_dp
!    f_ion_He0=0.0_dp
   

          phi%hv_h=phi_hv_he0*f_heat+phi_hv_he1   
          phi%h=phi%h+phi_hv_he0*f_ion_H/(frth0*hplanck)!/vol
          phi%he(0)=phi%he(0)+phi_hv_he0*f_ion_He0/(frthe0*hplanck)!/vol

  endif
!-----------------------------------------------------------
  end subroutine photoion
! =======================================================================
end module radiation
