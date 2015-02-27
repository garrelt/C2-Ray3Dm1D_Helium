module doric_module

  use precision, only: dp
  use abundances, only: abu_he
  use parameter, only: do_couple

  implicit none

contains

  subroutine doric (dt,ndens_electron,&
                    begin_xHI, begin_xHII, begin_xHeI, begin_xHeII, begin_xHeIII, &
                    avg_xHI, avg_xHII, avg_xHeI, avg_xHeII, avg_xHeIII, &
                    end_xHI, end_xHII, end_xHeI, end_xHeII, end_xHeIII, & 
                    temperature,photoionization_HI_rate,photoionization_HeI_rate,photoionization_HeII_rate,&
                    y,z,y2a,y2b)



    use cgsconstants, only: bh00,bhe00, bhe10, albpow,alcpow, &
         colh0,colhe, temphe,temph0,ev2fr
    use cgsphotoconstants, only: sigma_HI_at_ion_freq, sigma_HeI_at_ion_freq
    use file_admin, only: stdinput
    use clumping, only: clumping_factor
    use abundances, only: abu_he
    use type_definition, only: photrates, ionstates
    use tped, only: electrondens ! should this really be used inside doric?
    use c2ray_parameters, only: epsilon_value

    implicit none

    real(kind=dp),intent(in)    :: dt     !< time step
    real(kind=dp),intent(inout) :: ndens_electron    !< electron density
    real(kind=dp),intent(in) :: begin_xHI
    real(kind=dp),intent(in) :: begin_xHII
    real(kind=dp),intent(in) :: begin_xHeI
    real(kind=dp),intent(in) :: begin_xHeII
    real(kind=dp),intent(in) :: begin_xHeIII
    real(kind=dp),intent(inout) :: avg_xHI
    real(kind=dp),intent(inout) :: avg_xHII
    real(kind=dp),intent(inout) :: avg_xHeI
    real(kind=dp),intent(inout) :: avg_xHeII
    real(kind=dp),intent(inout) :: avg_xHeIII
    real(kind=dp),intent(inout) :: end_xHI
    real(kind=dp),intent(inout) :: end_xHII
    real(kind=dp),intent(inout) :: end_xHeI
    real(kind=dp),intent(inout) :: end_xHeII
    real(kind=dp),intent(inout) :: end_xHeIII
    real(kind=dp),intent(in) :: temperature
    real(kind=dp),intent(inout)   :: y
    real(kind=dp),intent(inout)   :: z  
    real(kind=dp),intent(inout) :: y2a
    real(kind=dp),intent(inout) :: y2b
    real(kind=dp),intent(in) :: photoionization_HI_rate
    real(kind=dp),intent(in) :: photoionization_HeI_rate
    real(kind=dp),intent(in) :: photoionization_HeII_rate
    real(kind=dp) ::alpha_HI_B,alpha_HeI_1,alpha_HeI_B,alpha_HeI_A,alpha_HeII_B, &
                    alpha_HeII_2,alpha_HI_A,alpha_HeII_A,alpha_HeII_1 
    real(kind=dp) :: p, f, w, l, m
    real(kind=dp) :: lambda_1, lambda_2, lambda_3
    real(kind=dp) :: lambda_1_dt, lambda_2_dt, lambda_3_dt
    real(kind=dp) :: exp_lambda_1_dt,exp_lambda_2_dt,exp_lambda_3_dt
    real(kind=dp) :: x2_1, x2_2, x3_1, x3_2
    real(kind=dp) :: c_1,c_2,c_3,Scoef,Rcoef,Tcoef, xp_1,xp_2,xp_3,Kcoef, &
                     Bcoef,Bminus,Bplus
    real(kind=dp) :: avg_factor_1, avg_factor_2, avg_factor_3
    real(kind=dp) :: A_11,A_12,A_13,A_22,A_23,A_33,two_A_32 
    real(kind=dp) :: g_1, g_2, A_32
    real(kind=dp) :: heliumfraction,  normfac ! calculate as abu_he/(1-abu_he)

   !> Hydrogen 0 A recombination parameter
   real(kind=dp) :: arech0
   !> Hydrogen 0 B recombination parameter
   real(kind=dp) :: brech0


   !> Helium   0 A recombination parameter
   real(kind=dp) :: areche0
   !> Helium   0 B recombination parameter
   real(kind=dp) :: breche0
   !> Helium   0 1 recombination parameter
   real(kind=dp) :: oreche0

   !> Helium   + A recombination parameter
   real(kind=dp) :: areche1
   !> Helium   + B recombination parameter
   real(kind=dp) :: breche1
   !> Helium   + 2 recombination parameter
   real(kind=dp) :: treche1


   !> H0 collisional ionization parameter at T=temp0
   real(kind=dp) :: colli_HI
   !> He0 collisional ionization parameter at T=temp0
   real(kind=dp) :: colli_HeI
   !> He1 collisional ionization parameter at T=temp0
   real(kind=dp) :: colli_HeII
   !> Fraction fo He++ -> He+ recombination photons that goes into 2photon decay
   real(kind=dp) :: v_factor
    real(kind=dp) :: sqrtt0,lambda_HI,lambda_HeI,lambda_HeII,dielectronic

     lambda_HI =  2.0_dp*(temph0/temperature)                                           
     arech0=1.269e-13*lambda_HI**1.503_dp/(1.0_dp+(lambda_HI/0.522)**0.470)**1.923   
     brech0=2.753e-14*lambda_HI**1.500_dp/(1.0_dp+(lambda_HI/2.740)**0.407)**2.242    
                                          
    lambda_HeII=2.0_dp*(temphe(1)/temperature)                                                  
    breche1=5.5060e-14_dp*lambda_HeII**1.5_dp/(1.0_dp+(lambda_HeII/2.740_dp)**0.407_dp)**2.242_dp 
    areche1=2.538e-13*lambda_HeII**1.503_dp/(1.0_dp+(lambda_HeII/0.522_dp)**0.470_dp)**1.923_dp
    treche1 = 3.4e-13_dp*(temperature/1.0e4_dp)**(-0.6_dp) 

    if (temperature < 9.e3_dp) then   ! was 1.5e4
        lambda_HeI=2.0_dp*(temph0/temperature)
        areche0=1.269e-13_dp*lambda_HeI**1.503_dp/(1.0_dp+(lambda_HeI/0.522)**0.470)**1.923
        breche0=2.753e-14_dp*lambda_HeI**1.500_dp/(1.0_dp+(lambda_HeI/2.740)**0.407)**2.242
    else
        lambda_HeI=   2.0_dp*(temphe(0)/temperature)
        dielectronic = 1.9e-3_dp*temperature**(-1.5_dp)* &
             exp(-4.7e5_dp/temperature)*(1.0_dp+0.3_dp*exp(-9.4e4_dp/temperature))
        areche0= 3.000e-14_dp*lambda_HeI**0.654_dp + dielectronic
        breche0= 1.260e-14_dp*lambda_HeI**0.750_dp + dielectronic
    endif
        oreche0=areche0-breche0

    sqrtt0 =sqrt(temperature)
    colli_HI =colh0*sqrtt0*exp(-temph0/temperature)
    colli_HeI=colhe(0)*sqrtt0*exp(-temphe(0)/temperature)
    colli_HeII=colhe(1)*sqrtt0*exp(-temphe(1)/temperature)

    v_factor=0.285_dp*(temperature/1.0e4_dp)**0.119_dp 

    ! see Osterbrock, 1989
    p = 0.96_dp 
    ! fraction of photons from 2-photon decay, energetic enough to ionize hydrogen
    l = 1.425_dp 
    ! fraction of photons from 2-photon decay, energetic enough to ionize neutral helium
    m = 0.737_dp 
    heliumfraction=abu_he/(1.0_dp-abu_he)    ! He/H
    f= max(min(10.0_dp*end_xHI,1.0_dp),0.01_dp) 
    ! These numbers come from Flower & Perinotto (1980)
    w = (l-m)+m*y  

    alpha_HI_B   = clumping_factor* brech0
    alpha_HI_A   = clumping_factor* arech0
    alpha_HeI_1  = clumping_factor* oreche0
    alpha_HeI_B  = clumping_factor* breche0    
    alpha_HeI_A  = clumping_factor* areche0
    alpha_HeII_B = clumping_factor* breche1
    alpha_HeII_A = clumping_factor* areche1
    alpha_HeII_2 = clumping_factor* treche1
    alpha_HeII_1 =  alpha_HeII_A - alpha_HeII_B

    ! elements of matrix A and vector g
    g_1=max(photoionization_HI_rate+ndens_electron*colli_HI,1.0e-200_dp)    
    g_2=max(photoionization_HeI_rate+ndens_electron*colli_HeI,1.0e-200_dp)
    A_32=max(photoionization_HeII_rate+ndens_electron*colli_HeII,1.0e-200_dp)

    A_11 = -(g_1+ndens_electron*alpha_HI_B)

    if (do_couple) then


    A_12 = (y*alpha_HeI_1+p*alpha_HeI_B)*ndens_electron*heliumfraction
    A_13 = ((f*z*(1.0_dp-v_factor) +v_factor*w)*alpha_HeII_B +alpha_HeII_2+ &
           (1.0_dp-y2a-y2b)*alpha_HeII_1)*heliumfraction*ndens_electron  
    A_22 = -g_2-A_32-ndens_electron*(alpha_HeI_A-(1.0_dp-y)*alpha_HeI_1)
    A_33 = -ndens_electron*(alpha_HeII_A-y2a*alpha_HeII_1)             
    A_23 = -g_2+ndens_electron*alpha_HeII_B*(f*(1.0_dp-z)*(1.0_dp-v_factor)+ &
           v_factor*(1.425_dp-w))-A_33+alpha_HeII_1*y2b*ndens_electron

    else

    A_12 = 0.0_dp
    A_13 = 0.0_dp
    A_22 = -(g_2+A_32+ndens_electron*alpha_HeI_B)
    A_33 = -ndens_electron*alpha_HeII_B       
    A_23 = -g_2-A_33
 
    endif

    ! this is just abbrivations
    two_A_32 = 2.0_dp*A_32

    ! These are just abbrivations 
    Bcoef = A_33-A_22
    Scoef = sqrt(Bcoef*Bcoef+4.0_dp*A_32*A_23)
    Kcoef = 1.0_dp/(A_23*A_32-A_33*A_22)
    Bminus = Bcoef-Scoef
    Bplus = Bcoef+Scoef

    ! These are the Eigenvalues of the Matrix
    lambda_1 = A_11
    lambda_2 = 0.5_dp*(A_33+A_22-Scoef)
    lambda_3 = 0.5_dp*(A_33+A_22+Scoef)

    ! elements of particular solution vector
    xp_1 = -1.0_dp/A_11*(g_1+(A_12*A_33-A_13*A_32)*(g_2*Kcoef))   
    xp_2 = g_2*(A_33*Kcoef)                        
    xp_3 = -g_2*(A_32*Kcoef)  

    ! elements of some column vectors
    x2_1 = -A_13/(A_11-lambda_2)+(A_12/two_A_32)*Bplus/(A_11-lambda_2)  
    x2_2 = (-Bplus)/(two_A_32)
    x3_1 = (-two_A_32 *A_13+A_12*(Bminus))/(two_A_32*(A_11-lambda_3))  
    x3_2 = (-Bminus)/(two_A_32)

    ! These are just abbrivations
    Rcoef = two_A_32 *(xp_2-begin_xHeII)
    Tcoef = xp_3-begin_xHeIII

    ! coefficients of homogeneous solutions
    c_1 = -xp_1+ (x3_1-x2_1)*(Rcoef/(2.0_dp*Scoef))+ &
          Tcoef*((Bplus*x3_1/(2.0_dp*Scoef)-Bminus*x2_1/(2.0_dp*Scoef)))+ &
          begin_xHII
    c_2 = (Rcoef+(Bminus)*Tcoef)/(2.0_dp*Scoef)
    c_3 = -(Rcoef+(Bplus)*Tcoef)/(2.0_dp*Scoef)

    ! arguments of exponential functions
    lambda_1_dt = dt*lambda_1
    lambda_2_dt = dt*lambda_2
    lambda_3_dt = dt*lambda_3

    ! exponential functions of homogeneous solutions
    exp_lambda_1_dt = exp(lambda_1_dt)
    exp_lambda_2_dt = exp(lambda_2_dt)
    exp_lambda_3_dt = exp(lambda_3_dt)

    ! ionization fractions at the end of the time-step
    end_xHII = c_1*exp_lambda_1_dt+c_2*exp_lambda_2_dt*x2_1+c_3*exp_lambda_3_dt*x3_1+xp_1 
    end_xHeII = c_2*exp_lambda_2_dt*x2_2+c_3*exp_lambda_3_dt*x3_2+xp_2
    end_xHeIII = c_2*exp_lambda_2_dt+c_3*exp_lambda_3_dt+xp_3
    end_xHI = 1.0_dp-end_xHII
    end_xHeI= 1.0_dp-end_xHeII-end_xHeIII

    if (end_xHI < epsilon_value) then
       end_xHI=epsilon_value
       end_xHII=1.0_dp-epsilon_value
    endif
    if (end_xHII < epsilon_value) then
       end_xHII =epsilon_value
       end_xHI=1.0_dp-epsilon_value
    endif

    if ((end_xHeI.le.epsilon_value).or.(end_xHeII.le.epsilon_value).or. &
         (end_xHeIII.le.epsilon_value))then
       if (end_xHeI .lt. epsilon_value) then
          end_xHeI=epsilon_value
       endif

       if (end_xHeII .lt. epsilon_value) then
          end_xHeII=epsilon_value      
       endif
       if (end_xHeIII .lt. epsilon_value) then
          end_xHeIII=epsilon_value
       endif

       normfac=end_xHeI+end_xHeII+end_xHeIII
       end_xHeI=end_xHeI/normfac
       end_xHeII=end_xHeII/normfac
       end_xHeIII=end_xHeIII/normfac
    endif

    if (abs(lambda_1_dt).lt.1.0e-8) then
       avg_factor_1=c_1
    else
       avg_factor_1 = c_1*(exp_lambda_1_dt- 1.0_dp)/lambda_1_dt
!if (isnan(avg_factor_1)) write(*,*) 'avg_factor_1'
    endif

    if (abs(lambda_2_dt).lt.1.0e-8) then
       avg_factor_2=c_2
    else
       avg_factor_2 = c_2*(exp_lambda_2_dt- 1.0_dp)/lambda_2_dt 
!if (isnan(avg_factor_2)) write(*,*) 'avg_factor_2'      
    endif

    if (abs(lambda_3_dt).lt.1.0e-8) then
       avg_factor_3=c_3
    else
       avg_factor_3 = c_3*(exp_lambda_3_dt- 1.0_dp)/lambda_3_dt
!if (isnan(avg_factor_3)) write(*,*) 'avg_factor_3'
    endif

    avg_xHII  = xp_1 + avg_factor_1 + x2_1*avg_factor_2 + x3_1*avg_factor_3
    avg_xHeII = xp_2 +                x2_2*avg_factor_2 + x3_2*avg_factor_3
    avg_xHeIII = xp_3 +                       avg_factor_2 +        avg_factor_3
    avg_xHI  = 1.0_dp-avg_xHII
    avg_xHeI = 1.0_dp-avg_xHeII-avg_xHeIII

    if (avg_xHII .lt. epsilon_value) then
       avg_xHII=epsilon_value
       avg_xHI=1.0_dp-epsilon_value
    endif
    if (avg_xHI .lt. epsilon_value) then
       avg_xHI=epsilon_value
       avg_xHII=1.0_dp-epsilon_value
    endif

    if ((avg_xHeI.le.epsilon_value).or.(avg_xHeII.le.epsilon_value).or. &
         (avg_xHeIII.le.epsilon_value))then

       if (avg_xHeII .lt. epsilon_value) avg_xHeII=epsilon_value   
       if (avg_xHeIII .lt. epsilon_value) avg_xHeIII=epsilon_value
       if (avg_xHeI .lt. epsilon_value) avg_xHeI=epsilon_value
       normfac=avg_xHeI+avg_xHeII+avg_xHeIII
       avg_xHeI=avg_xHeI/normfac
       avg_xHeII=avg_xHeII/normfac
       avg_xHeIII=avg_xHeIII/normfac
!if (isnan(avg_xHeI)) write(*,*) 'avg_xHeI'
!if (isnan(avg_xHeII)) write(*,*) 'avg_xHeII'
!if (isnan(avg_xHeIII)) write(*,*) 'avg_xHeIII'
    endif

    return
  end subroutine doric

  ! ===========================================================================

  subroutine prepare_doric_factors(xHI,xHeI,xHeII,number_density_atom,yfrac,zfrac,y2afrac,y2bfrac)

    use cgsphotoconstants, only: sigma_H_heth, sigma_H_heLya,sigma_He_heLya
    use cgsphotoconstants, only: sigma_HeI_at_ion_freq, sigma_He_he2
    use cgsphotoconstants, only: sigma_HeII_at_ion_freq,sigma_He_he2,sigma_H_he2

    implicit none

    real(kind=dp),intent(in) :: xHI 
    real(kind=dp),intent(in) :: xHeI 
    real(kind=dp),intent(in) :: xHeII 
    real(kind=dp),intent(in) :: number_density_atom
    real(kind=dp),intent(out) :: yfrac,zfrac
    real(kind=dp),intent(out) :: y2afrac,y2bfrac
    real(kind=dp) :: fHI
    real(kind=dp) :: fHeI
    real(kind=dp) :: fHeII
    real(kind=dp) :: tau_H_heth
    real(kind=dp) :: tau_He_heth
    real(kind=dp) :: tau_H_heLya
    real(kind=dp) :: tau_He_heLya 
    real(kind=dp) :: tau_He2_he2th
    real(kind=dp) :: tau_He_he2th
    real(kind=dp) :: tau_H_he2th

    fHI = xHI*number_density_atom*(1.0_dp-abu_he)
    fHeI = xHeI*number_density_atom*abu_he
    fHeII = xHeII*number_density_atom*abu_he

    tau_H_heth  = fHI*sigma_H_heth ! opt depth of HI at HeI ion threshold
    tau_He_heth = fHeI*sigma_HeI_at_ion_freq ! opt depth of HeI at HeI ion threshold
    tau_H_heLya = fHI*sigma_H_heLya ! opt depth of H  at he+Lya (40.817eV)
    tau_He_heLya= fHeI*sigma_He_heLya ! opt depth of He at he+Lya (40.817eV) 
    tau_H_he2th = fHI*sigma_H_he2 ! opt depth of H at HeII ion threshold
    tau_He_he2th = fHeI*sigma_He_he2 ! opt depth of HeI at HeII ion threshold
    tau_He2_he2th = fHeII*sigma_HeII_at_ion_freq ! opt depth of HeII at HeII ion threshold
    
    ! Ratios of these optical depths needed in doric
    yfrac= tau_H_heth /(tau_H_heth +tau_He_heth)
    zfrac= tau_H_heLya/(tau_H_heLya+tau_He_heLya)
    y2afrac=  tau_He2_he2th /(tau_He2_he2th +tau_He_he2th+tau_H_he2th)
    y2bfrac=  tau_He_he2th /(tau_He2_he2th +tau_He_he2th+tau_H_he2th)
    
  end subroutine prepare_doric_factors

  ! =======================================================================

  !> Calculates the column density (of hydrogen)
  !! for a cell of ionization fraction xh, length path, 
  !! and density ndenstime dependent ionization state for hydrogen
  function coldens(path,neufrac,ndens,abundance)

    ! Calculates the column density (of hydrogen)
    ! for a cell of ionization fraction xh, length dr, 
    ! and density ndens

    real(kind=dp) :: coldens
    real(kind=dp),intent(in) :: path     !< path length through cell
    real(kind=dp),intent(in) :: neufrac  !< neutral H/He fraction
    real(kind=dp),intent(in) :: ndens    !< total number density   
    real(kind=dp),intent(in) :: abundance !<the abundance of the element
    ! Column density over a distance dr (cell)
    coldens=neufrac*ndens*path*abundance

  end function coldens

  !=======================================================================

  !> Sets the boundary condition for the hydrogen column density
  function coldens_bndry_HI()

    ! Sets the boundary condition for the hydrogen column density

    use cgsphotoconstants, only: sigma_HI_at_ion_freq,sigma_HeI_at_ion_freq,sigma_HeII_at_ion_freq
    use parameter, only: boundary_tauHI

    real(kind=dp):: coldens_bndry_HI

    ! Column density at the left of the first cell
    coldens_bndry_HI=boundary_tauHI/sigma_HI_at_ion_freq

  end function coldens_bndry_HI


  function coldens_bndry_HeI()
    use cgsphotoconstants, only: sigma_HeI_at_ion_freq
    use parameter, only: boundary_tauHeI
    !use radiation_sizes, only: boundary_tauHeI
    real(kind=dp):: coldens_bndry_HeI  
    coldens_bndry_HeI=boundary_tauHeI/sigma_HeI_at_ion_freq
  end function coldens_bndry_HeI

  function coldens_bndry_HeII()
    use cgsphotoconstants, only: sigma_HeII_at_ion_freq
    use parameter, only: boundary_tauHeII
    !use radiation_sizes, only: boundary_tauHeII
    real(kind=dp):: coldens_bndry_HeII  
    coldens_bndry_HeII=boundary_tauHeII/sigma_HeII_at_ion_freq
  end function coldens_bndry_HeII


end module doric_module
  
