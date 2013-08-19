!>
!! \brief This module contains routines for calculating the ionization evolution of a single grid point.
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema, Martina M. Friedrich
!!
!! \b Date: 
!!
!!

module doric_module

  ! This module contains routines having to do with the calculation of
  ! the ionization evolution of a single grid point.
  ! It can used for Yguazu-a or non-hydro photo-ionization calculations.

  ! - doric : time dependent solution of the ionization equation
  ! - coldens : column density of a single cell.
  ! - coldensh_bndry : Boundary condition for H column density

  use precision, only: dp
  use abundances  

  implicit none

  ! No public data associated with this module

contains

  !=======================================================================

  !> calculates time dependent ionization state for hydrogen

  subroutine doric (dt,rhe,rhh,ion,phi,yfrac,zfrac,y2afrac,y2bfrac)

    ! calculates time dependent ionization state for hydrogen and helium

    ! Author: Garrelt Mellema, Martina M. Friedrich (helium)
    ! Date: 11-May-2005 (f90) (15 June 2004

    ! Version:
    ! Reimplemented helium, now in the correct way with two different 
    ! time evolution factors for the three helium states

    ! Modifications:
    ! 15-Jun-2004 (GM): test for (1.0d0-ee)/deltht to avoid fp errors.
    ! 11-May-2005 (GM): f90

    use cgsconstants, only: bh00,bhe00, bhe10, albpow,alcpow, &
         colh0,colhe, temphe,temph0,ev2fr,brech0,areche0,breche0, &
         oreche0,breche1,treche1,n_el_crit, &
         colli_HI,colli_HeI,colli_HeII, ini_rec_colion_factors,arech0,areche1,v
    use cgsphotoconstants, only: sigma_HI_at_ion_freq, sigma_HeI_at_ion_freq
    !use subgrid_clumping, only: clumping
    use file_admin, only: stdinput
    use material, only: clumping, ionstates
    use abundances, only: abu_he
    use radiation, only: photrates
    use tped, only: electrondens ! should this really be used inside doric?
    use c2ray_parameters, only: epsilon
    real(kind=dp),intent(in)    :: dt     !< time step
    real(kind=dp),intent(inout) :: rhe    !< electron density
    real(kind=dp),intent(in)    :: rhh    !< H density (or total density?) TOTAL
    !    real(kind=dp),intent(in)    :: coldensh0   	 !< H0 column density of the cell
    !    real(kind=dp),intent(in)    :: coldenshe0 	 !< He0 column density of the cell
    real(kind=dp),intent(inout)   :: yfrac
    real(kind=dp),intent(inout)   :: zfrac   
    real(kind=dp),intent(inout) :: y2afrac
    real(kind=dp),intent(inout) :: y2bfrac
    !   logical, intent(in) :: flag
    type(photrates),intent(in)  :: phi !< photo-ionization rates
    type(ionstates) :: ion !< ionzation states

    real(kind=dp) ::alpha_h_B,alpha_he_1,alpha_he_B,alpha_he_A,alpha_he2_B, &
         alpha_he2_2,alpha_h_A,alpha_he2_A,alpha_he2_1 
    real(kind=dp) :: pfrac,ffrac,wfrac
    real(kind=dp) :: lambda1, lambda2, lambda3
    real(kind=dp) :: lam1dt, lam2dt, lam3dt
    real(kind=dp) :: elam1dt,elam2dt,elam3dt
    real(kind=dp) :: eigv2x, eigv2y, eigv3x, eigv3y
    real(kind=dp) :: coef1,coef2,coef3,Scoef,Rcoef,Tcoef, rx,ry,rz,QHEPcoef, &
         Bcoef,BminusS,BplusS
    real(kind=dp) :: avg_factor_1, avg_factor_2, avg_factor_3
    real(kind=dp) :: Lmat,Mmat,Nmat,Pmat,Qmat,Emat,twoaihe1 
    real(kind=dp) :: aih0, aihe0, aihe1
    real(kind=dp) :: heliumfraction,  normfac ! calculate as abu_he/(1-abu_he)
    !***********************************************************

    pfrac= 0.96_dp ! see Osterbrock, 1989
    heliumfraction=abu_he/(1.0_dp-abu_he)    ! He/H
    ffrac= max(min(10.0_dp*ion%h(0),1.0_dp),0.01_dp) 
    !vfrac= 0.285_dp 		           ! taken from Hummer and Seaton, 1964 table five
    ! This is a function of temperature, moved to
    ! ini_rec_colion_factors in cgsconstants
    wfrac= (1.425_dp-0.737_dp)+0.737_dp*yfrac  ! These numbers come from Flower & Perinotto (1980)

    alpha_h_B   = clumping* brech0
    alpha_h_A   = clumping* arech0
    alpha_he_1  = clumping* oreche0
    alpha_he_B  = clumping* breche0    
    alpha_he_A  = clumping* areche0
    alpha_he2_B = clumping* breche1
    alpha_he2_A = clumping* areche1
    alpha_he2_2 = clumping* treche1
    alpha_he2_1 =  alpha_he2_A - alpha_he2_B

    aih0=max(phi%photo_cell_HI+rhe*colli_HI,1.0e-200_dp)       ! colli_HI and co are calculated in a routine 
    aihe0=max(phi%photo_cell_HeI+rhe*colli_HeI,1.0e-200_dp) ! inside cgsconstants. ini_rec_colion_factors   
    aihe1=max(phi%photo_cell_HeII+rhe*colli_HeII,1.0e-200_dp) ! that should be called every time temp0 is updated


    !***********************************************************  
    !         L M N           A
    ! Matrix= 0 P Q        g= G
    !         0 H E           0
    !
    ! solve differential Eq: dx/dt= Matrix*x + g
    !
    ! where x is the vector:  x_h(1), x_he(1), x_he(2)
    ! alpha_he2_2=0.0_dp;

    Lmat= -(aih0+rhe*alpha_h_B)
    Mmat=  ( yfrac*rhe*alpha_he_1  + pfrac*rhe*alpha_he_B )*heliumfraction
    Nmat=  ((ffrac*zfrac*(1.0_dp-v) +v*wfrac)*alpha_he2_B +alpha_he2_2 &
         +(1.0_dp-y2afrac-y2bfrac)*alpha_he2_1)*heliumfraction*rhe  
    Pmat=    - aihe0-aihe1-rhe*(alpha_he_A-(1.0_dp-yfrac)*alpha_he_1)
    Emat=  -rhe*(alpha_he2_A-y2afrac*alpha_he2_1)            
    Qmat=  -aihe0+rhe*alpha_he2_B*(ffrac*(1.0_dp-zfrac)*(1.0_dp-v)+ &
         v*(1.425_dp-wfrac))-Emat+alpha_he2_1*y2bfrac*rhe

    !! Hmat=  aihe1  !write aihe1 everywhere instead of Hmat 
    !! Amat=  aih0   !write aih0 everywhere instead of Amat
    !! Gmat=  aihe0  !write aihe0 everywhere instead of Gmat



    !*********************************************************************
    !THIS IS FOR TEST WITH A-coefficients only
    ! Lmat= -(aih0+rhe*alpha_h_A)
    ! Mmat=  0.0_dp
    ! Nmat=  0.0_dp 
    ! Pmat=  -(aihe0+aihe1+rhe*alpha_he_A)
    ! Emat=  -rhe*alpha_he2_A            
    ! Qmat=  -aihe0-Emat

    !THIS IS FOR TEST WITH B-coefficients only
    ! Lmat= -(aih0+rhe*alpha_h_B)
    ! Mmat=  0.0_dp
    ! Nmat=  0.0_dp 
    ! Pmat=  -(aihe0+aihe1+rhe*alpha_he_B)
    ! Emat=  -rhe*alpha_he2_B            
    ! Qmat=  -aihe0-Emat
    !*********************************************************************

    ! These are just abbrivations 
    Bcoef= Emat-Pmat

    Scoef   = sqrt(Bcoef*Bcoef+4.0_dp*aihe1*Qmat)
    QHEPcoef=1.0_dp/(Qmat*aihe1-Emat*Pmat)

    BminusS = Bcoef-Scoef
    BplusS  = Bcoef+Scoef


    ! These are the Eigenvalues of the Matrix
    lambda1 = Lmat
    lambda2 = 0.5_dp*(Emat+Pmat-Scoef)
    lambda3 = 0.5_dp*(Emat+Pmat+Scoef)


    ! This is the particular solution vector, which one gets with the Ansatz:
    ! dx/dt=0

    rx = -1.0_dp/Lmat*(aih0+(Mmat*Emat-Nmat*aihe1)*(aihe0*QHEPcoef))   
    ry = aihe0*(Emat*QHEPcoef)                        
    rz = -aihe0*(aihe1*QHEPcoef)                      


    ! lambda1, lambda2 and lambda3 are the Eigenvalues.
    ! The coerresponding Eigenvectors are: 
    ! [1 0 0], [eigv2x eigv2y 1] and [eigc3x eigv3y 1] 

    twoaihe1= 2.0_dp*aihe1

    eigv2x = -Nmat/(Lmat-lambda2)+(Mmat/twoaihe1)*BplusS/(Lmat-lambda2)  
    eigv3x = (-twoaihe1*Nmat+Mmat*(BminusS))/(twoaihe1*(Lmat-lambda3))  
    eigv2y = (-BplusS)/(twoaihe1)
    eigv3y = (-BminusS)/(twoaihe1)

    !These are just abbrivations
    Rcoef=twoaihe1*(ry-ion%he_old(1))
    Tcoef=rz-ion%he_old(2)

    !These are the coefficients for the general solution: 
    !x=coef1*Eigvec1*exp(t*lambda1)+
    !  coef2*Eigvec2*exp(t*lambda2)+
    !  coef3*Eigvec3*exp(t*lambda3)+
    !  [rx,ry,rz]

    !------------------------------------------------------------
    coef2 =  (Rcoef+(BminusS)*Tcoef)  /(2.0_dp*Scoef)
    coef3 = -(Rcoef+(BplusS)*Tcoef)  /(2.0_dp*Scoef)

    ! coef1 is rather complex, you can group the different terms in different ways. 
    ! You should try to keep the number small, so divide if possible number of equal size. 
    ! I found the formulation below to work in all the cases I tested. 
    ! Best option: 
    coef1 = -rx+ (eigv3x-eigv2x)*(Rcoef/(2.0_dp*Scoef))+ &
         Tcoef*((BplusS*eigv3x/(2.0_dp*Scoef)-BminusS*eigv2x/(2.0_dp*Scoef)))+ &
         ion%h_old(1)

    lam1dt = dt*lambda1
    lam2dt = dt*lambda2
    lam3dt = dt*lambda3

    elam1dt = exp(lam1dt)
    elam2dt = exp(lam2dt)
    elam3dt = exp(lam3dt)

    ion%h(1) = coef1*elam1dt+coef2*elam2dt*eigv2x+coef3*elam3dt*eigv3x + rx 
    ion%he(1)=               coef2*elam2dt*eigv2y+coef3*elam3dt*eigv3y + ry
    ion%he(2)=               coef2*elam2dt       +coef3*elam3dt        + rz
    ion%h(0) = 1.0_dp-ion%h(1)
    ion%he(0)= 1.0_dp-ion%he(1)-ion%he(2)

    !-----------------------------------------------------------

    ! determine neutral densities (take care of precision fluctuations)

    if (ion%h(0) < epsilon) then
       ion%h(0)=epsilon
       ion%h(1)=1.0_dp-epsilon
    endif
    if (ion%h(1) < epsilon) then
       ion%h(1) =epsilon
       ion%h(0)=1.0_dp-epsilon
    endif

    if ((ion%he(0).le.epsilon).or.(ion%he(1).le.epsilon).or. &
         (ion%he(2).le.epsilon))then
       if (ion%he(0) < epsilon) then
          ion%he(0)=epsilon
       endif

       if (ion%he(1) < epsilon) then
          ion%he(1)=epsilon      
       endif
       if (ion%he(2) < epsilon) then
          ion%he(2)=epsilon
       endif

       normfac=ion%he(0)+ion%he(1)+ion%he(2)
       ion%he(0)=ion%he(0)/normfac
       ion%he(1)=ion%he(1)/normfac
       ion%he(2)=ion%he(2)/normfac
    endif
    !------------------------------------------------------------

    !*****AVERAGE*******
    ! Determine average ionization fraction over the time step
    ! Mind fp fluctuations. (1.0-ee)/deltht should go to 1.0 for
    ! small deltht, but finite precision leads to values slightly
    ! above 1.0 and for very small values even to 0.0.

    if (abs(lam1dt).lt.1.0e-8) then
       avg_factor_1=coef1
    else
       avg_factor_1 = coef1*(elam1dt- 1.0_dp)/lam1dt
    endif

    if (abs(lam2dt).lt.1.0e-8) then
       avg_factor_2=coef2
    else
       avg_factor_2 = coef2*(elam2dt- 1.0_dp)/lam2dt       
    endif

    if (abs(lam3dt).lt.1.0e-8) then
       avg_factor_3=coef3
    else
       avg_factor_3 = coef3*(elam3dt- 1.0_dp)/lam3dt
    endif

    ion%h_av(1)  = rx + avg_factor_1 + eigv2x*avg_factor_2 + eigv3x*avg_factor_3
    ion%he_av(1) = ry +                eigv2y*avg_factor_2 + eigv3y*avg_factor_3
    ion%he_av(2) = rz +                       avg_factor_2 +        avg_factor_3
    ion%h_av(0)  = 1.0_dp-ion%h_av(1)
    ion%he_av(0) = 1.0_dp-ion%he_av(1)-ion%he_av(2)

    if (ion%h_av(1) < epsilon) then
       ion%h_av(1)=epsilon
       ion%h_av(0)=1.0_dp-epsilon
    endif
    if (ion%h_av(0) < epsilon) then
       ion%h_av(0)=epsilon
       ion%h_av(1)=1.0_dp-epsilon
    endif

    if ((ion%he_av(0).le.epsilon).or.(ion%he_av(1).le.epsilon).or. &
         (ion%he_av(2).le.epsilon))then

       if (ion%he_av(1) < epsilon) ion%he_av(1)=epsilon   
       if (ion%he_av(2) < epsilon) ion%he_av(2)=epsilon
       if (ion%he_av(0) < epsilon) ion%he_av(0)=epsilon
       normfac=ion%he_av(0)+ion%he_av(1)+ion%he_av(2)
       ion%he_av(0)=ion%he_av(0)/normfac
       ion%he_av(1)=ion%he_av(1)/normfac
       ion%he_av(2)=ion%he_av(2)/normfac
    endif

    return
  end subroutine doric

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
    use radiation, only: boundary_tauHI

    real(kind=dp):: coldens_bndry_HI

    ! Column density at the left of the first cell
    coldens_bndry_HI=boundary_tauHI/sigma_HI_at_ion_freq

  end function coldens_bndry_HI


  function coldens_bndry_HeI()
    use cgsphotoconstants, only: sigma_HeI_at_ion_freq
    use radiation, only: boundary_tauHeI
    real(kind=dp):: coldens_bndry_HeI  
    coldens_bndry_HeI=boundary_tauHeI/sigma_HeI_at_ion_freq
  end function coldens_bndry_HeI

  function coldens_bndry_HeII()
    use cgsphotoconstants, only: sigma_HeII_at_ion_freq
    use radiation, only: boundary_tauHeII
    real(kind=dp):: coldens_bndry_HeII  
    coldens_bndry_HeII=boundary_tauHeII/sigma_HeII_at_ion_freq
  end function coldens_bndry_HeII


end module doric_module
  
