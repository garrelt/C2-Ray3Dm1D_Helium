module AT_time_calculation
  
  use my_mpi
  use file_admin, only: logf
  use precision, only: dp
  use array, only: number_density_array,temperature_array
  use array, only: xHI_array,xHII_array
  use array, only: xHeI_array,xHeII_array,xHeIII_array
  use tped, only: electrondens
    use abundances, only: abu_he
  use array, only: source_Ifront_array, global_timestep_array, global_photoionization_HI_array,&
        global_photoionization_HeI_array,global_photoionization_HeII_array
  use parameter, only: mesh, time_partition
  use array, only: xfinal
  use astroconstants, only: YEAR
  use cgsconstants, only: bh00,bhe00, bhe10, albpow,alcpow, &
         colh0,colhe, temphe,temph0,ev2fr
  !use recombination_collision_factor, only: recombination_collision_factor_initialization
  use doric_module, only: prepare_doric_factors
  use clumping, only: clumping_factor
  use astroconstants, only: YEAR
  use parameter, only: do_AT_accurate_timestep

  implicit none

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AT_time_equation()
     
    implicit none

    real(kind=dp) :: dt
    real(kind=dp) :: x_end
    integer :: i,j,k
    real(kind=dp) :: partition,final_xHII

    real(kind=dp) :: ndens_electron
    real(kind=dp) :: yHI
    real(kind=dp) :: yHII
    real(kind=dp) :: yHeI
    real(kind=dp) :: yHeII
    real(kind=dp) :: yHeIII
    real(kind=dp) :: photoionization_HI_rate,photoionization_HeI_rate,photoionization_HeII_rate
    real(kind=dp) :: number_density,temperature
    real(kind=dp) :: bigdt, smalldt, averagedt
    real(kind=dp) :: y 
    real(kind=dp) :: z
    real(kind=dp) :: y2a 
    real(kind=dp) :: y2b 
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
    real(kind=dp) :: resultant_xHII
    real(kind=dp) :: range_xfinal
    real(kind=dp) :: prev_dt
    integer ::countcount
    integer :: I_internal,i_partition

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

   I_internal=1

   partition = 1.0/(time_partition+1.0)
   ! initial xHI = 1.0e-3
   !partition = 1.0e-3/(time_partition+1.0)

   do k = 1, mesh
     do j = 1, mesh
       do i = 1, mesh
         countcount = 0

         ! do timestep calculation for each photo-cells
         if (global_photoionization_HI_array(i,j,k).gt.0) then

           if (do_AT_accurate_timestep .eqv. .false.) then
            
             ! Fast method
             if (xHII_array(i,j,k).lt. (1.0-partition)) then
               final_xHII = xHII_array(i,j,k)+partition
               dt = (log((xHII_array(i,j,k)-1)/(final_xHII-1)))/global_photoionization_HI_array(i,j,k)
             else
               dt = 1.0e50_dp
             endif
             !<!<!< Hydrogen case
           else
             !<!<!< Helium case
             yHII = xHII_array(i,j,k)
             do while (yHII.gt.xfinal(I_internal).and.I_internal.le.size(xfinal))
               I_internal = I_internal+1
             enddo

             ! the target ionized fraction is a range
             ! range = [xfinal(I_internal), xfinal(I_internal)+range_xfinal] 
             !range_xfinal = 0.8/(size(xfinal)+1)
             range_xfinal = 1.0e-3*0.8/(size(xfinal)+1)
             if (I_internal.le.size(xfinal)) then

               yHI = xHI_array(i,j,k)
               yHeI = xHeI_array(i,j,k)
               yHeII = xHeII_array(i,j,k)
               yHeIII = xHeIII_array(i,j,k)
               temperature = temperature_array(i,j,k) 
               number_density = number_density_array(i,j,k)
               ndens_electron = electrondens(number_density,yHII, yHeII, yHeIII)
               photoionization_HI_rate = global_photoionization_HI_array(i,j,k)
               photoionization_HeI_rate = global_photoionization_HeI_array(i,j,k)
               photoionization_HeII_rate = global_photoionization_HeII_array(i,j,k)

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

               call prepare_doric_factors(yHI,yHeI,yHeII,number_density,y,z,y2a,y2b)

               smalldt = 0.0
               bigdt = 1/(arech0*number_density)
               averagedt = 0.5*(smalldt+bigdt)
               dt=averagedt 
               prev_dt=averagedt

               p = 0.96_dp
               l = 1.425_dp
               m = 0.737_dp 
               heliumfraction = abu_he/(1.0_dp-abu_he) 
               f = max(min(10.0_dp*yHI,1.0_dp),0.01_dp) 
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

               g_1=max(photoionization_HI_rate+ndens_electron*colli_HI,1.0e-200_dp)    
               g_2=max(photoionization_HeI_rate+ndens_electron*colli_HeI,1.0e-200_dp)
               A_32=max(photoionization_HeII_rate+ndens_electron*colli_HeII,1.0e-200_dp)

               A_11 = -(g_1+ndens_electron*alpha_HI_B)
               A_12 = (y*alpha_HeI_1+p*alpha_HeI_B)*ndens_electron*heliumfraction
               A_13 = ((f*z*(1.0_dp-v_factor) +v_factor*w)*alpha_HeII_B +alpha_HeII_2+ &
                    (1.0_dp-y2a-y2b)*alpha_HeII_1)*heliumfraction*ndens_electron  
               A_22 = -g_2-A_32-ndens_electron*(alpha_HeI_A-(1.0_dp-y)*alpha_HeI_1)
               A_33 = -ndens_electron*(alpha_HeII_A-y2a*alpha_HeII_1)            
               A_23 = -g_2+ndens_electron*alpha_HeII_B*(f*(1.0_dp-z)*(1.0_dp-v_factor)+ &
                    v_factor*(1.425_dp-w))-A_33+alpha_HeII_1*y2b*ndens_electron

               two_A_32 = 2.0_dp*A_32
               Bcoef = A_33-A_22
               Scoef = sqrt(Bcoef*Bcoef+4.0_dp*A_32*A_23)
               Kcoef = 1.0_dp/(A_23*A_32-A_33*A_22)
               Bminus = Bcoef-Scoef
               Bplus = Bcoef+Scoef

               lambda_1 = A_11
               lambda_2 = 0.5_dp*(A_33+A_22-Scoef)
               lambda_3 = 0.5_dp*(A_33+A_22+Scoef)

               xp_1 = -1.0_dp/A_11*(g_1+(A_12*A_33-A_13*A_32)*(g_2*Kcoef))   
               xp_2 = g_2*(A_33*Kcoef)                        
               xp_3 = -g_2*(A_32*Kcoef)  

               x2_1 = -A_13/(A_11-lambda_2)+(A_12/two_A_32)*Bplus/(A_11-lambda_2)  
               x2_2 = (-Bplus)/(two_A_32)
               x3_1 = (-two_A_32 *A_13+A_12*(Bminus))/(two_A_32*(A_11-lambda_3))  
               x3_2 = (-Bminus)/(two_A_32)

               Rcoef = two_A_32 *(xp_2-yHeII)
               Tcoef = xp_3-yHeIII

               c_1 = -xp_1+ (x3_1-x2_1)*(Rcoef/(2.0_dp*Scoef))+ &
                     Tcoef*((Bplus*x3_1/(2.0_dp*Scoef)-Bminus*x2_1/(2.0_dp*Scoef)))+ yHII
               c_2 = (Rcoef+(Bminus)*Tcoef)/(2.0_dp*Scoef)
               c_3 = -(Rcoef+(Bplus)*Tcoef)/(2.0_dp*Scoef)

               do 

                 ! arguments of exponential functions
                 lambda_1_dt = dt*lambda_1
                 lambda_2_dt = dt*lambda_2
                 lambda_3_dt = dt*lambda_3

                 ! exponential functions of homogeneous solutions
                 exp_lambda_1_dt = exp(lambda_1_dt)
                 exp_lambda_2_dt = exp(lambda_2_dt)
                 exp_lambda_3_dt = exp(lambda_3_dt  )

                 ! ionization fractions at the end of the time-step  
                 resultant_xHII = c_1*exp_lambda_1_dt+c_2*exp_lambda_2_dt*x2_1+c_3*exp_lambda_3_dt*x3_1+xp_1 

                 ! check if the time-step reaches the requirement
                 if (resultant_xHII .ge. xfinal(I_internal) .and. &
                     resultant_xHII .le. (xfinal(I_internal)+range_xfinal) ) then
                   exit
                 endif

                 countcount = countcount+1
  
                 if (resultant_xHII .lt. xfinal(I_internal) .and. countcount.eq.1) then
                   exit
                 endif

                 ! if the time-step is underestimated
                 if (resultant_xHII .le. xfinal(I_internal)) then

                   smalldt = averagedt
                   averagedt = 0.5*(smalldt+bigdt)
                   dt = averagedt

                 endif

                 ! if the time-step is overestimated
                 if (resultant_xHII .ge. (xfinal(I_internal)+range_xfinal)) then

                   bigdt = averagedt
                   averagedt = 0.5*(smalldt+bigdt)
                   dt = averagedt

                 endif    

               enddo

             else
               dt = 1.0e50_dp
             endif


           endif

           global_timestep_array(i,j,k) = dt

         endif
       enddo
     enddo
   enddo


  end subroutine AT_time_equation
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AT_find_adaptive_timestep(dt)

    implicit none

    real(kind=dp),intent(out) :: dt
    integer :: i,j,k
    integer :: pos_i,pos_j,pos_k

    dt = 1.1e40_dp
    do i = 1, mesh
      do j = 1, mesh
        do k = 1, mesh
          if ( (global_timestep_array(i,j,k).gt.0) .and. (global_timestep_array(i,j,k).lt.dt) )then

             dt = global_timestep_array(i,j,k)

          endif
        enddo
      enddo
    enddo

  end subroutine AT_find_adaptive_timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module AT_time_calculation
  
