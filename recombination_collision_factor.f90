module recombination_collision_factor

  use precision, only: dp
  use cgsconstants, only: arech0,brech0,areche0,breche0,oreche0,areche1,breche1,treche1,&
                          colh0,temph0,colli_HI,colli_HeI,colli_HeII,v_factor,temphe,colhe

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine recombination_collision_factor_initialization(temperature)
     
    implicit none

    real(kind=dp),intent(in) :: temperature
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

  end subroutine recombination_collision_factor_initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module recombination_collision_factor
