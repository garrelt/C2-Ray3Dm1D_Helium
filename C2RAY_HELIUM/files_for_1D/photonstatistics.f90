
MODULE photonstatistics

  
  ! This module handles the calculation of the photon statistics
  ! For: C2-Ray

  ! Author: Garrelt Mellema

  ! Date: 26-Sep-2006

  ! Photon statistics
  ! photon_loss is a sum over all sources, summing is done in evolve0d.
  ! For parallelization consider a reduction, or may also be dropped
  ! entirely.
  ! Photon_losst contains the threaded values

  USE PRECISION, ONLY: dp
  USE cgsconstants, ONLY: albpow,bh00,colh0,temph0
  USE cgsconstants, ONLY: alcpow,bhe00,bhe10,colhe,temphe
  USE sizes, ONLY: mesh
  USE grid, ONLY: vol
  USE abundances, ONLY: abu_he
  USE material, ONLY: ndens, xh, xhe, temper, clumping
  USE tped, ONLY: electrondens
  !use subgrid_clumping, only: clumping

  LOGICAL,PARAMETER :: do_photonstatistics=.TRUE.
  REAL(kind=dp) :: totrec
  REAL(kind=dp) :: totcollisions
  REAL(kind=dp) :: dh0,dhe0,dhe1,dhe2, dh1
  REAL(kind=dp) :: total_ion
  REAL(kind=dp) :: grtotal_ion
  REAL(kind=dp) :: photon_loss

  REAL(kind=dp),PRIVATE :: h0_before,h0_after,h1_before,h1_after 
  REAL(kind=dp),PRIVATE :: he0_before,he0_after,he1_before,he1_after
  REAL(kind=dp),PRIVATE :: he2_before,he2_after
  INTEGER,PRIVATE :: i,j,k



CONTAINS
  SUBROUTINE initialize_photonstatistics ()

    ! set total number of ionizing photons used to zero
    grtotal_ion=0.0

  END SUBROUTINE initialize_photonstatistics


  SUBROUTINE calculate_photon_statistics (dt)



    REAL(kind=dp),INTENT(in) :: dt


    ! Call the individual routines needed for this calculation

    CALL state_after () ! number of neutrals after integration
    CALL total_rates (dt) ! total photons used in balancing recombinations etc.
    CALL total_ionizations () ! final statistics
    
  END SUBROUTINE calculate_photon_statistics

  SUBROUTINE state_before ()


    ! Photon statistics: calculate the number of neutrals before integration
    h0_before=0.0
    h1_before=0.0
    he0_before=0.0
    he1_before=0.0
    he2_before=0.0
    DO i=1,mesh(1)
       h0_before=h0_before+vol(i)*ndens(i,1,1)*xh(i,0)*(1.0_dp-abu_he)
       h1_before=h1_before+vol(i)*ndens(i,1,1)*xh(i,1)*(1.0_dp-abu_he)
       he0_before=he0_before+vol(i)*ndens(i,1,1)*xhe(i,0)*abu_he
       he1_before=he1_before+vol(i)*ndens(i,1,1)*xhe(i,1)*abu_he
       he2_before=he2_before+vol(i)*ndens(i,1,1)*xhe(i,2)*abu_he

    ENDDO
    
  END SUBROUTINE state_before


  SUBROUTINE total_rates(dt)



    REAL(kind=dp),INTENT(in) :: dt


    REAL(kind=dp),DIMENSION(0:1) :: yh
    REAL(kind=dp),DIMENSION(0:2) :: yhe
 
    ! Photon statistics: Determine total number of recombinations/collisions
    ! Should match the code in doric_module

    totrec=0.0
    totcollisions=0.0
    DO i=1,mesh(1)
       yh(0)=xh(i,0)
       yh(1)=xh(i,1)
       yhe(0)=xhe(i,0)
       yhe(1)=xhe(i,1)
       yhe(2)=xhe(i,2)
       totrec=totrec+vol(i)*ndens(i,1,1)* electrondens(ndens(i,1,1),yh,yhe)*clumping* &
          (xh(i,1)*(1.0_dp-abu_he)*  & 
          1.0_dp/(1.0_dp/(bh00*(temper(i)/1e4)**albpow)+1.0_dp/(bh00*5.0_dp*(temper(i)/1e4)**(1.95_dp*albpow))) +&
	  xhe(i,1)*abu_he*  &
          (bhe00*(temper(i)/1e4)**alcpow) + &
	  xhe(i,2)*abu_he*  &
          1.0_dp/(1.0_dp/(bhe10*(temper(i)/1e4)**(0.95_dp*albpow))+1.0_dp/(bhe10*11.0_dp*(temper(i)/1e4)**(albpow*1.95_dp))))


       totcollisions=totcollisions+vol(i)*ndens(i,1,1)*(1.0_dp-abu_he)* &
            xh(i,0)*electrondens(ndens(i,1,1),yh,yhe)* &
            colh0*SQRT(temper(i))*EXP(-temph0/temper(i))+ &
	    vol(i)*ndens(i,1,1)*abu_he*   &
            xhe(i,0)*electrondens(ndens(i,1,1),yh,yhe)* &
            colhe(0)*SQRT(temper(i))*EXP(-temphe(0)/temper(i))+ &
	    vol(i)*ndens(i,1,1)*abu_he*   &
            xhe(i,1)*electrondens(ndens(i,1,1),yh,yhe)* &
            colhe(1)*SQRT(temper(i))*EXP(-temphe(1)/temper(i))
    ENDDO

    totrec=totrec*dt
    totcollisions=totcollisions*dt

  END SUBROUTINE total_rates
  

  SUBROUTINE state_after()

    
    ! Photon statistics: Calculate the number of neutrals after the integration
    h0_after=0.0
    h1_after=0.0
    he0_after=0.0
    he1_after=0.0
    he2_after=0.0

    DO i=1,mesh(1)
       h0_after=h0_after+vol(i)*ndens(i,1,1)*xh(i,0)*(1.0_dp-abu_he)
       h1_after=h1_after+vol(i)*ndens(i,1,1)*xh(i,1)*(1.0_dp-abu_he)
       he0_after=he0_after+vol(i)*ndens(i,1,1)*xhe(i,0)*abu_he
       he1_after=he1_after+vol(i)*ndens(i,1,1)*xhe(i,1)*abu_he
       he2_after=he2_after+vol(i)*ndens(i,1,1)*xhe(i,2)*abu_he
    ENDDO
    
  END SUBROUTINE state_after
  

  SUBROUTINE total_ionizations ()

    
    ! Photon statistics: Total number of new ionizations
    dh0=(h0_before-h0_after)
    dh1=(h1_before-h1_after)
    dhe0=(he0_before-he0_after)
    dhe1=(he1_before-he1_after)
    dhe2=(he2_before-he2_after)
    total_ion=totrec+dh0+dhe0+dhe1
    
  END SUBROUTINE total_ionizations

END MODULE photonstatistics
