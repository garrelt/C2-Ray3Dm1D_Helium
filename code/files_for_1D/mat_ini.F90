!>
!! \brief This module contains data and routines for handling 
!! the material properties on the grid (1D)
!!
!! These properties are; density, temperature, clumping, ionization fractions
!! 
!! \b Author: Garrelt Mellema
!!
!! \b Date: 
!!
!! \b Version:  1D test problems\n
!! Problem 1: constant density (Strömgren problem)\n
!! Problem 2: 1/r density\n
!! Problem 3: 1/r^2 density (with core of radius r_core)\n
!! Problem 4: cosmological constant density (Shapiro & Giroux problem)\n


module material

  ! This module contains the grid data and routines for initializing them.
  ! These are
  !  - mat_ini : initializes temperature and ionization fractions at start

  ! Version: 1D test problems

  ! Problem 1: constant density (Strömgren problem)
  ! Problem 2: 1/r density
  ! Problem 3: 1/r^2 density (with core of radius r_core)
  ! Problem 4: cosmological constant density (Shapiro & Giroux problem)

  use precision, only: dp,si
  use cgsconstants, only: bh00,bhe00,bhe10,ini_rec_colion_factors,m_p
  use cgsconstants, only: brech0,breche0,breche1
  use astroconstants, only: YEAR
  use sizes, only: mesh
  use file_admin, only: stdinput, file_input
  use my_mpi
  use grid, only: x,vol
  use c2ray_parameters, only: epsilon
  use abundances, only: abu_h,abu_he,mu
  use cosmology_parameters, only: Omega0, H0
  !use cosmology, only: cosmology_init,H0,t0,zred_t0

  implicit none

  ! ndens - number density (cm^-3) of a cell
  ! temper - temperature (K) of a cell
  ! xh - ionization fractions for one cell
  real(kind=dp) :: ndens(mesh(1),1,1) !< number density (cm^-3) of a cell
  real(kind=dp) :: temper(mesh(1)) !< temperature (K) of a cell
  real(kind=dp) :: xh(mesh(1),0:1) !< ionization fractions for one cell
  real(kind=dp) :: xhe(mesh(1),0:2)!< ionization fraction He for one cell
  real(kind=dp) :: clumping !< global clumping factor
  real(kind=dp) :: r_core !< core radius (for problems 2 and 3) 
  real(kind=dp) :: dens_core !< core density (for problems 2 and 3)
  integer :: testnum !< number of test problem (1 to 4)
  logical :: isothermal !< is the run isothermal?
  real(kind=dp),dimension(3) :: gamma_uvb !< UV background for HI, HeI, HeII
  ! needed for analytical solution of cosmological Ifront
  real(kind=dp) :: t1 !< parameter for analytical solution of test 4 
  real(kind=dp) :: t0_t !< parameter for analytical solution of test 4 
  real(kind=dp) :: eta !< parameter for analytical solution of test 4  
  real(kind=dp) :: zred00 
  real(kind=dp),public :: n_LLS  ! just because cosmology needs it in the 3D version
  ! b) Model Songaila & Cowie (2010)
  real(kind=dp) :: y_LLS 


!*TEST******************************************************
   type ionstates    
     real(kind=dp) :: h(0:1)          !< H  ionization fractions        
     real(kind=dp) :: he(0:2)         !< He ionization fractions           
     real(kind=dp) :: h_av(0:1)       !< average H  ionization fractions        
     real(kind=dp) :: he_av(0:2)      !< average He ionization fractions 
     real(kind=dp) :: h_old(0:1)      !< H  ionization fractions from last time step
     real(kind=dp) :: he_old(0:2)     !< He ionization fractions from last time step
  end type ionstates
 


!***********************************************************

#ifdef MPI
  integer,private :: ierror !< MPI error flag
#endif

contains

  ! ============================================================================

  !> Initializes material properties on grid\n
  !! \b Author: Garrelt Mellema\n
  !! \b Date: 20-Aug-2006 (f77 21-May-2005 (derives from mat_ini_cosmo2.f))\n
  !! \b Version:
  !! - 1D\n
  !! - Four different test problems\n
  !! - Initially completely neutral\n

  subroutine mat_ini (restart)

    ! Initializes material properties on grid

    ! Author: Garrelt Mellema

    ! Date: 20-Aug-2006 (f77 21-May-2005 (derives from mat_ini_cosmo2.f))

    ! Version: 
    ! - 1D
    ! - Four different test problems
    ! - Initially completely neutral

    integer,intent(out) :: restart !< will be /= 0 if a restart is intended

    integer :: i,n ! loop counters
    real(kind=dp) :: dens_val
    real(kind=dp) :: temper_val
    real(kind=dp) :: alpha
    character(len=1) :: answer

    real(kind=dp),dimension(3) :: xions

    type(ionstates) :: ion

    ! restart
    restart=0 ! no restart by default

    ! Ask for input
    if (rank == 0) then
       if (.not.file_input) then
          write(*,'(A,$)') 'Which test? (1-4): '
       endif
       read(stdinput,*) testnum

#ifdef MPI
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(testnum,1,MPI_INTEGER,0,MPI_COMM_NEW,ierror)
#endif
    endif

    ! Set alpha according to test problem
    select case (testnum)
    case(1,4) 
       alpha=0.0
    case(2) 
       alpha=-1.0
    case(3) 
       alpha=-2.0
    end select

    if (rank == 0) then
       if (testnum == 1 .or. testnum == 4) then
          if (.not.file_input) write(*,'(A,$)') 'Enter density (cm^-3): '
          read(stdinput,*) dens_val
       elseif (testnum == 2 .or. testnum == 3) then
          if (.not.file_input) write(*,'(A,$)') 'Enter reference (core) radius (cm): '
          read(stdinput,*) r_core
          if (.not.file_input) write(*,'(A,$)') 'Enter density at reference (core)', &
               ' radius(cm^-3): '
          read(stdinput,*) dens_val
       endif
       
       if (.not.file_input) write(*,'(A,$)') 'Enter clumping factor: '
       read(stdinput,*) clumping
       if (.not.file_input) write(*,'(A,$)') 'Enter initial temperature (K): '
       read(stdinput,*) temper_val
       if (.not.file_input) write(*,'(A,$)') 'Isothermal? (y/n): '
       read(stdinput,*) answer
       ! Isothermal?
       if (answer == 'y' .or. answer == 'Y') then
          isothermal=.true.
       else
          isothermal=.false.
       endif
       if (.not.file_input) write(*,'(A,$)') 'Ionizing background (HI,HeI,HeII) (s^-1): '
       read(stdinput,*) gamma_uvb
       call ini_rec_colion_factors(temper_val) !initialize the collisional ion and recomb rates for inital temp

    endif
#ifdef MPI
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(dens_val,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    if (testnum == 2.or.testnum == 3) &
         call MPI_BCAST(r_core,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(clumping,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(temper_val,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(isothermal,1,MPI_LOGICAL,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(gamma_uvb,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
#endif

    ! For test problem 4: cosmological parameters
    if (testnum == 4) then
       ! Ask for cosmological parameters
       if (rank == 0) then
          write(*,'(A,$)') 'Initial redshift?'
          read(stdinput,*) zred00
          write(*,*) 'redshift=', zred00
       endif
#ifdef MPI
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(zred00,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
#endif

!       call cosmology_init (.true.)
!    else
!       call cosmology_init (.false.)
    endif

    ! Assign density and temperature to grid
    
    select case (testnum)
    case(1)
       do i=1,mesh(1)
          ndens(i,1,1)=dens_val
          temper(i)=temper_val
       enddo
       
    case(2,3)
       dens_core=dens_val
       do i=1,mesh(1)
          !     This is an attempt to make the initial conditions more
          !     correct: the density of a cell is the integrated density
          !     (mass) divided by the volume. Seems to give a worse fit
          !     to the analytical solution (slightly)
          !              rl=r(i)-0.5*dr
          !              rr=r(i)+0.5*dr
          !              ndens(i)=4.0*pi*dens_val*r_core**(-alpha)/vol(i)*
          !     $            (rr**(3.0+alpha)-rl**(3.0+alpha))/(3.0+alpha)
          
          !     This is just a straight sampling of the density distribution
          !     using the value at the cell centre.
          if (testnum == 3.and.x(i) <= r_core) then
             ! Flat core for test 3
             ndens(i,1,1)=dens_val
          else
             ndens(i,1,1)=dens_val*(x(i)/r_core)**alpha
          endif
          temper(i)=temper_val
       enddo
       
    case(4)
       ! For cosmological simulations, mean IGM
       dens_core=dens_val    ! save initial density value
       ! Parameters needed for the analytical solution
       ! Recombination time at z0
       t1 = 1./(bh00*clumping*dens_core) 
       t0_t = 2.0_dp*(1.+zred00)**(-1.5)/(3.*H0*sqrt(Omega0))
       eta = t0_t/t1*(1.+zred00)**3
       
       ! Set the density to the comoving value 
       ! (as is the spatial coordinate)
       ! evol_cosmo will set it to proper values.
       do i=1,mesh(1)
          ndens(i,1,1)=dens_val
          temper(i)=temper_val
       enddo
       dens_val=dens_val*(1.+zred00)**3 !otherwise recombination time 
       !scale below would be wrong
    end select
    
    ! Assign ionization fractions
    ! Use Gamma_UVB_H for this if it is not zero
    if (gamma_uvb(1) > 0.0) then
       do i=1,mesh(1)
          call find_ionfractions_from_uvb(i,ndens(i,1,1), xions)
          xh(i,0)=xions(1)
          xh(i,1)=1.0-xh(i,0)
          xhe(i,0)=xions(2)
          xhe(i,1)=xions(3)
          xhe(i,2)=1.0-(xhe(i,0)+xhe(i,1))
       enddo
    else
       do i=1,mesh(1)
          xh(i,0) =1.0_dp!1.0_dp-epsilon!-1.0e-14
          xh(i,1) =0.0_dp!epsilon!1.0e-14
          xhe(i,0)=1.0_dp-2.0_dp*epsilon!1.0_dp-epsilon!1.0_dp-2.0e-14    !1.0_dp-2.4e-9
          xhe(i,1)=epsilon!1.0e-14    !1.2e-9
          xhe(i,2)=epsilon!1.0e-14   !1.2e-9            
       enddo
    endif
    
    ! Report recombination time scale (in case of screen input)
    if (.not.file_input) write(*,'(A,1pe10.3,A)') 'Recombination time scale: ', &
         1.0/(dens_val*clumping*bh00*YEAR),' years'

    
  end subroutine mat_ini

  subroutine find_ionfractions_from_uvb (ii,nnd,xions)

    real(kind=dp),parameter :: convergence=0.01
    integer,intent(in) :: ii
    real(kind=dp),intent(in) :: nnd
    real(kind=dp),dimension(3),intent(out) :: xions

    real(kind=dp) :: rech2
    real(kind=dp) :: reche2
    real(kind=dp) :: reche3
    real(kind=dp) :: fe
    real(kind=dp) :: fe_prev

    rech2 = nnd * clumping * brech0
    reche2 = nnd * clumping* breche0
    reche3 = nnd * clumping* breche1
    fe=1.0
    
    ! Iterate to find the proper electron density (fe)
    do 
       xions(1)=fe*rech2/(gamma_uvb(1)+fe*rech2)
       xions(2)=fe*reche2/(gamma_uvb(2)*(1.0+gamma_uvb(3)/(fe*reche3)) + &
            fe*reche2)
       xions(3)=(1.0-xions(2))/(1.0+gamma_uvb(3)/(fe*reche3))
       fe_prev=fe
       fe=abu_h*(1.0-xions(1))+abu_he*(2.0-(2.0*xions(2)+xions(3)))
       if (ii == 1) then
          write(*,*) xions
          write(*,*) gamma_uvb(2:3)
          write(*,*) reche2,reche3
          write(*,*) fe_prev,fe
       endif
       if (abs(fe-fe_prev)/fe_prev < convergence) exit
    enddo

  end subroutine find_ionfractions_from_uvb

end module material
