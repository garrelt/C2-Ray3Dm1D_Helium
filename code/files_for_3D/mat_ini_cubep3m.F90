module material

  ! This module contains the grid data and routines for initializing them.
  ! These are
  !  - mat_ini : initializes temperature and ionization fractions at start
  !  - dens_ini : initializes the density field (from PMFAST output)
  !  - xfrac_ini : initializes ionization fractions (in case of restart).

  use precision, only: dp,si
  use sizes, only: mesh
  use file_admin, only: stdinput, logf, results_dir, file_input
  use my_mpi
  use grid, only: dr,vol,sim_volume
  use cgsconstants, only: m_p,ini_rec_colion_factors, c
  use cgsphotoconstants, only: sigma_HI_at_ion_freq
  use astroconstants, only: M_solar, Mpc
  use cosmology_parameters, only: Omega_B, Omega0, rho_crit_0, h, H0
  use nbody, only: nbody_type, M_grid, M_particle, id_str, dir_dens, NumZred, Zred_array, dir_clump, dir_LLS
  use nbody, only: densityformat, densityaccess, densityheader
  use nbody, only: clumpingformat, clumpingaccess, clumpingheader
  use nbody, only: LLSformat, LLSaccess, LLSheader  
  use nbody, only: density_convert_particle, density_convert_grid, density_unit
  use abundances, only: mu, abu_he
  use c2ray_parameters, only: type_of_clumping,clumping_factor, cosmological, type_of_LLS

  implicit none

  ! ndens - number density (cm^-3) of a cell
  ! SINGLE PRECISION! Be careful when passing this as argument to
  ! functions and subroutines.
  ! temper - temperature of grid cell under consideration
  ! temperature_grid - temperature (K) of entire grid (only used if not isothermal)
  ! temper_val - initial temperature and the one used if isothermal  
  real(kind=dp) :: temper
  real(kind=dp) :: temper_val
  real(kind=si),dimension(:,:,:,:),allocatable :: temperature_grid
  ! xh - ionization fractions for one cell
  real(kind=dp),dimension(:,:,:),allocatable :: ndens
  real(kind=dp),dimension(:,:,:,:),allocatable :: xh
  real(kind=dp),dimension(:,:,:,:),allocatable :: xhe
  logical isothermal
  real,public :: clumping
  real,dimension(:,:,:),allocatable :: clumping_grid
  real(kind=dp) :: avg_dens !< average density
  character(len=512) :: clumping_fit_file
  public :: set_clumping, clumping_point
  ! LLS data
  real(kind=dp),parameter :: opdepth_LL = 2.0 !< typical optical depth of LLS
  real(kind=dp),parameter :: N_1 = opdepth_LL / sigma_HI_at_ion_freq !< typical column density of LLS
  real(kind=dp),public :: n_LLS
  real(kind=dp),public :: coldensh_LLS = 0.0_dp ! Column density of LLSs per cell
  real(kind=dp),public :: mfp_LLS_pMpc
  real,dimension(:,:,:),allocatable :: LLS_grid
  ! LLS parameters
  ! a) Model Prochaska et al. (2010)
  !real(kind=dp),parameter :: C_LLS = 1.9
  !real(kind=dp),parameter :: z_x = 3.7
  !real(kind=dp),parameter,public :: y_LLS = 5.1
  ! b) Model Songaila & Cowie (2010)
  real(kind=dp),parameter :: C_LLS = 2.84
  real(kind=dp),parameter :: z_x = 3.5
  real(kind=dp),parameter,public :: y_LLS = 2.04
  real(kind=dp),parameter :: beta=1.28
  ! c) Model McQuinn et al. (2011)
  !real(kind=dp),parameter :: C_LLS = 2.34
  !real(kind=dp),parameter :: z_x = 3.5
  !real(kind=dp),parameter,public :: y_LLS = 2.85
  !real(kind=dp),parameter :: beta=1.3
  
  public :: set_LLS, LLS_point  
  
#ifdef MPI
  integer,private :: mympierror
#endif

   type ionstates    
     real(kind=dp) :: h(0:1)          !< H  ionization fractions        
     real(kind=dp) :: he(0:2)         !< He ionization fractions           
     real(kind=dp) :: h_av(0:1)       !< average H  ionization fractions        
     real(kind=dp) :: he_av(0:2)      !< average He ionization fractions 
     real(kind=dp) :: h_old(0:1)      !< H  ionization fractions from last time step
     real(kind=dp) :: he_old(0:2)     !< He ionization fractions from last time step
  end type ionstates

contains

  ! ============================================================================

  subroutine mat_ini (restart, nz0, ierror)

    ! Initializes material properties on grid

    ! Authors: Garrelt Mellema, Ilian Iliev

    ! Date: 30-Jan-2008 (20-Aug-2006 (f77 21-May-2005 (derives from mat_ini_cosmo2.f)))

    ! Version: 
    ! - Three-dimensional. 
    ! - Cosmological density field (read from file)
    ! - Initially completely neutral
    ! - Source points set elsewhere (multiple sources)

    ! History:
    ! - integrated with new cosmological approach developed for
    !    test 4.
    ! - adapted for multiple cosmological sources
    ! - f90 version with MPI

    integer,intent(out) :: restart ! will be /= 0 if a restart is intended
    integer,intent(out) :: nz0    ! nz0 is the number of the starting slice
    integer,intent(out) :: ierror ! will be /=0 if an error occurred

    integer :: i,j,k,n ! loop counters
    character(len=1) :: answer, isothermal_answer

    ierror=0
    ! Check consistency with nbody module
    if (nbody_type /= "cubep3m") then
       write(logf,*) "Error: wrong material module was compiled."
       write(logf,*) "       Expected cubep3m, but got", &
            trim(adjustl(nbody_type))
       ierror=1
    else
       ! restart
       restart=0 ! no restart by default

       if (rank == 0) then
          ! Ask for temperature, restart. Read in values
          if (.not.file_input) &
               write(*,"(A,$)") "Enter initial temperature (K): "
          read(stdinput,*) temper_val
          
          if (.not.file_input) &
               write(*,"(A,$)") "isothermal? y/n: "
          read(stdinput,*) isothermal_answer
          if (isothermal_answer.eq.'n')then
             isothermal=.false.
          elseif (isothermal_answer.eq.'y') then
             isothermal=.true.
          else 
             write(*,"(A,$)") "mistake you should write y or n" 
             write(*,"(A,$)") ' '
             call exit(0)                
          endif
                    
          if (.not.file_input) write(*,"(A,$)") "Restart (y/n)? : "
          read(stdinput,*) answer
          if (answer == "y" .or. answer == "Y") then
             restart=1
             ! Apparently confusing if this question is only asked
             ! conditionally
             !if (.not.file_input) &
             !     write(*,"(A,$)") "Restart at midpoint (y/n)? : "
             !read(stdinput,*) answer
             !if (answer == "y" .or. answer == "Y") restart=2
          endif
          if (.not.file_input) &
               write(*,"(A,$)") "Restart at midpoint (y/n)? : "
          read(stdinput,*) answer
          if (answer == "y" .or. answer == "Y") restart=2
          if (.not.file_input) write(*,"(A,$)") "Number of starting slice: "
          read(stdinput,*) nz0
          if (.not.file_input) &
               write(*,"(A,$)") "Clumping fit table file: "
          read(stdinput,*) clumping_fit_file
          !print*,"Number of starting slice: ",nz0
       endif
#ifdef MPI       
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(isothermal,1,MPI_LOGICAL, 0,MPI_COMM_NEW,mympierror)       
       call MPI_BCAST(temper_val,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(restart,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(nz0,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
       
       !*** initialize the collisional ion and recomb rates for inital temp
       call ini_rec_colion_factors(temper_val)        ! Scalar version for constant temperature

       
       ! Allocate density array
       allocate(ndens(mesh(1),mesh(2),mesh(3)))
       ! Assign dummy density to the grid
       ! This should be overwritten later (in dens_ini)
       ndens(:,:,:)=1.0
       
       ! Allocate temperature array and initialize if the run is not
       ! isothermal
       if (.not.isothermal) then
          allocate(temperature_grid(mesh(1),mesh(2),mesh(3),0:2))
          temperature_grid(:,:,:,:)=temper_val
       endif
       ! Report on temperature situation
       if (rank == 0) then
          if (isothermal) then
             write(logf,"(A)") "Thermal conditions: isothermal"
          else
             write(logf,"(A)") &
                  "Thermal conditions: applying heating and cooling"
          endif
       endif 
 
       
       ! Allocate ionization fraction arrays
       allocate(xh(mesh(1),mesh(2),mesh(3),0:1))
       allocate(xhe(mesh(1),mesh(2),mesh(3),0:2))
       ! Assign ionization fractions. For z = 40 - 20 RECFAST gives 
       ! an average ionization fraction of about 2e-4 for H and
       ! 1e-15 for He. We use this here.
       ! In case of a restart this will be overwritten in xfrac_ini
       xh(:,:,:,1)=2e-4
       xh(:,:,:,0)=1.0_dp-xh(:,:,:,1)
       xhe(:,:,:,1)=1e-15_dp
       xhe(:,:,:,0)=1.0_dp-xhe(:,:,:,1)
       xhe(:,:,:,2)=0.0
       ! Initialize LLS parametets
       call LLS_init ()       
    endif

  end subroutine mat_ini

  ! ===========================================================================

  subroutine dens_ini (zred_now,nz)

    ! Initializes density on the grid (at redshift zred_now)

    ! Authors: Garrelt Mellema, Ilian Iliev

    ! Date: 30-Jan-2008 (20-Aug-2006 (19-May-2005 (8-mar-2005, 23-Nov-2004, 02-Jun-2004)))

    ! Version: 
    ! - Three-dimensional. 
    ! - Cosmological density field (read from file, depends on redshift)
    ! - Source points set elsewhere (multiple sources)
    ! - PMFAST input
    ! - MPI

    real(kind=dp),intent(in) :: zred_now
    integer,intent(in) :: nz ! number in the list of redshifts (for
                             ! compatibility reasons)
    
    integer :: i,j,k,n,nfile ! loop counters
    character(len=512):: dens_file
    character(len=6) :: zred_str
    character(len=1) :: nfile_str
    real(kind=dp) :: convert ! conversion factor
    real(kind=dp) :: summed_density
    real(kind=si),dimension(:),allocatable :: partial_sum_ndens
    integer :: m1,m2,m3

    ! density in file is in 4B reals, read in via this array
    real(kind=si),dimension(:,:,:),allocatable :: ndens_real

    if (rank == 0) then
       ! construct filename
       write(zred_str,"(f6.3)") zred_now
       dens_file=trim(adjustl(dir_dens))// &
            trim(adjustl(zred_str))// &
            !"rho_"//trim(adjustl(id_str))//".dat"
            "n_all.dat"
       write(unit=logf,fmt="(4A)") "Reading ",id_str, &
            " density input from ",trim(dens_file)
       
       ! Open density file: note that it is in `binary" form
       open(unit=20,file=dens_file,form=densityformat, &
            access=densityaccess,status="old")
       
       ! Read in data
       if (densityheader) then
          read(20) m1,m2,m3
          if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
             write(logf,*) "Warning: file with densities unusable"
             write(logf,*) "mesh found in file: ",m1,m2,m3
             stop
          endif
       endif
       ! Read in data and store it in ndens
       ! Allocate array needed to read in data
       allocate(ndens_real(mesh(1),mesh(2),mesh(3)))
       read(20) ndens_real
       ndens(:,:,:)=ndens_real(:,:,:)
          
       ! close file
       close(20)
       
       ! Deallocate array needed for reading in the data.
       deallocate(ndens_real)
    endif
#ifdef MPI       
    ! Distribute the density to the other nodes
    call MPI_BCAST(ndens,mesh(1)*mesh(2)*mesh(3),MPI_DOUBLE_PRECISION,0,&
         MPI_COMM_NEW,mympierror)
#endif
       
    ! The original values in terms of the mean density
    ! Below is the conversion factor
    ! vol is redshift dependent, that is why we need to recalculate this
    ! M_particle and M_grid should be in g
    select case(density_unit)
    case ("grid")
       convert=density_convert_grid 
    case ("particle")
       convert=density_convert_particle
    case ("M0Mpc3")
       convert=M_solar/Mpc**3*h**2*Omega_B/Omega0/(mu*m_p)
    end select
    
    ! Check separately for cosmological.
    if (cosmological) then
       convert=convert*(1.0+zred_now)**3
    endif
      
    ! Assign density to the grid
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             ndens(i,j,k)=ndens(i,j,k)*convert
             ! take care of empty cells (0.1 particles)
             if (ndens(i,j,k)<= 0.0) ndens(i,j,k)=0.1*convert
          enddo
       enddo
    enddo
    
    ! Find summed and average density
    ! GM/110208: Make this a partial sum (since we are dealing
    ! with bigger meshes now.
    allocate(partial_sum_ndens(mesh(3)))
    do k=1,mesh(3)
       partial_sum_ndens(k)=sum(ndens(:,:,k))
    enddo
    avg_dens=sum(partial_sum_ndens)/ &
         (real(mesh(1),dp)*real(mesh(2),dp)*real(mesh(3),dp))
    deallocate(partial_sum_ndens)
    
    
    
    ! Report density field properties
    if (rank == 0) then
       write(logf,*) "Raw density diagnostics (cm^-3)"
       write(logf,*) "minimum density: ",minval(ndens)
       write(logf,*) "maximum density: ",maxval(ndens)
       write(logf,"(A,es10.3,A)") "Average density = ",avg_dens," cm^-3"
       write(logf,"(A,es10.3,A)") "Theoretical value = ", &
            rho_crit_0*Omega_B/(mu*m_p)*(1.0+zred_now)**3, &
            " cm^-3" 
       write(logf,"(A,es10.3,A)") "(at z=0 : ", &
            rho_crit_0*Omega_B/(mu*m_p), &
            " cm^-3)"
    endif

  end subroutine dens_ini

  ! ===========================================================================

  subroutine xfrac_ini (zred_now)

    ! Initializes ionization fractions on the grid (at redshift zred_now).
    ! They are read from an "Ifront3" file which should have been created
    ! by an earlier run. This file is read in restart mode.

    ! Author: Garrelt Mellema

    ! Date: 19-May-2005

    real(kind=dp),intent(in) :: zred_now
    
    character(len=512) :: xfrac_file
    character(len=512) :: xfrac_file_He1
    character(len=512) :: xfrac_file_He2
    character(len=6) :: zred_str
    integer :: m1,m2,m3
    ! Array needed to read in 4B reals
    real(kind=dp),dimension(:,:,:),allocatable :: xh1_real
    real(kind=dp),dimension(:,:,:),allocatable :: xhe1_real
    real(kind=dp),dimension(:,:,:),allocatable :: xhe2_real
    !real(kind=si),dimension(:,:,:),allocatable :: xh1_real

    if (rank == 0) then
       allocate(xh1_real(mesh(1),mesh(2),mesh(3)))
       allocate(xhe1_real(mesh(1),mesh(2),mesh(3)))
       allocate(xhe2_real(mesh(1),mesh(2),mesh(3)))
       write(zred_str,"(f6.3)") zred_now
!       xfrac_file= "./xfrac3d_"//trim(adjustl(zred_str))//".bin"
       xfrac_file= trim(adjustl(results_dir))// &
            !"Ifront3_"//trim(adjustl(zred_str))//".bin"
            "xfrac3d_"//trim(adjustl(zred_str))//".bin"
         xfrac_file_He1= trim(adjustl(results_dir))// &
            "xfrac3dHe1_"//trim(adjustl(zred_str))//".bin"   
         xfrac_file_He2= trim(adjustl(results_dir))// &
            "xfrac3dHe2_"//trim(adjustl(zred_str))//".bin"

       write(unit=logf,fmt="(2A)") "Reading ionization fractions from ", &
            trim(xfrac_file), "and", trim(xfrac_file_He1), "and",  trim(xfrac_file_He2)
       ! Open ionization fractions file
       open(unit=20,file=xfrac_file,form="unformatted",status="old")
       
       ! Read in data
       read(20) m1,m2,m3
       if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
          write(logf,*) "Warning: file with ionization fractions unusable"
          write(logf,*) "mesh found in file: ",m1,m2,m3
       else
          read(20) xh1_real
          ! To avoid xh(0)=0.0, we add a nominal ionization fraction of 10^-12
          xh(:,:,:,1)=xh1_real(:,:,:)
          !xh(:,:,:,1)=real(xh1_real(:,:,:),dp)
          xh(:,:,:,0)=1.0_dp-xh(:,:,:,1)
       endif
       ! close file
       close(20)
       deallocate(xh1_real)

!**** He1
       ! Open ionization fractions file
       open(unit=20,file=xfrac_file_He1,form="unformatted",status="old")       
       ! Read in data
       read(20) m1,m2,m3
       if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
          write(logf,*) "Warning: file with ionization fractions unusable"
          write(logf,*) "mesh found in file: ",m1,m2,m3
       else
          read(20) xhe1_real
          ! To avoid xh(0)=0.0, we add a nominal ionization fraction of 10^-12
          xhe(:,:,:,1)=xhe1_real(:,:,:)
       endif
       ! close file
       close(20)
       deallocate(xhe1_real)
!**** He2
       ! Open ionization fractions file
       open(unit=20,file=xfrac_file_He2,form="unformatted",status="old")       
       ! Read in data
       read(20) m1,m2,m3
       if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
          write(logf,*) "Warning: file with ionization fractions unusable"
          write(logf,*) "mesh found in file: ",m1,m2,m3
       else
          read(20) xhe2_real
          ! To avoid xh(0)=0.0, we add a nominal ionization fraction of 10^-12
          xhe(:,:,:,2)=xhe2_real(:,:,:)
          xhe(:,:,:,0)=1.0_dp-xhe(:,:,:,1)-xhe(:,:,:,2)
       endif
       ! close file
       close(20)
       deallocate(xhe2_real)
    endif

#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(xh,mesh(1)*mesh(2)*mesh(3)*2,MPI_DOUBLE_PRECISION,0,&
         MPI_COMM_NEW,mympierror)
    call MPI_BCAST(xhe,mesh(1)*mesh(2)*mesh(3)*3,MPI_DOUBLE_PRECISION,0,&
         MPI_COMM_NEW,mympierror)
#endif
    
  end subroutine xfrac_ini

  ! ===========================================================================
  
  subroutine protect_ionization_fractions(xfrac,lowfrac,highfrac,ifraction, &
       name)

    integer,intent(in) :: lowfrac
    integer,intent(in) :: highfrac
    integer,intent(in) :: ifraction
    real(kind=dp),dimension(mesh(1),mesh(2),mesh(3),lowfrac:highfrac), &
         intent(inout) :: xfrac
    character(len=*) :: name

    integer :: i,j,k

    ! Check the fractions for negative values
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             if (xfrac(i,j,k,ifraction) < 0.0d0) then
#ifdef MPILOG
                write(logf,*) name,' < 0 at ',i,j,k,xfrac(i,j,k,ifraction)
#endif
                xfrac(i,j,k,ifraction)=0.0_dp
             endif
             if (xfrac(i,j,k,ifraction) > 1.0d0) then
#ifdef MPILOG
                write(logf,*) name,' > 1 at ',i,j,k,xfrac(i,j,k,ifraction)
#endif
                xfrac(i,j,k,ifraction)=1.0_dp
             endif
          enddo
       enddo
    enddo
    
  end subroutine protect_ionization_fractions

  ! ===========================================================================

  subroutine temper_ini (zred_now)

    ! Initializes temperature on the grid (at redshift zred_now).
    ! They are read from an "Temper3D" file which should have been created
    ! by an earlier run. This file is read when in restart mode.

    ! Author: Garrelt Mellema

    ! Date: 02-June-2011

    real(kind=dp),intent(in) :: zred_now
    
    character(len=512) :: temper_file
    character(len=6) :: zred_str
    integer :: m1,m2,m3

    if (isothermal) then
       if (rank == 0) write(logf,"(A)") &
            "Incorrect call to temper_ini in isothermal case"
    else
       if (rank == 0) then
          write(zred_str,"(f6.3)") zred_now
          temper_file= trim(adjustl(results_dir))// &
               "Temper3D_"//trim(adjustl(zred_str))//".bin"
          
          write(unit=logf,fmt="(2A)") "Reading temperature from ", &
               trim(temper_file)
          ! Open temperature file
          open(unit=20,file=temper_file,form="unformatted",status="old")
          
          ! Read in data
          read(20) m1,m2,m3
          if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
             write(logf,*) "WARNING: file with temperatures unusable, as"
             write(logf,*) "mesh found in file: ",m1,m2,m3
          else
             read(20) temperature_grid(:,:,:,0)
          endif
          
          ! Fill the other parts of the temperature grid array
          ! See evolve for their use
          temperature_grid(:,:,:,1)=temperature_grid(:,:,:,0)
          temperature_grid(:,:,:,2)=temperature_grid(:,:,:,0)

          ! close file
          close(20)
       endif
       
#ifdef MPI       
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(temperature_grid,mesh(1)*mesh(2)*mesh(3)*3,MPI_REAL,0,&
            MPI_COMM_NEW,mympierror)
#endif
    endif

  end subroutine temper_ini

  ! ===========================================================================

  subroutine get_temperature_point (i,j,k,temper_inter,av_temper,temper)

    ! Puts value of temperature (from grid or initial condition value)
    ! in the module variable temper

    integer,intent(in) :: i,j,k
    real(kind=dp),intent(out) :: temper,av_temper,temper_inter

    if (isothermal) then
       temper = dble(temper_val)
       av_temper=dble(temper_val)
       temper_inter=dble(temper_val)
    else
       temper_inter = dble(temperature_grid(i,j,k,0))
       av_temper = dble(temperature_grid(i,j,k,1))   
       temper = dble(temperature_grid(i,j,k,2))            
    endif


  end subroutine get_temperature_point

  ! ===========================================================================

  subroutine set_temperature_point (i,j,k,temper_inter,av_temper)
    
    ! Puts value of module variable temper back in temperature grid
    ! (if running not isothermal)

    integer,intent(in) :: i,j,k
    real(kind=dp),intent(in) :: temper_inter,av_temper
    
    if (.not.isothermal) temperature_grid(i,j,k,0)=real(temper_inter)
    if (.not.isothermal) temperature_grid(i,j,k,1)=real(av_temper)    
    
  end subroutine set_temperature_point


  subroutine set_final_temperature_point ()
    
    ! Puts value of module variable temper back in temperature grid
    ! (if running not isothermal)

    !integer,intent(in) :: i,j,k
    !real(kind=dp),intent(in) :: temper
    
    if (.not.isothermal) temperature_grid(:,:,:,2)=temperature_grid(:,:,:,0)
    
  end subroutine set_final_temperature_point  

  ! ===========================================================================  
  subroutine set_clumping(z)

    !! 1: constant clumping (with clumping_factor)\n
    !! 2: 3.5Mpc PM, WMAP1 clumping\n
    !! 3: 3.5Mpc PM, WMAP3 clumping\n
    !! 4: 1 Mpc P3M\n
    !! 5: position dependent clumping

    real(kind=dp),intent(in) :: z

    select case (type_of_clumping)
    case(1)
       clumping = clumping_factor
    case(2) 
       clumping = 27.466*exp(-0.114*z+0.001328*z*z)
    case(3)
       clumping = 26.2917*exp(-0.1822*z+0.003505*z*z)
    case(4)
       clumping = 17.57*exp(-0.101*z+0.0011*z*z)
    case(5)
       call clumping_init (z)
    end select

    if (rank == 0) write(logf,*) "Setting (mean) global clumping factor to ", &
         clumping,"(type ", type_of_clumping,")"
    
  end subroutine set_clumping

  ! ===========================================================================

  subroutine clumping_point (i,j,k)

    integer,intent(in) :: i,j,k

    if (type_of_clumping /= 5) then
       write(logf,*) &
         "Error: calling position dependent, but array is not initialized."
    else
       clumping = clumping_grid(i,j,k)
    endif

  end subroutine clumping_point

  ! ===========================================================================

  subroutine clumping_init (zred_now)

    ! Initializes position dependent clumping (at redshift zred_now)

    ! Author: Ilian Iliev, modified from code by Garrelt Mellema

    ! Date: 03-Mar-2008 (5-Mar-2007( 20-Aug-2006, 19-May-2005 (8-mar-2005, 23-Nov-2004, 02-Jun-2004))

    ! Version: 
    ! - Three-dimensional. 
    ! - gas clumping field (read from file, depends on redshift)
    ! - Source points set elsewhere (multiple sources)
    ! - PMFAST input
    ! - MPI

    real(kind=dp),intent(in) :: zred_now
    
    integer :: i,j,k,n,nfile ! loop counters
    integer :: m1,m2,m3 ! size of mesh in clumping file (header)
    character(len=512):: clump_file
    character(len=6) :: zred_str
    character(len=1) :: nfile_str
    real(kind=dp) :: z_read
    real(kind=dp) :: a0
    real(kind=dp) :: a1
    real(kind=dp) :: a2
    real(kind=dp) :: error
    real(kind=dp) :: avg_dens !< average density
    integer :: io_status

    ! clumping in file is in 4B reals, read in via this array
    !real(kind=si),dimension(:,:,:),allocatable :: clumping_real

    if (.not.(allocated(clumping_grid))) &
         allocate(clumping_grid(mesh(1),mesh(2),mesh(3)))

    if (rank == 0) then
       ! construct filename
       write(zred_str,"(f6.3)") zred_now
       write(30,*) "Reading ",id_str," input"
       clump_file=trim(adjustl(dir_clump))// &
            trim(adjustl(zred_str))// &
            "n_all.dat"
            !"c_all.dat"
            !"clump_"//trim(adjustl(id_str))//".dat"
       ! Clumping file should be density with halos included.
       ! This most likely has an identical file name to density
       ! with halos excluded, but resides in a different directory.
       ! GM, 2009-12-09

       write(unit=logf,fmt="(4A)") "Reading ",id_str, &
            " clumping input from ",trim(clump_file)
       ! Open clumping file: note that the format is determined
       ! by the values of clumpingformat and clumping access,
       ! both set in the nbody module.
       open(unit=20,file=clump_file,form=clumpingformat, &
            access=clumpingaccess,status="old")
       
       ! Read in data
       ! Read in header if there is one
       if (clumpingheader) then
          read(20) m1,m2,m3
          if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
             write(logf,*) "Warning: file with clumping factors unusable"
             write(logf,*) "mesh found in file: ",m1,m2,m3
             stop
          endif
       endif
       ! Read in data and store it in clumping_grid
       read(20) clumping_grid
       write(logf,*) 'Clumping data read'
       !clumping_grid(:,:,:)=clumping_real(:,:,:)
       
       ! close file
       close(20)
       
       ! Normalize to average density
       avg_dens=sum(clumping_grid)/ &
            (real(mesh(1),dp)*real(mesh(2),dp)*real(mesh(3),dp))
       clumping_grid=clumping_grid/avg_dens
       clumping_grid=log10(clumping_grid*clumping_grid)
       ! Report on data: min, max, total
       ! Disabled (GM 2009-12-09): if we read in a density field
       ! there is no sense in reporting these.
       ! assign mean to clumping for reporting in set_clumping
       !clumping=sum(clumping_grid)/(mesh(1)*mesh(2)*mesh(3))
       !write(logf,*) "Statistics BEFORE applying clumping fit"
       !write(logf,*) "minimum: ",minval(clumping_grid)
       !write(logf,*) "maximum: ",maxval(clumping_grid)
       !write(logf,*) "average clumping: ",clumping

       ! Read clumping fit (for clumping based on smaller scales
       ! This file contains the parameters for the fit
       ! y=a0+a1*x+a2*x²
       ! where x is the alog10(<n²>_int) and y is alog10(<n²>_Jeans).
       open(unit=21,file=trim(adjustl(dir_clump))//clumping_fit_file, &
          status="old",form="formatted")
       z_read=0.0
       io_status=0
       ! Read in lines until the correct redshift is found.
       do while(z_read-zred_now > 1e-2 .and. io_status == 0)
          read(unit=21,iostat=io_status) z_read,a0,a1,a2,error
       enddo
       close(unit=21)
       if (io_status /= 0) then
          ! Something went wrong, assume no clumping and warn
          a0=0.0
          a1=1.0
          a2=0.0
          write(logf,*) "WARNING: no clumping data found in ", & 
               trim(adjustl(clumping_fit_file))," for redshift ",zred_now
       endif
       clumping_grid(:,:,:)=10.0**(a0+ &
            (a1-1.0)*clumping_grid(:,:,:)+ &
            a2*clumping_grid(:,:,:)*clumping_grid(:,:,:))
    endif

#ifdef MPI       
    ! Distribute the density to the other nodes
    call MPI_BCAST(clumping_grid,mesh(1)*mesh(2)*mesh(3), &
         MPI_REAL,0,MPI_COMM_NEW,mympierror)
#endif
       
    ! Report on data: min, max, total
    ! assign mean to clumping for reporting in set_clumping
    clumping=sum(clumping_grid)/(mesh(1)*mesh(2)*mesh(3))
    if (rank == 0) then
       write(logf,*) "Statistics AFTER applying clumping fit"
       write(logf,*) "clumping fit for redshift ",z_read
       write(logf,*) "minimum: ",minval(clumping_grid)
       write(logf,*) "maximum: ",maxval(clumping_grid)
       write(logf,*) "average clumping: ",clumping
    endif

  end subroutine clumping_init
  ! ===========================================================================                        
  subroutine LLS_init ()
    
#ifdef IFORT
    !For gamma function
    use ISO_C_BINDING
#endif
    
#ifdef IFORT
    !For gamma function
    interface
       real(c_double) function tgamma (y) bind(c)
         use iso_c_binding
         real(c_double), value :: y
       end function tgamma
    end interface
#endif

    if (type_of_LLS == 1) then
       ! 1/distance between LLSs expressed in grid cells (z=0)
       n_LLS = C_LLS * (1.0/(1.0 + z_x)) ** y_LLS * dr(1) * H0*sqrt(Omega0) / c

       !n_LLS=n_LLS * ((1.0 + zred)/(1.0+zred_prev))** (y_LLS+1.5)
       !mfp=c/((1.0+z) * Hz * C_LLS * ((1.0 + z)/(1.0 + z_x)) ** y_LLS )
       ! Add the beta correction as explained in Songaila & Cowie (2010).
       ! This corrects for the fact that not all LLS have the same
       ! column density. beta is the slope of the distribution function
       ! of LLS over column densities.
       ! This expression needs the gamma function. For the intel compiler
       ! this is tgamma(). For other compilers (Fortran 2008 standard) this
       ! is gamma().
#ifdef IFORT    
       n_LLS=n_LLS*tgamma(2.0-beta)/(opdepth_LL**(1.0-beta))
#else
       n_LLS=n_LLS*gamma(2.0-beta)/(opdepth_LL**(1.0-beta))
#endif
    else
       ! Set distance between LLS (and mean free path) to infinity
       ! If type_of_LLS = 2 this will be overwritten by cell specific
       ! values.
       n_LLS=0.0d0
    endif

  end subroutine LLS_init

  ! ===========================================================================

  subroutine set_LLS (z)

    !! Two cases:\n
    !! 1: constant LLS optical depth per cell\n
    !! 2: position dependent LLS optical depth

    real(kind=dp),intent(in) :: z

    select case (type_of_LLS)
    case(1)
       ! Column density per cell due to LLSs
       coldensh_LLS = N_1 * n_LLS
       mfp_LLS_pMpc=dr(1)/n_LLS
    case(2) 
       call read_lls_grid (z)
    end select

    if (rank == 0) then
       write(logf,*) "Average optical depth per cell due to LLSs: ", &
            coldensh_LLS*sigma_HI_at_ion_freq,"(type ", type_of_LLS,")"
       write(logf,*) "Mean free path (pMpc): ", mfp_LLS_pMpc
    endif
    
  end subroutine set_LLS

  ! ===========================================================================

  subroutine LLS_point (i,j,k)

    integer,intent(in) :: i,j,k

    if (type_of_LLS /= 2) then
       write(logf,*) &
         "Error: calling position dependent LLS, but array is not initialized."
    else
       coldensh_LLS = LLS_grid(i,j,k)
    endif

  end subroutine LLS_point

  ! ===========================================================================

  subroutine read_LLS_grid (zred_now)

    ! Initializes position dependent LLS optical depth (at redshift zred_now)

    ! Author: Garrelt Mellema

    ! Date: 16-Mar-2011

    ! Version: 

    real(kind=dp),intent(in) :: zred_now
    
    integer :: m1,m2,m3 ! size of mesh in cross section file (header)
    character(len=512):: LLS_file
    character(len=6) :: zred_str

    ! clumping in file is in 4B reals, read in via this array
    !real(kind=si),dimension(:,:,:),allocatable :: clumping_real

    if (.not.(allocated(LLS_grid))) &
         allocate(LLS_grid(mesh(1),mesh(2),mesh(3)))

    if (rank == 0) then
       ! construct filename
       write(zred_str,"(f6.3)") zred_now
       LLS_file=trim(adjustl(dir_LLS))// &
            trim(adjustl(zred_str))// &
            "cross_section.bin"

       write(unit=logf,fmt="(4A)") "Reading ",id_str, &
            " clumping input from ",trim(LLS_file)

       ! Open clumping file: note that the format is determined
       ! by the values of clumpingformat and clumping access,
       ! both set in the nbody module.
       open(unit=22,file=LLS_file,form=LLSformat, &
            access=LLSaccess,status="old")
       
       ! Read in data
       ! Read in header if there is one
       if (LLSheader) then
          read(22) m1,m2,m3
          if (m1 /= mesh(1).or.m2 /= mesh(2).or.m3 /= mesh(3)) then
             write(logf,*) "Warning: file with LLS cross sections unusable"
             write(logf,*) "mesh found in file: ",m1,m2,m3
             stop
          endif
       endif
       ! Read in data and store it in clumping_grid
       read(22) LLS_grid
       write(logf,*) 'LLS data read'
       
       ! close file
       close(unit=22)
       
       ! Calculate mean free path
       ! Make sure sim_volume is in proper length units
       mfp_LLS_pMpc=sim_volume/(sum(LLS_grid)*Mpc*(1.0+zred_now)**3)

       ! Convert to column density
       LLS_grid=LLS_grid/vol ! 1/(mean free path) = n_LLS
       LLS_grid=N_1 * LLS_grid
    endif

#ifdef MPI       
    ! Distribute the density to the other nodes
    call MPI_BCAST(LLS_grid,mesh(1)*mesh(2)*mesh(3), &
         MPI_REAL,0,MPI_COMM_NEW,mympierror)
#endif
       
    ! Report on data: min, max, total
    ! assign mean to clumping for reporting in set_clumping
    coldensh_LLS=sum(LLS_grid)/(mesh(1)*mesh(2)*mesh(3))
    if (rank == 0) then
       write(logf,*) "Statistics on LLS column density"
       write(logf,*) "minimum: ",minval(LLS_grid)
       write(logf,*) "maximum: ",maxval(LLS_grid)
       write(logf,*) "average clumping: ",coldensh_LLS
       write(logf,*) "mean free path: ",mfp_LLS_pMpc
    endif

  end subroutine read_LLS_grid
end module material
