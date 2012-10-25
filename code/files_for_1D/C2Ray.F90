!>
!! \brief Main program for C2Ray-1D
!!
!! C2Ray-1D does a one-dimensional photo-ionization calculation 
!! one of a series of test problems.\n
!! The main programme calls a number of initialization routines
!! and then enters the main integration loop, which ends when one
!! of the stopping conditions is met. At specified times it calls
!! the output module routines to produce output.
!! After the integration loop ends, a number of closing down routines
!! are called and the programme stops.
!!
!! \b Author: Garrelt Mellema \n
!!
!! \b Date: 23-Sep-2006
!<
Program C2Ray

  ! Author: Garrelt Mellema

  ! Date: 23-Sep-2006

  ! Goal:
  ! One dimensional photo-ionization calculation for a series of
  ! test problems.
  
  ! Version notes:
  ! - Does not include hydrodynamics
  ! - Assumes time step

  ! Needs following modules
  use precision, only: dp
  use clocks, only: setup_clocks, update_clocks, report_clocks
  use file_admin, only: stdinput, logf, file_input, flag_for_file_input
  use astroconstants, only: YEAR
  use my_mpi, only: mpi_setup, mpi_end, rank
  use output_module, only: setup_output,output,close_down
  use grid, only: grid_ini
  use radiation, only: rad_ini,sourcetype
  use cosmology, only: cosmology_init, redshift_evol, &
       time2zred, zred2time, zred,  cosmo_evol
  use c2ray_parameters, only: cosmological
  !use cosmological_evolution, only: cosmo_evol
  use material, only: mat_ini, testnum,zred00,isothermal
  use times, only: time_ini, number_timesteps
  use sourceprops, only: source_properties_ini
  use evolve, only: evolve1D
  use times, only: end_time,dt,output_time, timess, outputsteps
  use radiative_cooling, only: setup_cool
#ifdef XLF
  ! Modules for the xlf (IBM) compiler
  USE XLFUTILITY, only: iargc, getarg, flush => flush_
#endif

  implicit none

  ! Integer variables
  integer :: nstep !< time step counter
  integer :: restart !< restart if not zero (not used in 1D code)

  ! Time variables
  real(kind=dp) :: sim_time !< actual time (s)
  real(kind=dp) :: next_output_time !< time of next output (s)
  real(kind=dp) :: actual_dt !< actual time step (s)

  !> Input file
  character(len=512) :: inputfile

  ! Initialize clocks (cpu and wall)
  call setup_clocks

  ! Set up MPI structure (compatibility mode) & open log file
  call mpi_setup()

  ! Set up input stream (either standard input or from file given
  ! by first argument)
  if (rank == 0) then
     write(logf,*) "screen input or file input?"
     flush(logf)
     if (COMMAND_ARGUMENT_COUNT () > 0) then
        call GET_COMMAND_ARGUMENT(1,inputfile)
        write(logf,*) "reading input from ",trim(adjustl(inputfile))
        open(unit=stdinput,file=inputfile,status="old")
        call flag_for_file_input(.true.)
     else
        write(logf,*) "reading input from command line"
     endif
     flush(logf)
  endif

  ! Initialize output
  call setup_output ()

  ! Initialize grid
  call grid_ini ()

  ! Initialize the material properties (calls also ini_rec_colion_factors)
  call mat_ini (restart)

  if (.not.isothermal) call setup_cool () ! moved here from radiation

  ! Initialize photo-ionization calculation
  call rad_ini( )

  ! Initialize source property
  call source_properties_ini (sourcetype)

  ! Initialize time step parameters
  call time_ini ()

  ! Set time to zero
  sim_time=0.0
  next_output_time=0.0

  ! Initialize and update cosmology (transform from comoving to proper values)
  if (cosmological) then
     write(logf,*) "Cosmology initialization 1 (redshift, time, redshift0): ", &
          zred, sim_time,zred00
     call cosmology_init(zred00,sim_time)
     call redshift_evol(sim_time)
     write(logf,*) "Cosmology initialization 2 (redshift, time): ", &
           zred, sim_time
     call cosmo_evol( )
     write(logf,*) "Cosmology initialization 3 (redshift, time): ", zred
  endif

  ! Record initial conditions
  call output(0,0.0d0,timess(1),end_time)

  ! Loop until end time is reached
  do nstep=1,number_timesteps!-1

 	sim_time=timess(nstep) 
  	dt=sim_time-timess(nstep-1)
  	actual_dt=dt
     ! Write output
     !	if (mod(nstep,outputsteps).lt.1.0) then
     !   	call output(nstep,sim_time,dt,end_time)
     !   	write(50,*) sim_time/YEAR
     !   endif

!     ! Write output
!     if (abs(sim_time-next_output_time).le.1e-6*time) then
!        call output(sim_time,dt,end_time)
!        next_output_time=next_output_time+output_time
!     endif

     ! Make sure you produce output at the correct time
     ! dt=YEAR*10.0**(min(5.0,(-2.0+real(nstep)/1e5*10.0)))
!     actual_dt=min(abs(next_output_time-sim_time),dt)

!     nstep=nstep+1
     ! Report time and time step
     write(logf,'(A,2(1pe10.3,1x),A)') 'Time, dt:', &
          sim_time/YEAR,actual_dt/YEAR,' (years)'
!     write(48,*) '0', nstep 
     ! For cosmological simulations evolve proper quantities
     if (cosmological) then
        call redshift_evol(sim_time+0.5*actual_dt)
        call cosmo_evol()
     endif

     ! Take one time step
     call evolve1D(actual_dt)

     	if (mod(nstep,outputsteps).lt.1.0) then
        	call output(nstep,sim_time,dt,end_time)
        	write(50,*) sim_time/YEAR
        endif

     sim_time=sim_time+actual_dt     ! Update time  
!     if (abs(sim_time-end_time).lt.1e-6*end_time) exit     

     ! Update clock counters (cpu + wall, to avoid overflowing the counter)
     call update_clocks ()

  enddo
   nstep=number_timesteps

  ! Scale to the current redshift
  if (cosmological) then
     call redshift_evol(sim_time)
     call cosmo_evol()
  endif
!if (abs(sim_time-end_time).gt.100) then
!  call output(nstep,sim_time,dt,end_time)! Write final output
!  write(50,*) sim_time/YEAR
!endif

  ! Clean up some stuff
  call close_down ()

  ! Report clocks (cpu and wall)
  call report_clocks ()

  ! End the run
  call mpi_end ()

end Program C2Ray
