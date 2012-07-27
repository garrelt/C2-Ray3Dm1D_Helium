!>
!! \brief Main program for C2Ray-1D
!!
!! \b Author: Garrelt Mellema \n
!! \b Date: 23-Sep-2006
!<
Program C2Ray

  ! Author: Garrelt Mellema

  ! Date: 23-Sep-2006

  ! Goal:
  ! One dimensional photo-ionization calculation for a series of
  ! test problems.
  
  ! Does not include hydrodynamics
  ! Assumes constant time step

  ! Needs following `modules'
  ! c2ray_parameters : all the tunable parameters for the code
  ! my_mpi : sets up the MPI parallel environment
  ! output_module : output routines
  ! grid : sets up the grid
  ! radiation : radiation tools
  ! pmfast : interface to PMFAST output
  ! cosmology : cosmological utilities
  ! material : material properties
  ! times : time and time step utilities
  ! sourceprops : source properties
  ! evolve : evolve grid in time

  use precision, only: dp
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
  use file_admin, only:stdinput
  use times, only: end_time,dt,output_time, timess, outputsteps
  use radiative_cooling, only: setup_cool
#ifdef XLF
  USE XLFUTILITY, only: iargc, getarg, flush => flush_
#endif

  implicit none

  ! CPU time variables
  real :: tstart 			!< Start time for CPU report
  real :: tend 				!< End time for CPU report
  					! Wall clock time variables
  integer :: cntr1 			!< Start time wall clock
  integer :: cntr2 			!< End time wall clock
  integer :: countspersec 		!< counts per second (for wall clock time)   
  integer :: nstep
  integer :: restart,nz,flag
  					! Time variables
  real(kind=dp) :: time 		!< actual time (s)
  real(kind=dp) :: next_output_time 	!< time of next output (s)
  real(kind=dp) :: actual_dt 		!< actual time step (s)
  real(kind=dp) :: test, test2
  character(len=512) :: inputfile	!> Input file
  call cpu_time(tstart)		! Initialize cpu timer
  call system_clock(cntr1)		! Initialize wall cock times
  call mpi_setup()			! Set up MPI structure
  if (iargc() > 0) then 		! Set up input stream (either standard input or from file given
     call getarg(1,inputfile)  		! by first argument)
     if (rank == 0) then
        write(*,*) 'reading input from ',trim(adjustl(inputfile))
        open(unit=stdinput,file=inputfile)
     endif
  endif


  call setup_output()     		! Initialize output
  call grid_ini()         		! Initialize grid
  call mat_ini (restart)  		! Initialize the material properties (calls also ini_rec_colion_factors)
  if (.not.isothermal) call setup_cool () ! moved here from radiation
  call rad_ini( )         		! Initialize photo-ionization calculation
  call source_properties_ini (sourcetype)         ! Initialize source property
  call time_ini ()        		! Initialize time step parameters
  time=0.0                		! Set time to zero
  next_output_time=0.0 !*
  if (cosmological) then  		! Initialize cosmology     
     write(*,*) zred, time,zred00
     call cosmology_init(zred00,time)
     call redshift_evol(time)
     write(*,*) zred, time
     call cosmo_evol( )
     write(*,*) zred
  endif

  do  nstep=1,number_timesteps!-1

 	time=timess(nstep) 
  	dt=time-timess(nstep-1)
  	actual_dt=dt
     ! Write output
     !	if (mod(nstep,outputsteps).lt.1.0) then
     !   	call output(nstep,time,dt,end_time)
     !   	write(50,*) time/YEAR
     !   endif

!     ! Write output
!     if (abs(time-next_output_time).le.1e-6*time) then
!        call output(time,dt,end_time)
!        next_output_time=next_output_time+output_time
!     endif

     ! Make sure you produce output at the correct time
     ! dt=YEAR*10.0**(min(5.0,(-2.0+real(nstep)/1e5*10.0)))
!     actual_dt=min(abs(next_output_time-time),dt)

!     nstep=nstep+1
     ! Report time and time step
     write(30,'(A,2(1pe10.3,1x),A)') 'Time, dt:', &
          time/YEAR,actual_dt/YEAR,' (years)'
!     write(48,*) '0', nstep 
     ! For cosmological simulations evolve proper quantities
     if (cosmological) then
        call redshift_evol(time+0.5*actual_dt)
        call cosmo_evol()
     endif
     ! Take one time step

     call evolve1D(actual_dt) ! Take one time step

     	if (mod(nstep,outputsteps).lt.1.0) then
        	call output(nstep,time,dt,end_time)
        	write(50,*) time/YEAR
        endif

     time=time+actual_dt     ! Update time  
!     if (abs(time-end_time).lt.1e-6*end_time) exit     
  enddo
   nstep=number_timesteps

  if (cosmological) then  ! Scale to the current redshift
     call redshift_evol(time)
     call cosmo_evol()
  endif
!if (abs(time-end_time).gt.100) then
!  call output(nstep,time,dt,end_time)! Write final output
!  write(50,*) time/YEAR
!endif
!  write(*,*)'fit,~1,table=',mycountertest,mycountertest2,mycountertest3
  call close_down ()! Clean up some stuff
  call cpu_time(tend)! Find out CPU time
  call system_clock(cntr2,countspersec)
  write(30,*) 'CPU time: ',tend-tstart,' s'
  write(30,*) 'Wall clock time: ',(cntr2-cntr1)/countspersec,' s'
  call mpi_end ()  ! End the run

end Program C2Ray
