module clocks

  use file_admin, only: logf, timefile, results_dir
  use my_mpi, only: rank
  use parameter, only: sim_time, next_output_time, end_time
  use astroconstants, only: YEAR

  implicit none

  ! start time for call to CPU routine
  real :: start_cpu_time 
  ! end time for call to CPU routine
  real :: end_cpu_time 
  ! the CPU hours
  integer :: cpu_hour = 0 
  ! the CPU minutes
  integer :: cpu_minute = 0 
  ! the CPU seconds
  real :: cpu_second = 0.0 
  
  ! Wall clock time variables
  ! (use 8 byte integers to avoid the clock reset)
  integer(kind=8) :: start_wallclock_count !< Start time for call to wall clock routine
  integer(kind=8) :: end_wallclock_count !< End time for call to wall clock routine
  integer(kind=8) :: absolute_wallclock_count !< Asbolute start time wall clock

  integer(kind=8) :: counts_per_sec !< counts per second (for wall clock time)

  integer :: clock_hours=0 !< accumulates the wall clock hours
  integer :: clock_minutes=0 !< accumulates the wall clock minutes
  real :: clock_seconds=0.0 !< accumulates the wall clock seconds
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set up all the clocks and initialize the time recording file
  subroutine setup_clocks
    
    ! name of the time file
    character(len=512) :: filename

    call setup_cpuclock()
    call setup_wallclock()

    ! Initialize the time recording file - timefile.log
    if (rank .eq. 0) then
       filename = trim(adjustl(trim(adjustl(results_dir))//"timefile.log"))
       open(unit=timefile,file=filename,status="unknown",action="write")
       write(unit=timefile,fmt="(A)") "Timings file for C2-Ray run"
    endif

  end subroutine setup_clocks
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set up cpu clock
  subroutine setup_cpuclock
    
    ! Initialize cpu timer
    call cpu_time(start_cpu_time)
    
  end subroutine setup_cpuclock
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set up wall clock
  subroutine setup_wallclock
    
    ! Initialize wall cock timer
    call system_clock(start_wallclock_count)
    
    ! record the wall clock time when the program begins
    absolute_wallclock_count = start_wallclock_count
    
  end subroutine setup_wallclock
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Update all the clocks
  subroutine update_clocks
    
    call update_cpuclock()
    call update_wallclock()
    
  end subroutine update_clocks
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Update CPU clock
  subroutine update_cpuclock
    
    ! Find out intermediate CPU time 
    call cpu_time(end_cpu_time)
    cpu_second = cpu_second+real(end_cpu_time-start_cpu_time)
    start_cpu_time = end_cpu_time
    cpu_minute = cpu_minute+int(cpu_second)/60
    cpu_second = mod(cpu_second,60.0)
    cpu_hour = cpu_hour+cpu_minute/60
    cpu_minute = mod(cpu_minute,60)

  end subroutine update_cpuclock
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Update wall clock
  subroutine update_wallclock
    
    ! Find out intermediate wall clock time    
    call system_clock(end_wallclock_count,counts_per_sec)
    clock_seconds = clock_seconds+real(end_wallclock_count-start_wallclock_count)/real(counts_per_sec)
    start_wallclock_count = end_wallclock_count
    clock_minutes = clock_minutes+int(clock_seconds)/60
    clock_seconds = mod(clock_seconds,60.0)
    clock_hours = clock_hours+clock_minutes/60
    clock_minutes = mod(clock_minutes,60)
   
  end subroutine update_wallclock
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Return number of seconds wall clock since start of program
  real function timestamp_wallclock ()
    
    call system_clock(end_wallclock_count,counts_per_sec)
    timestamp_wallclock = real(end_wallclock_count-absolute_wallclock_count)/real(counts_per_sec)

  end function timestamp_wallclock
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Report all the clocks
  subroutine report_clocks
    
    call report_cpuclock()
    call report_wallclock()
    
  end subroutine report_clocks
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Report CPU clock to log file
  subroutine report_cpuclock
    
    call update_cpuclock()
    if (rank .eq. 0) then
       write(logf,*) "CPU time: ", cpu_hour, ' hours', cpu_minute, ' minutes', &
                                   cpu_second,' seconds.'
       write(*,*) "CPU time: ", cpu_hour, ' hours', cpu_minute,' minutes', &
                                cpu_second, ' seconds.'
    endif

  end subroutine report_cpuclock
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Report wall clock to log file
  subroutine report_wallclock
    
    call update_wallclock()
    if (rank .eq. 0) then
       write(logf,*) "Wall clock time: ", clock_hours, ' hours', &
                     clock_minutes, ' minutes', clock_seconds, ' seconds.'
       write(*,*) "Wall clock time: ",clock_hours, ' hours', &
                  clock_minutes, ' minutes' ,clock_seconds, ' seconds.'
    endif

  end subroutine report_wallclock

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Report simulation time, next output time and end time
  subroutine report_simulation_time

    if (rank .eq. 0) then
       write(logf,*) "Simulation time = ",sim_time/YEAR,"years"
       write(logf,*) "Next output time = ",next_output_time/YEAR,"years"
       write(logf,*) "End time = ",end_time/YEAR,"years"
    endif

  end subroutine report_simulation_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module clocks
