!>
!! \brief This module handles the time variables
!!
!! Module for C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2008-05-29 (no date before)
!<
module times

  ! This module handles the time variables

  use precision, only: dp
  use my_mpi
  use file_admin, only: stdinput
  use astroconstants, only: YEAR

  implicit none

!  private

  real(kind=dp),public :: end_time !< End time of simulation
  real(kind=dp),public :: dtlog , dt!< Time step
  real(kind=dp),public :: output_time !< After how many outputs to write an output
  real(kind=dp),dimension(:),allocatable  :: timess
  integer :: i, number_timesteps, outputsteps, j
contains

  ! =======================================================================

  !>
  !!  Initializes the module variables
  !<
  subroutine time_ini ( )
    
    ! Initializes number of time steps per frame (integration and output)

    ! Author: Garrelt Mellema

    ! Date: 20-Aug-2006 (f77: 10-Mar-2004)

    real(kind=dp) :: end_time_yrs,output_time_yrs,startvalue
#ifdef MPI
    integer :: ierror
#endif

    if (rank == 0) then
       ! Ask for end time
       write(*,'(A,$)') 'Enter end time of calculation (years): '
       read(stdinput,*) end_time_yrs

       !  Convert to seconds
       end_time=end_time_yrs*YEAR
       !write(52,*) end_time/YEAR
       ! Ask for number of time steps
       write(*,'(A,$)') 'Enter number of time steps: '
       read(stdinput,*) number_timesteps
      ! write(52,*) number_timesteps
       ! Ask for interval between outputs
       write(*,'(A,$)') 'Enter steps between outputs '
       read(stdinput,*) outputsteps
       outputsteps=int(outputsteps)
       ! Convert to seconds
!       output_time=output_time_yrs*YEAR
!       write(52,*) output_time_yrs
    endif
    allocate(timess(0:number_timesteps))!-1))
#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(end_time,1,MPI_DOUBLEPRECISION,0,&
         MPI_COMM_NEW,ierror)
    call MPI_BCAST(output_time,1,MPI_DOUBLEPRECISION,0,MPI_COMM_NEW,ierror)
#endif
  
!!-----------------------------------------------  
!!   ! Set value of time step  log distribution
!!   write(*,*) 'end_time=',end_time
!   write(52,*) end_time/YEAR
!   startvalue=1*YEAR
!   dtlog=(log10(end_time)-log10(startvalue))/(real(number_timesteps-1))!
!!   write(*,*) 'dtlog=', dtlog
!   write(52,*) outputsteps
!!   write(*,*) 'number_timesteps=', number_timesteps
!   write(52,*) number_timesteps
!   open(unit=78,file='Times.out',status='unknown')
!   timess(0)=0.0_dp
!   do i=1,(number_timesteps)
!     timess(i)=startvalue*10**(i*dtlog)
!   write(78,*)  timess(i)/YEAR
!   enddo
!   close(78)   
!-----------------------------------------------

!------------------------------------------------
!   !Set value of time step     ^5 distribution
!    write(*,*) 'end_time=',end_time/YEAR
!    write(52,*) end_time/YEAR   
!    dtlog=((end_time)**0.2)/real(number_timesteps+1)  !-YEAR
!    write(*,*) 'dtlog=', dtlog
!    write(52,*) outputsteps
!    write(*,*) 'number_timesteps=', number_timesteps
!    write(52,*) number_timesteps
!    do i=0,number_timesteps
!    dt_new(i)=((i+1)*dtlog)**5   !+YEAR
!    enddo
!-----------------------------------------------


!---! Set value of time step     linear distribution
   write(*,*) 'end_time=',end_time/YEAR
   write(52,*) end_time/YEAR   
   dtlog=(end_time/real(number_timesteps)) !-YEAR
   write(*,*) 'dtlog=', dtlog
   write(52,*) outputsteps
   write(*,*) 'number_timesteps=', number_timesteps
   write(52,*) number_timesteps
   do i=0,number_timesteps!-1
   timess(i)=i*dtlog   !+YEAR
   enddo

!---!--------------------------------------------

    dt=end_time/real(number_timesteps)
  end subroutine time_ini

end module times
