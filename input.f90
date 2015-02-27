module input

  use my_mpi
  use file_admin, only: stdinput,logf,flag_for_file_input,file_input
  use parameter, only: input_source_temperature,input_IGM_temperature,&
                       number_timesteps,number_outputs
  use c2ray_parameters, only: T_eff_nominal


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine setup_input()

    implicit none

    character(len=512) :: inputfile

    if (rank .eq. 0) then
      write(logf,*) "Beginning of input setup"
      flush(logf)
      if (COMMAND_ARGUMENT_COUNT () .gt. 0) then
        call GET_COMMAND_ARGUMENT(1,inputfile)
        write(logf,*) "Reading input from ", trim(adjustl(inputfile))
        open(unit=stdinput,file=inputfile)
        call flag_for_file_input(.true.)
      else
        write(logf,*) "Reading input from command line"
      endif
      write(logf,*) "End of input setup"

      write(logf,*)      
      flush(logf)
    endif

  end subroutine setup_input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine input_data()

    implicit none

    integer :: mympierror

    if (rank .eq. 0) then

      write(logf,*) "Beginning of input data"

      ! Ask for source temperature
      if (.not.file_input) then
        write(*,"(A,$)") "Enter source temperature (K): "
      endif
      read(stdinput,*) input_source_temperature

      if (input_source_temperature.lt.2000 .or. input_source_temperature.gt.200000) then
        write(*,*) 'Source temperature is out of range (2000K < T < 200000K)'
        write(*,*) 'The sources will be set to 50000K'
        input_source_temperature = T_eff_nominal
      endif
      write(logf,*) "Source temperature = ",input_source_temperature,"K"

      ! Ask for IGM temperature
      if (.not.file_input) then
        write(*,"(A,$)") "Enter initial IGM temperature (K): "
      endif
      read(stdinput,*) input_IGM_temperature
      write(logf,*) "Initial IGM temperature = ",input_IGM_temperature,"K"

      ! Ask for number of time steps
      if (.not.file_input) then
        write(*,'(A,$)') 'Enter number of time steps between slices: '
      endif
      read(stdinput,*) number_timesteps
      write(logf,*) "Number of time steps between slices = ", number_timesteps

      ! Ask for interval between outputs
      if (.not.file_input) then
        write(*,'(A,$)') 'Enter number of outputs between slices: '
      endif
      read(stdinput,*) number_outputs
      write(logf,*) "Number of outputs between slices = ", number_outputs

      write(logf,*) "End of input data"
      write(logf,*) 
      flush(logf)

    endif

#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(input_source_temperature,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(input_IGM_temperature,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(number_timesteps,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(number_outputs,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

  end subroutine input_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module input

