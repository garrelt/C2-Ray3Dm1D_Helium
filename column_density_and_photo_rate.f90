module column_density_and_photo_rate

  use my_mpi
  use parameter, only: actual_dt,LTE_size_of_m
  use sourceprops, only: number_of_source, srcpos
  use global_column_density, only: global_get_column_density
  use global_photo_rate, only: global_get_photo_rate
  use array, only: global_LTE_array,LTE_m_to_i

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine master_distributing_job_to_slaves()

    implicit none

    integer :: source_index
    integer :: sources_done,whomax,who,answer
    integer,dimension(:),allocatable :: counts
    integer :: mympierror

#ifdef MPI

    sources_done = 0         
    source_index = 0
    allocate(counts(0:nprocs-1))
  
    whomax = min(number_of_source,nprocs-1)

    do who = 1,whomax

      if (source_index .le. number_of_source) then
        source_index = source_index+1
        call MPI_Send (source_index, 1, MPI_INTEGER, who, 1, MPI_COMM_NEW, mympierror)
      endif

    enddo
    
    do while (sources_done .lt. number_of_source)

      call MPI_Recv (answer,     & ! address of receive buffer
           1,		   & ! number of items to receive
           MPI_INTEGER,	   & ! type of data
           MPI_ANY_SOURCE,  & ! can receive from any other
           1,		   & ! tag
           MPI_COMM_NEW,	   & ! communicator
           mympi_status,	   & ! status
           mympierror)
       
       who = mympi_status(MPI_SOURCE) ! find out who sent us the answer
       sources_done = sources_done+1 ! and the number of sources done
       
       if (source_index .lt. number_of_source) then

          source_index = source_index+1
          call MPI_Send (source_index, 1, MPI_INTEGER, &
               who,		&	
               1,		&	
               MPI_COMM_NEW, &
               mympierror)
       endif

    enddo
    
    do who = 1,nprocs-1

       call MPI_Send (0, 1, MPI_INTEGER, &
            who,			  &
            2,			  & ! tag 
            MPI_COMM_NEW,	          &
            mympierror)
       
       ! the slave will send to master the number of calculations
       ! that have been performed. 
       ! We put this number in the counts array.
       
       call MPI_Recv (counts(who), & ! address of receive buffer
            1,                & ! number of items to receive
            MPI_INTEGER,      & ! type of data 
            who,              & ! receive from process who 
            7,                & ! tag 
            MPI_COMM_NEW,     & ! communicator 
            mympi_status,     & ! status
            mympierror)
    enddo

    deallocate(counts)

#endif

  end subroutine master_distributing_job_to_slaves

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine slave_get_column_density_and_photo_rate(case_of_evolution,ADP_help_do_LTE_ray_tracing) 

    implicit none

    character,intent(in) :: case_of_evolution
    logical,intent(in) :: ADP_help_do_LTE_ray_tracing
    integer :: local_count
    integer :: source_index
    integer :: mympierror

#ifdef MPI

    local_count = 0
    call MPI_Recv (source_index,  & ! address of receive buffer
         1,    & ! number of items to receive
         MPI_INTEGER,  & ! type of data
         0,		  & ! can receive from master only
         MPI_ANY_TAG,  & ! can expect two values, so
         ! we use the wildcard MPI_ANY_TAG 
         ! here
         MPI_COMM_NEW, & ! communicator
         mympi_status, & ! status
         mympierror)
    
    ! if tag equals 2, then skip the calculations   
    if (mympi_status(MPI_TAG) /= 2) then

      do 

        if (case_of_evolution.eq.'L') then

          ! make sure the source point is in LTE
          if (global_LTE_array(srcpos(1,source_index),srcpos(2,source_index),srcpos(3,source_index)).eq.1) then ! maybe redundant
            call global_get_column_density(source_index,'L',ADP_help_do_LTE_ray_tracing) 
            call global_get_photo_rate(source_index,'L',ADP_help_do_LTE_ray_tracing) 
          endif

        elseif (case_of_evolution.eq.'A') then

          ! Do the source at hand
          call global_get_column_density(source_index,'A',ADP_help_do_LTE_ray_tracing)
          call global_get_photo_rate(source_index,'A',ADP_help_do_LTE_ray_tracing) 

        endif

        ! Update local counter
        local_count = local_count+1
          
        ! Send 'answer'
        call MPI_Send (local_count, 1,  & ! sending one int 
             MPI_INTEGER, 0, & ! to master
             1,              & ! tag
             MPI_COMM_NEW,   & ! communicator
             mympierror)
          
        ! Receive new source number
        call MPI_Recv (source_index,     & ! address of receive buffer
             1,            & ! number of items to receive
             MPI_INTEGER,  & ! type of data
             0,            & ! can receive from master only
             MPI_ANY_TAG,  & !  can expect two values, so
             !  we use the wildcard MPI_ANY_TAG 
             !  here
             MPI_COMM_NEW, & ! communicator
             mympi_status, & ! status
             mympierror)
          
        ! leave this loop if tag equals 2
        if (mympi_status(MPI_TAG) == 2) then
          exit 
        endif
      enddo
    endif
    
    call MPI_Send (local_count,  &
         1,           & 
         MPI_INTEGER, & ! sending one int
         0,           & ! to master
         7,           & ! tag
         MPI_COMM_NEW,& ! communicator
         mympierror)

#endif

  end subroutine slave_get_column_density_and_photo_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine processor_get_column_density_and_photo_rate(case_of_evolution,ADP_help_do_LTE_ray_tracing)
    
    implicit none

    character,intent(in) :: case_of_evolution
    logical,intent(in) :: ADP_help_do_LTE_ray_tracing
    integer :: ns,source_index
    integer :: mympierror

#ifdef MPI

      call MPI_Barrier(MPI_COMM_NEW, mympierror)

#endif

    if (case_of_evolution.eq.'L') then 
        
      do ns = 1+rank,number_of_source,nprocs
        if (global_LTE_array(srcpos(1,ns),srcpos(2,ns),srcpos(3,ns)).ge.0) then !maybe redundant, need to check
          call global_get_column_density(ns,'L',ADP_help_do_LTE_ray_tracing) 
          call global_get_photo_rate(ns,'L',ADP_help_do_LTE_ray_tracing) 
        endif

      enddo

    elseif (case_of_evolution.eq.'A') then

      do source_index = 1+rank,number_of_source,nprocs

        call global_get_column_density(source_index,'A',ADP_help_do_LTE_ray_tracing)
        call global_get_photo_rate(source_index,'A',ADP_help_do_LTE_ray_tracing) 

      enddo

    endif

#ifdef MPI

      call MPI_Barrier(MPI_COMM_NEW, mympierror)

#endif

  end subroutine processor_get_column_density_and_photo_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module column_density_and_photo_rate
