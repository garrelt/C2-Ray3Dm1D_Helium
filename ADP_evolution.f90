module ADP_evolution

  use precision, only: dp
  use my_mpi 
  use file_admin, only: logf, timefile, iterdump, results_dir, dump_dir
  use clocks, only: timestamp_wallclock
  use c2ray_parameters, only: convergence_fraction
  use parameter, only: mesh, ADP_size_of_m
  use array, only: temperature_array
  use array, only: xHI_array,xHII_array
  use array, only: xHeI_array,xHeII_array,xHeIII_array
  use array, only: temperature_array
  use array, only: intermediate_average_xHI,intermediate_average_xHII
  use array, only: intermediate_end_xHI,intermediate_end_xHII
  use array, only: intermediate_average_xHeI,intermediate_average_xHeII,intermediate_average_xHeIII
  use array, only: intermediate_end_xHeI,intermediate_end_xHeII,intermediate_end_xHeIII
  use array, only: intermediate_end_temperature,intermediate_average_temperature
  use array, only: ADP_m_to_i, ADP_m_to_j, ADP_m_to_k
  use array, only: ADP_MPI_intermediate_average_xHI, &
                   ADP_MPI_intermediate_average_xHII, &
                   ADP_MPI_intermediate_average_xHeI, &
                   ADP_MPI_intermediate_average_xHeII, &
                   ADP_MPI_intermediate_average_xHeIII, &
                   ADP_MPI_intermediate_end_xHI, &
                   ADP_MPI_intermediate_end_xHII, &
                   ADP_MPI_intermediate_end_xHeI, &
                   ADP_MPI_intermediate_end_xHeII, &
                   ADP_MPI_intermediate_end_xHeIII, &
                   ADP_MPI_intermediate_end_temperature, &
                   ADP_MPI_intermediate_average_temperature
  use array, only: photoionization_HI_array
  use array, only: photoionization_HeI_array
  use array, only: photoionization_HeII_array
  use array, only: heating_array
  use array, only: buffer
  use array, only: intermediate_end_temperature
  use parameter, only: number_of_source
  use local_evolution, only: evolution_of_cell
  use parameter, only: iteration_counter,convergence_failure,conv_criterion
  use ADP_convergence, only: ADP_convergence_test
  use ADP_copy, only: ADP_copy_final_answer


  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ADP_evolve(dt)

    implicit none

    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    integer :: convergence_failure_of_processor
    integer :: i,j,k,m  
    integer,dimension(1:3) :: pos
    integer :: number_of_job_per_processor, mod_number
    integer :: first_job, last_job
    integer :: temp_convergence
    integer :: mympierror

    convergence_failure = 0 
    convergence_failure_of_processor = 0

    number_of_job_per_processor = ADP_size_of_m/nprocs
    mod_number = mod(ADP_size_of_m,nprocs)

    first_job = rank*number_of_job_per_processor+1 + min(rank,mod_number)
    last_job = (rank+1)*number_of_job_per_processor + min(rank+1,mod_number)

#ifdef MPI

      call MPI_Barrier(MPI_COMM_NEW, mympierror)

#endif

!$OMP PARALLEL PRIVATE(i,j,k,pos)

!$OMP DO SCHEDULE (DYNAMIC,1) 

    do m = first_job,last_job
      i = ADP_m_to_i(m)
      j = ADP_m_to_j(m)
      k = ADP_m_to_k(m)
      pos = (/ i,j,k /)
      ADP_MPI_intermediate_average_xHI(i,j,k) = intermediate_average_xHI(i,j,k)
      ADP_MPI_intermediate_average_xHII(i,j,k) = intermediate_average_xHII(i,j,k)
      ADP_MPI_intermediate_average_xHeI(i,j,k) = intermediate_average_xHeI(i,j,k)
      ADP_MPI_intermediate_average_xHeII(i,j,k) = intermediate_average_xHeII(i,j,k)
      ADP_MPI_intermediate_average_xHeIII(i,j,k) = intermediate_average_xHeIII(i,j,k)
      ADP_MPI_intermediate_average_temperature(i,j,k) = intermediate_average_temperature(i,j,k)
      ADP_MPI_intermediate_end_xHI(i,j,k) = intermediate_end_xHI(i,j,k)
      ADP_MPI_intermediate_end_xHII(i,j,k) = intermediate_end_xHII(i,j,k)
      ADP_MPI_intermediate_end_xHeI(i,j,k) = intermediate_end_xHeI(i,j,k)
      ADP_MPI_intermediate_end_xHeII(i,j,k) = intermediate_end_xHeII(i,j,k)
      ADP_MPI_intermediate_end_xHeIII(i,j,k) = intermediate_end_xHeIII(i,j,k)
      ADP_MPI_intermediate_end_temperature(i,j,k) = intermediate_end_temperature(i,j,k)
      call evolution_of_cell(dt,pos,convergence_failure_of_processor,'A')
    enddo  

!$OMP END DO

!$OMP END PARALLEL

#ifdef MPI

    call MPI_Barrier(MPI_COMM_NEW, mympierror)

    call MPI_ALLREDUCE(convergence_failure_of_processor,temp_convergence,1, MPI_INTEGER, MPI_SUM, MPI_COMM_NEW, mympierror)        
    convergence_failure = temp_convergence

    call MPI_ALLREDUCE(ADP_MPI_intermediate_average_xHI, buffer, mesh*mesh*mesh, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror) 
    do m = 1,ADP_size_of_m
      intermediate_average_xHI(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m)) = buffer(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m))
    enddo       

    call MPI_ALLREDUCE(ADP_MPI_intermediate_average_xHII, buffer, mesh*mesh*mesh, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_NEW, mympierror)        
    do m = 1,ADP_size_of_m
      intermediate_average_xHII(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m)) = buffer(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m))
    enddo       

    call MPI_ALLREDUCE(ADP_MPI_intermediate_average_xHeI, buffer, mesh*mesh*mesh, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror) 
    do m = 1,ADP_size_of_m
      intermediate_average_xHeI(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m)) = buffer(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m))
    enddo 

    call MPI_ALLREDUCE(ADP_MPI_intermediate_average_xHeII, buffer, mesh*mesh*mesh, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror) 
    do m = 1,ADP_size_of_m
      intermediate_average_xHeII(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m)) = buffer(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m))
    enddo 

    call MPI_ALLREDUCE(ADP_MPI_intermediate_average_xHeIII, buffer, mesh*mesh*mesh, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror) 
    do m = 1,ADP_size_of_m
      intermediate_average_xHeIII(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m)) = buffer(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m))
    enddo 

    call MPI_ALLREDUCE(ADP_MPI_intermediate_average_temperature,buffer,mesh*mesh*mesh, &
                       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_NEW,mympierror)        
    do m = 1,ADP_size_of_m
     intermediate_average_temperature(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m)) = buffer(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m))
    enddo       

    call MPI_ALLREDUCE(ADP_MPI_intermediate_end_xHI, buffer, mesh*mesh*mesh, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)        
    do m = 1,ADP_size_of_m
      intermediate_end_xHI(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m)) = buffer(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m))
    enddo       


    call MPI_ALLREDUCE(ADP_MPI_intermediate_end_xHII, buffer, mesh*mesh*mesh, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)        
    do m = 1,ADP_size_of_m
      intermediate_end_xHII(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m)) = buffer(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m))
    enddo       

    call MPI_ALLREDUCE(ADP_MPI_intermediate_end_xHeI, buffer, mesh*mesh*mesh, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)        
    do m = 1,ADP_size_of_m
      intermediate_end_xHeI(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m)) = buffer(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m))
    enddo    

    call MPI_ALLREDUCE(ADP_MPI_intermediate_end_xHeII, buffer, mesh*mesh*mesh, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)        
    do m = 1,ADP_size_of_m
      intermediate_end_xHeII(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m)) = buffer(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m))
    enddo    

    call MPI_ALLREDUCE(ADP_MPI_intermediate_end_xHeIII, buffer, mesh*mesh*mesh, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)        
    do m = 1,ADP_size_of_m
      intermediate_end_xHeIII(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m)) = buffer(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m))
    enddo    

    call MPI_ALLREDUCE(ADP_MPI_intermediate_end_temperature,buffer,mesh*mesh*mesh, &
                       MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_NEW,mympierror)        
    do m = 1,ADP_size_of_m
      intermediate_end_temperature(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m)) = buffer(ADP_m_to_i(m),ADP_m_to_j(m),ADP_m_to_k(m))
    enddo       

    call MPI_Barrier(MPI_COMM_NEW, mympierror)

#endif


  end subroutine ADP_evolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ADP_evolution
