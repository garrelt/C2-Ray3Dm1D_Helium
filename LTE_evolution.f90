module LTE_evolution

  use precision, only: dp 
  use my_mpi
  use array, only: global_LTE_array, &
                   LTE_m_to_i,LTE_m_to_j,LTE_m_to_k, &
                   LTE_MPI_intermediate_average_xHI,&
                   LTE_MPI_intermediate_average_xHII,&
                   LTE_MPI_intermediate_average_xHeI,&
                   LTE_MPI_intermediate_average_xHeII,&
                   LTE_MPI_intermediate_average_xHeIII,&
                   LTE_MPI_intermediate_average_temperature, &
                   LTE_MPI_intermediate_end_xHI,&
                   LTE_MPI_intermediate_end_xHII,&
                   LTE_MPI_intermediate_end_xHeI,&
                   LTE_MPI_intermediate_end_xHeII,&
                   LTE_MPI_intermediate_end_xHeIII,&
                   LTE_MPI_intermediate_end_temperature, &
                   xHI_array, xHII_array, temperature_array
  use parameter, only: mesh, LTE_size_of_m
  use array, only: intermediate_average_xHI,intermediate_average_xHII
  use array, only: intermediate_end_xHI,intermediate_end_xHII
  use array, only: intermediate_average_xHeI,intermediate_average_xHeII,intermediate_average_xHeIII
  use array, only: intermediate_end_xHeI,intermediate_end_xHeII,intermediate_end_xHeIII
  use array, only: intermediate_end_temperature,intermediate_average_temperature
  use array, only: buffer,global_evolved_LTE_array
  use parameter, only: LTE_iteration_counter,LTE_convergence_failure
  use local_evolution, only: evolution_of_cell

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LTE_evolve()

    implicit none

    integer :: i,j,k,m 
    integer,dimension(1:3) :: pos
    real(kind=dp) :: dt 
    integer :: convergence_failure_of_processor
    integer :: number_of_job_per_processor, mod_number
    integer :: first_job, last_job
    integer :: dummy
    integer :: mympierror
    integer :: temp_convergence


    dummy = -1

      number_of_job_per_processor = LTE_size_of_m/nprocs
      mod_number = mod(LTE_size_of_m,nprocs)

      first_job = rank*number_of_job_per_processor+1 + min(rank,mod_number)
      last_job = (rank+1)*number_of_job_per_processor + min(rank+1,mod_number)

#ifdef MPI

      call MPI_Barrier(MPI_COMM_NEW, mympierror)

#endif

!$OMP PARALLEL PRIVATE(i,j,k,pos,dt)

!$OMP DO SCHEDULE (DYNAMIC,1) 

      do m = first_job,last_job

        i = LTE_m_to_i(m)
        j = LTE_m_to_j(m)
        k = LTE_m_to_k(m)
        pos = (/ i,j,k  /)
        dt = global_LTE_array(i,j,k)

        call evolution_of_cell(dt,pos,dummy,'L')

      enddo

!$OMP END DO

!$OMP END PARALLEL

#ifdef MPI

      call MPI_Barrier(MPI_COMM_NEW, mympierror) 

      call MPI_ALLREDUCE(LTE_MPI_intermediate_end_xHI, buffer, mesh*mesh*mesh, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_NEW, mympierror)        
      do m = 1,LTE_size_of_m
        intermediate_end_xHI(LTE_m_to_i(m),LTE_m_to_j(m),LTE_m_to_k(m)) = buffer(LTE_m_to_i(m),LTE_m_to_j(m),LTE_m_to_k(m))
      enddo       

      call MPI_ALLREDUCE(LTE_MPI_intermediate_end_xHII, buffer, mesh*mesh*mesh, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_NEW, mympierror)        
      do m = 1,LTE_size_of_m
        intermediate_end_xHII(LTE_m_to_i(m),LTE_m_to_j(m),LTE_m_to_k(m)) = buffer(LTE_m_to_i(m),LTE_m_to_j(m),LTE_m_to_k(m))
      enddo       

      call MPI_ALLREDUCE(LTE_MPI_intermediate_end_xHeI, buffer, mesh*mesh*mesh, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_NEW, mympierror)        
      do m = 1,LTE_size_of_m
        intermediate_end_xHeI(LTE_m_to_i(m),LTE_m_to_j(m),LTE_m_to_k(m)) = buffer(LTE_m_to_i(m),LTE_m_to_j(m),LTE_m_to_k(m))
      enddo       

      call MPI_ALLREDUCE(LTE_MPI_intermediate_end_xHeII, buffer, mesh*mesh*mesh, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_NEW, mympierror)        
      do m = 1,LTE_size_of_m
        intermediate_end_xHeII(LTE_m_to_i(m),LTE_m_to_j(m),LTE_m_to_k(m)) = buffer(LTE_m_to_i(m),LTE_m_to_j(m),LTE_m_to_k(m))
      enddo       

      call MPI_ALLREDUCE(LTE_MPI_intermediate_end_xHeIII, buffer, mesh*mesh*mesh, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_NEW, mympierror)        
      do m = 1,LTE_size_of_m
        intermediate_end_xHeIII(LTE_m_to_i(m),LTE_m_to_j(m),LTE_m_to_k(m)) = buffer(LTE_m_to_i(m),LTE_m_to_j(m),LTE_m_to_k(m))
      enddo       

      call MPI_ALLREDUCE(LTE_MPI_intermediate_end_temperature,buffer,mesh*mesh*mesh,MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_NEW,mympierror)        
      do m = 1,LTE_size_of_m
  intermediate_end_temperature(LTE_m_to_i(m),LTE_m_to_j(m),LTE_m_to_k(m)) = buffer(LTE_m_to_i(m),LTE_m_to_j(m),LTE_m_to_k(m))
      enddo       

      call MPI_Barrier(MPI_COMM_NEW, mympierror)

#endif


  end subroutine LTE_evolve

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module LTE_evolution
