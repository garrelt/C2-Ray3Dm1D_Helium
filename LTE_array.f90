module LTE_array

  use my_mpi
  use file_admin, only: logf
  use parameter, only: mesh, number_of_source, LTE_size_of_m
  use array, only: source_class_index,photoionization_HI_array, &
                   photoionization_HeI_array,photoionization_HeII_array,&
                   heating_array, global_LTE_array, equivalent_class_index_array, &
                   global_evolved_LTE_array, &
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
                   LTE_m_to_i, LTE_m_to_j, LTE_m_to_k, &
                   LTE_photon_loss_array
  use sourceprops, only: number_of_source
  use parameter, only: mesh,conv_criterion
  use parameter, only:  LTE_iteration_counter,LTE_convergence_failure
  use c2ray_parameters, only: convergence_fraction
  use array, only: intermediate_average_xHI,intermediate_average_xHII, &
                   intermediate_end_xHI, intermediate_end_xHII, &
                   intermediate_average_xHeI,intermediate_average_xHeI, &
                   intermediate_average_xHeI, intermediate_end_xHeI, &
                   intermediate_end_xHeII, intermediate_end_xHeIII, &
                   intermediate_end_temperature, intermediate_average_temperature
  use array, only: xHI_array,xHII_array,xHeI_array,xHeII_array,xHeIII_array,temperature_array,&
                    global_evolved_LTE_array
  use array, only: photoionization_HI_array, &
                   photoionization_HeI_array, &
                   photoionization_HeII_array, &
                   heating_array

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LTE_array_initialization_1()

    implicit none

    integer :: checkstat

    if (.not.allocated(source_class_index)) then
      allocate(source_class_index(1:number_of_source),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'source_class_index checkstat',checkstat
    endif 

    !global_LTE_array = 0
    !photoionization_HI_array = 0.0
    !photoionization_HeI_array = 0.0
    !photoionization_HeII_array = 0.0
    !heating_array = 0.0
    equivalent_class_index_array = 0
    source_class_index = -1

    LTE_photon_loss_array = 0.0

  end subroutine LTE_array_initialization_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LTE_array_initialization_2()

    implicit none

    integer :: i,j,k
    integer :: counter
    integer :: checkstat

      if (.not.allocated(LTE_MPI_intermediate_average_xHI)) then
        allocate(LTE_MPI_intermediate_average_xHI(1:mesh,1:mesh,1:mesh),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_average_xHI checkstat',checkstat
      endif 

      if (.not.allocated(LTE_MPI_intermediate_average_xHII)) then
        allocate(LTE_MPI_intermediate_average_xHII(1:mesh,1:mesh,1:mesh),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_average_xHII checkstat',checkstat
      endif 

      if (.not.allocated(LTE_MPI_intermediate_average_xHeI)) then
        allocate(LTE_MPI_intermediate_average_xHeI(1:mesh,1:mesh,1:mesh),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_average_xHeI checkstat',checkstat
      endif 

      if (.not.allocated(LTE_MPI_intermediate_average_xHeII)) then
        allocate(LTE_MPI_intermediate_average_xHeII(1:mesh,1:mesh,1:mesh),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_average_xHeII checkstat',checkstat
      endif 

      if (.not.allocated(LTE_MPI_intermediate_average_xHeIII)) then
        allocate(LTE_MPI_intermediate_average_xHeIII(1:mesh,1:mesh,1:mesh),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_average_xHeIII checkstat',checkstat
       endif 

      if (.not.allocated(LTE_MPI_intermediate_average_temperature)) then
        allocate(LTE_MPI_intermediate_average_temperature(1:mesh,1:mesh,1:mesh),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_average_temperature checkstat',checkstat
      endif 

      if (.not.allocated(LTE_MPI_intermediate_end_xHI)) then
        allocate(LTE_MPI_intermediate_end_xHI(1:mesh,1:mesh,1:mesh),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_end_xHI checkstat',checkstat
      endif 

      if (.not.allocated(LTE_MPI_intermediate_end_xHII)) then
        allocate(LTE_MPI_intermediate_end_xHII(1:mesh,1:mesh,1:mesh),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_end_xHII checkstat',checkstat
      endif 

      if (.not.allocated(LTE_MPI_intermediate_end_xHeI)) then
        allocate(LTE_MPI_intermediate_end_xHeI(1:mesh,1:mesh,1:mesh),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_end_xHeI checkstat',checkstat
      endif 

      if (.not.allocated(LTE_MPI_intermediate_end_xHeII)) then
        allocate(LTE_MPI_intermediate_end_xHeII(1:mesh,1:mesh,1:mesh),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_end_xHeII checkstat',checkstat
      endif 

      if (.not.allocated(LTE_MPI_intermediate_end_xHeIII)) then
        allocate(LTE_MPI_intermediate_end_xHeIII(1:mesh,1:mesh,1:mesh),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_end_xHeIII checkstat',checkstat
      endif 

      if (.not.allocated(LTE_MPI_intermediate_end_temperature)) then
        allocate(LTE_MPI_intermediate_end_temperature(1:mesh,1:mesh,1:mesh),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_end_temperature checkstat',checkstat
      endif 

!averge no use
      LTE_MPI_intermediate_average_xHI = 0.0
      LTE_MPI_intermediate_average_xHII = 0.0
      LTE_MPI_intermediate_average_xHeI = 0.0
      LTE_MPI_intermediate_average_xHeII = 0.0
      LTE_MPI_intermediate_average_xHeIII = 0.0
      LTE_MPI_intermediate_average_temperature = 0.0
      LTE_MPI_intermediate_end_xHI = 0.0
      LTE_MPI_intermediate_end_xHII = 0.0
      LTE_MPI_intermediate_end_xHeI = 0.0
      LTE_MPI_intermediate_end_xHeII = 0.0
      LTE_MPI_intermediate_end_xHeIII = 0.0
      LTE_MPI_intermediate_end_temperature = 0.0

      !LTE_size_of_m = sum(global_LTE_array)-sum(global_evolved_LTE_array)
      LTE_size_of_m = 0
      do k = 1,mesh
        do j = 1,mesh
          do i = 1,mesh
            if (global_LTE_array(i,j,k).gt.0) LTE_size_of_m = LTE_size_of_m+1
          enddo
        enddo
      enddo

!write(*,*) 'in LTE_array global_LTE_array(15,15,15) is ',global_LTE_array(15,15,15)
write(*,*)'rank',rank,'report LTE_size_of_m is',LTE_size_of_m

      if (.not.allocated(LTE_m_to_i)) then
        allocate(LTE_m_to_i(1:LTE_size_of_m),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_m_to_i checkstat',checkstat
      endif 

      if (.not.allocated(LTE_m_to_j)) then
        allocate(LTE_m_to_j(1:LTE_size_of_m),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_m_to_j checkstat',checkstat
      endif 

      if (.not.allocated(LTE_m_to_k)) then
        allocate(LTE_m_to_k(1:LTE_size_of_m),STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_m_to_k checkstat',checkstat
      endif 



      counter = 1

      do k = 1,mesh
        do j = 1,mesh
          do i = 1,mesh
            if (global_LTE_array(i,j,k).gt.0) then 
              LTE_m_to_i(counter) = i
              LTE_m_to_j(counter) = j
              LTE_m_to_k(counter) = k
              counter = counter+1
            endif
          enddo
        enddo
      enddo


    if (rank .eq. 0) then
      !write(*,*) 'LTE ',counter-1
      write(logf,*) "LTE evolution cells = ",counter-1
      flush(logf) 
    endif

  end subroutine LTE_array_initialization_2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LTE_reset_evolved_LTE_field()	

    implicit none

    global_LTE_array = -1
    global_evolved_LTE_array = 0

  end subroutine LTE_reset_evolved_LTE_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LTE_array_destruction2()

    implicit none
    integer :: checkstat

      !write(logf,*) "deallocate(LTE_MPI_intermediate_average_xHI)"
      !flush(logf) 
      if (allocated(LTE_MPI_intermediate_average_xHI)) then
        deallocate(LTE_MPI_intermediate_average_xHI,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_average_xHI checkstat',checkstat
      endif

      !write(logf,*) "deallocate(LTE_MPI_intermediate_average_xHII)"
      !flush(logf) 
      if (allocated(LTE_MPI_intermediate_average_xHII)) then 
        deallocate(LTE_MPI_intermediate_average_xHII,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_average_xHII checkstat',checkstat
      endif

      !write(logf,*) "deallocate(LTE_MPI_intermediate_average_xHeI)"
      !flush(logf)  
      if (allocated(LTE_MPI_intermediate_average_xHeI)) then 
        deallocate(LTE_MPI_intermediate_average_xHeI,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_average_xHeI checkstat',checkstat
      endif

      !write(logf,*) "eallocate(LTE_MPI_intermediate_average_xHeII)"
      !flush(logf)  
      if (allocated(LTE_MPI_intermediate_average_xHeII)) then 
        deallocate(LTE_MPI_intermediate_average_xHeII,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_average_xHeII checkstat',checkstat
      endif

      !write(logf,*) "deallocate(LTE_MPI_intermediate_average_xHeIII)"
      !flush(logf)  
      if (allocated(LTE_MPI_intermediate_average_xHeIII)) then 
        deallocate(LTE_MPI_intermediate_average_xHeIII,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_average_xHeIII checkstat',checkstat
      endif

      !write(logf,*) "deallocate(LTE_MPI_intermediate_average_temperature)"
      !flush(logf)  
      if (allocated(LTE_MPI_intermediate_average_temperature)) then 
        deallocate(LTE_MPI_intermediate_average_temperature,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_average_temperature checkstat',checkstat
      endif

      !write(logf,*) "deallocate(LTE_MPI_intermediate_end_xHI)"
      !flush(logf)  
      if (allocated(LTE_MPI_intermediate_end_xHI)) then 
        deallocate(LTE_MPI_intermediate_end_xHI,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_end_xHI checkstat',checkstat
      endif

      !write(logf,*) "deallocate(LTE_MPI_intermediate_end_xHII)"
      !flush(logf)  
      if (allocated(LTE_MPI_intermediate_end_xHII)) then 
        deallocate(LTE_MPI_intermediate_end_xHII,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_end_xHII checkstat',checkstat
      endif

      !write(logf,*) "deallocate(LTE_MPI_intermediate_end_xHeI)"
      !flush(logf)  
      if (allocated(LTE_MPI_intermediate_end_xHeI)) then 
        deallocate(LTE_MPI_intermediate_end_xHeI,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_end_xHeI checkstat',checkstat
      endif

      !write(logf,*) "deallocate(LTE_MPI_intermediate_end_xHeII)"
      !flush(logf)  
      if (allocated(LTE_MPI_intermediate_end_xHeII)) then 
        deallocate(LTE_MPI_intermediate_end_xHeII,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_end_xHeII checkstat',checkstat
      endif

      !write(logf,*) "deallocate(LTE_MPI_intermediate_end_xHeIII)"
      !flush(logf)  
      if (allocated(LTE_MPI_intermediate_end_xHeIII)) then 
        deallocate(LTE_MPI_intermediate_end_xHeIII,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_end_xHeIII checkstat',checkstat
      endif

      !write(logf,*) "eallocate(LTE_MPI_intermediate_end_temperature)"
      !flush(logf)  
      if (allocated(LTE_MPI_intermediate_end_temperature)) then 
        deallocate(LTE_MPI_intermediate_end_temperature,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_MPI_intermediate_end_temperature checkstat',checkstat
      endif

      if (allocated(LTE_m_to_i)) then 
        deallocate(LTE_m_to_i,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_m_to_i checkstat',checkstat
      endif

      if (allocated(LTE_m_to_j)) then 
        deallocate(LTE_m_to_j,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_m_to_j checkstat',checkstat
      endif

      if (allocated(LTE_m_to_k)) then 
        deallocate(LTE_m_to_k,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'LTE_m_to_k checkstat',checkstat
      endif

    if (rank .eq. 0) then
      write(logf,*) "Finish LTE evolution"
      flush(logf) 
    endif

  end subroutine LTE_array_destruction2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LTE_array_destruction1()

    implicit none
    integer :: checkstat

    if (allocated(source_class_index)) then 
      deallocate(source_class_index,STAT=checkstat)
      if (checkstat.gt.0) write(*,*) 'rank',rank,'source_class_index checkstat',checkstat
    endif

global_LTE_array = -1.0
    !photoionization_HI_array = 0.0
    !photoionization_HeI_array = 0.0 
    !photoionization_HeII_array = 0.0  
    !heating_array = 0.0    
!write(*,*) 'after LTE, the time is ',global_LTE_array(5,5,5)
  end subroutine LTE_array_destruction1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module LTE_array
