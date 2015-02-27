module ADP_array

  use my_mpi
  use file_admin, only: logf
  use array, only: xHI_array,xHII_array,xHeI_array,xHeII_array,xHeIII_array,temperature_array
  use array, only: intermediate_average_xHI,intermediate_average_xHII, &
                   intermediate_end_xHI, intermediate_end_xHII, &
                   intermediate_average_xHeI,intermediate_average_xHeII, &
                   intermediate_average_xHeIII, intermediate_end_xHeI, &
                   intermediate_end_xHeII, intermediate_end_xHeIII, &
                   intermediate_end_temperature, intermediate_average_temperature, &
                   photoionization_HI_array, &
                   photoionization_HeI_array,photoionization_HeII_array,&
                   heating_array
  use parameter, only:  iteration_counter,convergence_failure
  use parameter, only: mesh,conv_criterion, ADP_size_of_m
  use c2ray_parameters, only: convergence_fraction
  use parameter, only: number_of_source, do_LTE
  use array, only: photoionization_HI_array, &
                   photoionization_HeI_array, &
                   photoionization_HeII_array, &
                   heating_array, &
                   ADP_MPI_intermediate_average_xHI, &
                   ADP_MPI_intermediate_average_xHII, &
                   ADP_MPI_intermediate_average_xHeI, &
                   ADP_MPI_intermediate_average_xHeII, &
                   ADP_MPI_intermediate_average_xHeIII, &
                   ADP_MPI_intermediate_average_temperature, &
                   ADP_MPI_intermediate_end_xHI, &
                   ADP_MPI_intermediate_end_xHII, &
                   ADP_MPI_intermediate_end_xHeI, &
                   ADP_MPI_intermediate_end_xHeII, &
                   ADP_MPI_intermediate_end_xHeIII, &
                   ADP_MPI_intermediate_end_temperature, &
                   global_LTE_array, &
                   ADP_m_to_i, ADP_m_to_j, ADP_m_to_k


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ADP_array_initialization_1

    implicit none

    integer :: i,j,k, counter,LTE_size_of_m

      counter = 1

      intermediate_average_xHI(:,:,:) = xHI_array(:,:,:)
      intermediate_average_xHII(:,:,:) = xHII_array(:,:,:)
      intermediate_average_xHeI(:,:,:) = xHeI_array(:,:,:)
      intermediate_average_xHeII(:,:,:) = xHeII_array(:,:,:)
      intermediate_average_xHeIII(:,:,:) = xHeIII_array(:,:,:)
      intermediate_average_temperature(:,:,:) = temperature_array(:,:,:)
      intermediate_end_xHI(:,:,:) = xHI_array(:,:,:)
      intermediate_end_xHII(:,:,:) = xHII_array(:,:,:)
      intermediate_end_xHeI(:,:,:) = xHeI_array(:,:,:)
      intermediate_end_xHeII(:,:,:) = xHeII_array(:,:,:)
      intermediate_end_xHeIII(:,:,:) = xHeIII_array(:,:,:)
      intermediate_end_temperature(:,:,:) = temperature_array(:,:,:)

      allocate(ADP_MPI_intermediate_average_xHI(1:mesh,1:mesh,1:mesh))
      allocate(ADP_MPI_intermediate_average_xHII(1:mesh,1:mesh,1:mesh))
      allocate(ADP_MPI_intermediate_average_xHeI(1:mesh,1:mesh,1:mesh))
      allocate(ADP_MPI_intermediate_average_xHeII(1:mesh,1:mesh,1:mesh))
      allocate(ADP_MPI_intermediate_average_xHeIII(1:mesh,1:mesh,1:mesh))
      allocate(ADP_MPI_intermediate_average_temperature(1:mesh,1:mesh,1:mesh))
      allocate(ADP_MPI_intermediate_end_xHI(1:mesh,1:mesh,1:mesh)) 
      allocate(ADP_MPI_intermediate_end_xHII(1:mesh,1:mesh,1:mesh))
      allocate(ADP_MPI_intermediate_end_xHeI(1:mesh,1:mesh,1:mesh)) 
      allocate(ADP_MPI_intermediate_end_xHeII(1:mesh,1:mesh,1:mesh)) 
      allocate(ADP_MPI_intermediate_end_xHeIII(1:mesh,1:mesh,1:mesh)) 
      allocate(ADP_MPI_intermediate_end_temperature(1:mesh,1:mesh,1:mesh))

      ADP_MPI_intermediate_average_xHI = 0
      ADP_MPI_intermediate_average_xHII = 0
      ADP_MPI_intermediate_average_xHeI = 0
      ADP_MPI_intermediate_average_xHeII = 0
      ADP_MPI_intermediate_average_xHeIII = 0
      ADP_MPI_intermediate_average_temperature = 0
      ADP_MPI_intermediate_end_xHI = 0 
      ADP_MPI_intermediate_end_xHII = 0
      ADP_MPI_intermediate_end_xHeI = 0 
      ADP_MPI_intermediate_end_xHeII = 0 
      ADP_MPI_intermediate_end_xHeIII = 0 
      ADP_MPI_intermediate_end_temperature = 0

      iteration_counter = 0
      convergence_failure = mesh*mesh*mesh 

      LTE_size_of_m = 0
      do k = 1,mesh
        do j = 1,mesh
          do i = 1,mesh
            if (global_LTE_array(i,j,k).gt.0) LTE_size_of_m = LTE_size_of_m+1
          enddo
        enddo
      enddo

if (do_LTE .eqv. .true.) then
      ! what if mesh*mesh*mesh-sum(global_LTE_array) = 0 ?
      allocate(ADP_m_to_i(1:mesh*mesh*mesh-LTE_size_of_m))
      allocate(ADP_m_to_j(1:mesh*mesh*mesh-LTE_size_of_m))
      allocate(ADP_m_to_k(1:mesh*mesh*mesh-LTE_size_of_m))
      ADP_size_of_m =  size(ADP_m_to_i)

      do k = 1,mesh
        do j = 1,mesh
          do i = 1,mesh
            if (global_LTE_array(i,j,k).lt.0) then 
              ADP_m_to_i(counter) = i
              ADP_m_to_j(counter) = j
              ADP_m_to_k(counter) = k
              counter = counter+1
            endif
          enddo
        enddo
      enddo
else

      allocate(ADP_m_to_i(1:mesh*mesh*mesh))
      allocate(ADP_m_to_j(1:mesh*mesh*mesh))
      allocate(ADP_m_to_k(1:mesh*mesh*mesh))
      ADP_size_of_m =  size(ADP_m_to_i)

      do k = 1,mesh
        do j = 1,mesh
          do i = 1,mesh

              ADP_m_to_i(counter) = i
              ADP_m_to_j(counter) = j
              ADP_m_to_k(counter) = k
              counter = counter+1

          enddo
        enddo
      enddo

endif

    if (rank .eq. 0) then
      !write(*,*) 'ADP ',counter-1
      write(logf,*) "ADP evolution cells = ",ADP_size_of_m
      write(*,*) "ADP evolution cells = ",ADP_size_of_m
      flush(logf) 
    endif

      conv_criterion = min(int(convergence_fraction*mesh*mesh*mesh),(number_of_source-1)/3 )
      conv_criterion = max(conv_criterion,1)

  end subroutine ADP_array_initialization_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ADP_reset_photo_array()

    implicit none

    !integer :: i,j,k

    !do i=1,mesh
    !do j=1,mesh
    !do k=1,mesh

    !if (global_LTE_array(i,j,k) .lt. 0.0) then

    !  photoionization_HI_array(i,j,k) = 0.0
    !  photoionization_HeI_array(i,j,k) = 0.0 
    !  photoionization_HeII_array(i,j,k) = 0.0  
    !  heating_array(i,j,k) = 0.0

    !endif

    !enddo
    !enddo
    !enddo

      photoionization_HI_array = 0.0
      photoionization_HeI_array = 0.0 
      photoionization_HeII_array = 0.0  
      heating_array = 0.0
  end subroutine ADP_reset_photo_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ADP_array_destruction	

    implicit none

    deallocate(ADP_MPI_intermediate_average_xHI)
    deallocate(ADP_MPI_intermediate_average_xHII)
    deallocate(ADP_MPI_intermediate_average_xHeI)
    deallocate(ADP_MPI_intermediate_average_xHeII)
    deallocate(ADP_MPI_intermediate_average_xHeIII)
    deallocate(ADP_MPI_intermediate_average_temperature)
    deallocate(ADP_MPI_intermediate_end_xHI) 
    deallocate(ADP_MPI_intermediate_end_xHII)
    deallocate(ADP_MPI_intermediate_end_xHeI) 
    deallocate(ADP_MPI_intermediate_end_xHeII) 
    deallocate(ADP_MPI_intermediate_end_xHeIII) 
    deallocate(ADP_MPI_intermediate_end_temperature)
    deallocate(ADP_m_to_i)
    deallocate(ADP_m_to_j)
    deallocate(ADP_m_to_k)

  end subroutine ADP_array_destruction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ADP_array
