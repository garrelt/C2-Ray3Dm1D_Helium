module ADP_copy

  use my_mpi
  use file_admin, only: logf
  use parameter, only: mesh
  use array, only: temperature_array
  use array, only: xHI_array,xHII_array
  use array, only: xHeI_array,xHeII_array,xHeIII_array
  use array, only: intermediate_end_xHI,intermediate_end_xHII
  use array, only: intermediate_end_xHeI,intermediate_end_xHeII,intermediate_end_xHeIII
  use array, only: intermediate_end_temperature
  use array, only: global_LTE_array

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ADP_copy_final_answer()

    implicit none

    integer :: i,j,k

    do k = 1,mesh
      do j = 1,mesh
        do i = 1,mesh
          if (global_LTE_array(i,j,k).lt.0) then
 
!if (isnan(intermediate_end_xHI(i,j,k))) write(*,*) 'intermediate_end_xHI',i,j,k
!if (isnan(intermediate_end_xHII(i,j,k))) write(*,*) 'intermediate_end_xHII',i,j,k
!if (isnan(intermediate_end_xHeI(i,j,k))) write(*,*) 'intermediate_end_xHeI',i,j,k
!if (isnan(intermediate_end_xHeII(i,j,k))) write(*,*) 'intermediate_end_xHeII',i,j,k
!if (isnan(intermediate_end_xHeIII(i,j,k))) write(*,*) 'intermediate_end_xHeIII',i,j,k
!if (isnan(intermediate_end_temperature(i,j,k))) write(*,*) 'intermediate_end_temperature',i,j,k
          xHI_array(i,j,k)=intermediate_end_xHI(i,j,k)
          xHII_array(i,j,k)=intermediate_end_xHII(i,j,k)
          xHeI_array(i,j,k)=intermediate_end_xHeI(i,j,k)
          xHeII_array(i,j,k)=intermediate_end_xHeII(i,j,k)
          xHeIII_array(i,j,k)=intermediate_end_xHeIII(i,j,k)
          temperature_array(i,j,k) = intermediate_end_temperature(i,j,k)

          endif
        enddo
      enddo
    enddo

!write(*,*) 'after ADP evolution xHI_array(5,5,5)',xHI_array(5,5,5)
!write(*,*) 'after ADP evolution xHeI_array(5,5,5)',xHeI_array(5,5,5)
!write(*,*) 'after ADP evolution xHeII_array(5,5,5)',xHeII_array(5,5,5)
!write(*,*) 'after ADP evolution xHeIII_array(5,5,5)',xHeIII_array(5,5,5)

    if (rank .eq. 0) then
      write(logf,*) "Finsih ADP evolution"
      flush(logf) 
    endif

  end subroutine ADP_copy_final_answer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ADP_copy
