module LTE_copy

  use parameter, only: mesh,LTE_size_of_m
  use array, only: xHI_array, xHII_array, xHeI_array, xHeII_array, xHeIII_array, temperature_array, &
                   global_LTE_array, intermediate_end_xHI, &
                   intermediate_end_xHII, intermediate_end_xHeI, intermediate_end_xHeII, &
                   intermediate_end_xHeIII,intermediate_end_temperature, &
                   global_LTE_array, global_evolved_LTE_array,LTE_m_to_i,LTE_m_to_j,LTE_m_to_k


  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LTE_copy_final_answer()

    implicit none

    integer :: i,j,k,m

    !do k = 1,mesh
    !  do j = 1,mesh
    !    do i = 1,mesh
    !      if (global_LTE_array(i,j,k).eq.1 .and. &
    !          global_evolved_LTE_array(i,j,k).eq.0) then

    !        xHI_array(i,j,k) = intermediate_end_xHI(i,j,k)
    !        xHII_array(i,j,k) = intermediate_end_xHII(i,j,k)
    !        xHeI_array(i,j,k) = intermediate_end_xHeI(i,j,k)
    !        xHeII_array(i,j,k) = intermediate_end_xHeII(i,j,k)
    !        xHeIII_array(i,j,k) = intermediate_end_xHeIII(i,j,k)
    !        temperature_array(i,j,k) = intermediate_end_temperature(i,j,k)

    !      endif
    !    enddo
    !  enddo
    !enddo

    do m = 1,LTE_size_of_m
      i = LTE_m_to_i(m)
      j = LTE_m_to_j(m)
      k = LTE_m_to_k(m)
      xHI_array(i,j,k) = intermediate_end_xHI(i,j,k)
      xHII_array(i,j,k) = intermediate_end_xHII(i,j,k)
      xHeI_array(i,j,k) = intermediate_end_xHeI(i,j,k)
      xHeII_array(i,j,k) = intermediate_end_xHeII(i,j,k)
      xHeIII_array(i,j,k) = intermediate_end_xHeIII(i,j,k)
      temperature_array(i,j,k) = intermediate_end_temperature(i,j,k)
    enddo

!write(*,*) 'after LTE evolution xHI_array(5,5,5)',xHI_array(5,5,5)
!write(*,*) 'after LTE evolution xHeI_array(5,5,5)',xHeI_array(5,5,5)
!write(*,*) 'after LTE evolution xHeII_array(5,5,5)',xHeII_array(5,5,5)
!write(*,*) 'after LTE evolution xHeIII_array(5,5,5)',xHeIII_array(5,5,5)

  end subroutine LTE_copy_final_answer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module LTE_copy
