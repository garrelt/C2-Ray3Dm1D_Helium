module AT_ifront

  use precision, only: dp
  use parameter, only: mesh, do_periodic
  use array, only: global_Ifront_array, ionized_array, xHII_array, &
                   ifront_plus_ionized_array, non_screened_Ifront_array
  use parameter, only: time_partition
  use sourceprops, only: srcpos

  implicit none 
	
contains
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AT_Ifront_around_the_source(bubble,ns)
	
    implicit none

    logical,intent(out) :: bubble
    integer :: ns

    bubble = .true.
    if (xHII_array(srcpos(1,ns),srcpos(2,ns),srcpos(3,ns)).lt.real(time_partition)/real(time_partition+1)) then
    ! initial xHI = 1e-3
    !if (xHII_array(srcpos(1,ns),srcpos(2,ns),srcpos(3,ns)).lt.1.0-1.0e-3*real(time_partition)/real(time_partition+1)) then
      bubble = .false.
      non_screened_Ifront_array(0,0,0) = 1
    endif

  end subroutine AT_Ifront_around_the_source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! The i-front definition here can be different from the i-front definition in LTE_field.f90
  subroutine AT_Ifront_array_creation()

    implicit none
			
    integer :: i,j,k
    real(kind=dp) :: min_xHII
	
    !min_xHII = real(time_partition)/real(time_partition+1)
	
    ! initial xHI = 1e-3
    min_xHII = 1.0-1.0/real(time_partition+1)
   
    do i = 1, mesh
      do j = 1, mesh
        do k = 1, mesh

          if (xHII_array(i,j,k).ge.min_xHII) then
            ionized_array(i,j,k) = 1
            call AT_Ifront_plus_ionized_array_creation(i,j,k)
          endif

        enddo
      enddo
    enddo

    global_Ifront_array = ifront_plus_ionized_array - ionized_array

  end subroutine AT_Ifront_array_creation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine AT_Ifront_plus_ionized_array_creation(i,j,k)

    implicit none

    integer,intent(in) :: i, j, k
    integer :: i_minus_one, i_plus_one
    integer :: j_minus_one, j_plus_one
    integer :: k_minus_one, k_plus_one
		
    i_minus_one = i-1
    j_minus_one = j-1
    k_minus_one = k-1
    i_plus_one = i+1
    j_plus_one = j+1
    k_plus_one = k+1

    !check if boundary of the grid is met 
    if (i.eq.1) then
      if (do_periodic) then
        i_minus_one = mesh
      else
        i_minus_one = i
      endif
    endif

    if (j.eq.1) then
      if (do_periodic) then
        j_minus_one = mesh
      else
        j_minus_one = j
      endif
    endif

    if (k.eq.1) then
      if (do_periodic) then
        k_minus_one = mesh
      else
        k_minus_one = k
      endif
    endif

    if (i.eq.mesh) then
      if (do_periodic) then 
        i_plus_one = 1
      else
        i_plus_one = i
      endif
    endif

    if (j.eq.mesh) then
      if (do_periodic) then
        j_plus_one = 1
      else
        j_plus_one = j
      endif
    endif

    if (k.eq.mesh) then
      if (do_periodic) then
        k_plus_one = 1
      else
        k_plus_one = k
      endif
    endif

    ifront_plus_ionized_array(i_minus_one,j_minus_one,k_minus_one) = 1
    ifront_plus_ionized_array(i_minus_one,j_minus_one,k) = 1
    ifront_plus_ionized_array(i_minus_one,j_minus_one,k_plus_one) = 1
    ifront_plus_ionized_array(i_minus_one,j,k_minus_one) = 1
    ifront_plus_ionized_array(i_minus_one,j,k) = 1
    ifront_plus_ionized_array(i_minus_one,j,k_plus_one) = 1
    ifront_plus_ionized_array(i_minus_one,j_plus_one,k_minus_one) = 1
    ifront_plus_ionized_array(i_minus_one,j_plus_one,k) = 1
    ifront_plus_ionized_array(i_minus_one,j_plus_one,k_plus_one) = 1
		
    ifront_plus_ionized_array(i,j_minus_one,k_minus_one) = 1
    ifront_plus_ionized_array(i,j_minus_one,k) = 1
    ifront_plus_ionized_array(i,j_minus_one,k_plus_one) = 1
    ifront_plus_ionized_array(i,j,k_minus_one) = 1
    ifront_plus_ionized_array(i,j,k) = 1
    ifront_plus_ionized_array(i,j,k_plus_one) = 1
    ifront_plus_ionized_array(i,j_plus_one,k_minus_one) = 1
    ifront_plus_ionized_array(i,j_plus_one,k) = 1
    ifront_plus_ionized_array(i,j_plus_one,k_plus_one) = 1
		
    ifront_plus_ionized_array(i_plus_one,j_minus_one,k_minus_one) = 1
    ifront_plus_ionized_array(i_plus_one,j_minus_one,k) = 1
    ifront_plus_ionized_array(i_plus_one,j_minus_one,k_plus_one) = 1
    ifront_plus_ionized_array(i_plus_one,j,k_minus_one) = 1
    ifront_plus_ionized_array(i_plus_one,j,k) = 1
    ifront_plus_ionized_array(i_plus_one,j,k_plus_one) = 1
    ifront_plus_ionized_array(i_plus_one,j_plus_one,k_minus_one) = 1
    ifront_plus_ionized_array(i_plus_one,j_plus_one,k) = 1
    ifront_plus_ionized_array(i_plus_one,j_plus_one,k_plus_one) = 1

  end subroutine AT_Ifront_plus_ionized_array_creation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module AT_ifront
