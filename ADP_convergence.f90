module ADP_convergence

  use my_mpi
  use file_admin, only: logf
  use c2ray_parameters, only: convergence_fraction
  use parameter, only: conv_criterion , convergence_failure, iteration_counter
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ADP_convergence_test(converged,iterated)

    implicit none

    logical, intent(out) :: converged
    logical, intent(out) :: iterated

     ! Convergence test
     if (convergence_failure .lt. conv_criterion) then
       converged = .true.
     else
       converged = .false.
     endif

     if (iteration_counter .gt. 100) then
       iterated = .true.
     else
       iterated = .false.
     endif

     ! Iteration loop counter
     iteration_counter = iteration_counter+1

    if (rank .eq. 0) then
      if (converged) write(logf,*) "converged with iteration counter ", iteration_counter
      if (iterated) write(logf,*) "iterated with iteration counter ", iteration_counter
      flush(logf) 
    endif

  end subroutine ADP_convergence_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ADP_convergence
