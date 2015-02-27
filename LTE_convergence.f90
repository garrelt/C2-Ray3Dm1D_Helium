module LTE_convergence

  use c2ray_parameters, only: convergence_fraction
  use parameter, only: conv_criterion , LTE_convergence_failure, LTE_iteration_counter
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine LTE_convergence_test(converged,iterated)

    implicit none

    logical, intent(out) :: converged
    logical, intent(out) :: iterated

     ! Convergence test
     if (LTE_convergence_failure .le. conv_criterion) then
       converged = .true.
       !write(*,*) "Multiple sources convergence reached"
     else
       converged = .false.
     endif

     if (LTE_iteration_counter .gt. 100) then
       iterated = .true.
       write(*,*) 'Multiple sources not converging'
     else
       iterated = .false.
     endif

     ! Iteration loop counter
     LTE_iteration_counter = LTE_iteration_counter+1

  end subroutine LTE_convergence_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module LTE_convergence
