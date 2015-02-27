module mathconstant

  use precision, only: dp

  implicit none

  ! the number pi
  real(kind=dp), parameter :: pi = 3.141592654
  real(kind=dp), parameter :: four_pi = 4*3.141592654 
  real(kind=dp), parameter :: four_over_three_pi = 4.0*pi/3.0

end module mathconstant
