module sizes

  ! Module for C2-Ray
  ! Author: Garrelt Mellema
  ! Date: 2006-08-20
  ! This module is also accepted by the F compiler

  ! This module contains basic parameter definitions
  
  ! use precision
  implicit none

  private
  integer,parameter,public :: Ndim=1

  !Size of the mesh for spatial coordinate.
  integer,parameter,public :: mesh(1)=150

end module sizes


