module density

  use my_mpi
  use file_admin, only: logf
  use precision, only: dp
  use parameter, only: mesh
  use cgsconstants, only: m_p
  use cosmology_parameters, only: Omega_B, rho_crit_0 
  use abundances, only: mu
  use array, only: number_density_array
  use parameter, only: mesh

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine density_initialization(zred_now)

    implicit none

    real(kind=dp),intent(in) :: zred_now
    integer :: i,j,k
    real(kind=dp) :: avg_dens
    integer :: mesh_x,mesh_y,mesh_z
    integer :: mympierror

    ! Assign density to the grid (average density at this redshift)
    !avg_dens = rho_crit_0*Omega_B/(mu*m_p)*(1.0+zred_now)**3
    !avg_dens = rho_crit_0*Omega_B/(mu*m_p)*(1.0+9.0)**3 ! Test 1
    !do k = 1,mesh
    !  do j = 1,mesh
    !    do i = 1,mesh
    !      number_density_array(i,j,k) = avg_dens
    !    enddo
    !  enddo
    !enddo

    !if (rank .eq. 0) then
    !  write(logf,*) "Cosmological density field = ",avg_dens," cm^-3" 
    !  flush(logf) 
    !endif

    if (rank .eq. 0) then

      open(unit=51,file="density_test4.bin",form="unformatted",status="old")
      read(51) mesh_x,mesh_y,mesh_z
      read(51) (((number_density_array(i,j,k),i=1,mesh_x),j=1,mesh_y),k=1,mesh_z)
      close(51)
      write(logf,*) "Reading density file from density_test4.bin" 
      flush(logf) 
    endif

    call MPI_BCAST(number_density_array,mesh*mesh*mesh,MPI_REAL,0,MPI_COMM_NEW,mympierror)

  end subroutine density_initialization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module density
