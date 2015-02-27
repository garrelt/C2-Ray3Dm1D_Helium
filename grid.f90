module grid

  use precision, only: dp
  use astroconstants, only: Mpc
  use cosmology_parameters, only: h
  use my_mpi
  use file_admin, only: logf
  use parameter, only : cell_size, cell_volume, grid_volume, boxsize, mesh

  implicit none

  ! spatial coordinate x
  real(kind=dp),dimension(:),allocatable :: x_coordinate 
  ! spatial coordinate y
  real(kind=dp),dimension(:),allocatable :: y_coordinate 
  ! spatial coordinate z
  real(kind=dp),dimension(:),allocatable :: z_coordinate 
  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initializes grid properties
  subroutine grid_initialization()

    implicit none

    integer :: i,j,k
    integer :: alloc_status
    real(kind=dp) :: xgrid,ygrid,zgrid
    integer :: ierror

    if (rank .eq. 0) then
      write(logf,*) "Beginning of grid initialization" 
    endif

    ! Simulation volume (comoving)
    grid_volume = (boxsize/h*Mpc)**3
    
    ! Ask for grid size (if rank 0 and not set in nbody module)
    if (rank .eq. 0) then
      xgrid = boxsize
      ygrid = boxsize
      zgrid = boxsize
      ! Report
      write(logf,*) "Box size = ", boxsize, " Mpc/h (comoving)"
      write(logf,*) "Resolution = ", mesh,"^3"

    endif
    
#ifdef MPI
    call MPI_BCAST(xgrid,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(ygrid,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
    call MPI_BCAST(zgrid,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,ierror)
#endif
    
    xgrid = xgrid*Mpc/h
    ygrid = ygrid*Mpc/h
    zgrid = zgrid*Mpc/h
    
    cell_size = xgrid/real(mesh)

    ! allocate coordinates
    allocate(x_coordinate(mesh),stat=alloc_status)
    allocate(y_coordinate(mesh),stat=alloc_status)
    allocate(z_coordinate(mesh),stat=alloc_status)

    ! coordinates of a cell
    do i=1,mesh
      x_coordinate(i) = (real(i)-0.5)*cell_size
    enddo

    do j=1,mesh
      y_coordinate(j) = (real(j)-0.5)*cell_size
    enddo

    do k=1,mesh
      z_coordinate(k) = (real(k)-0.5)*cell_size
    enddo
    
    cell_volume = cell_size*cell_size*cell_size

    if (rank .eq. 0) then
      write(logf,*) "End of grid initialization" 
      write(logf,*)      
      flush(logf) 
    endif

  end subroutine grid_initialization
  
end module grid
