!!This module contains data and routines for MPI parallelization

module my_mpi

  use file_admin, only: logf, results_dir

#ifdef MY_OPENMP
  USE OMP_LIB, only: omp_get_num_threads, omp_get_thread_num
#endif

#ifdef IFORT
  USE IFPORT, only: hostnm
#endif

  implicit none

  include 'mpif.h'

  ! rank of the processor
  integer :: rank = 0  
          
  ! number of the processors
  integer :: nprocs = 1   
          
  ! number of threads per processor
  integer,public :: nthreads = 1  
    
  ! the MPI communicator
  integer,public :: MPI_COMM_NEW 
  
  ! MPI status array
  integer,public,dimension(MPI_STATUS_SIZE) :: mympi_status

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_setup()

    ! name of the log file
    character(len=512) :: filename  
      
    ! control variable for MPI
    integer :: ierror

    ! thread 
    integer :: thread_id

    character(len=100) :: hostname

    ! Initialize MPI
    call MPI_INIT(ierror) 
 
    ! Find processor rank
    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror) 

    ! Find total number of processors
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierror)

    if (rank.eq.0) write(*,*) 'number of MPI processors = ',nprocs
    write(*,*) 'my rank = ', rank

    ! Only rank 0 processor can recognize the C2Ray.log
    if (rank.eq.0) then
      filename = trim(adjustl(trim(adjustl(results_dir))//"C2Ray.log"))
      open(unit=logf,file=filename,status="unknown",action="write")

      write(logf,"(A)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(logf,"(A)") "!!                                                                 !!"
      write(logf,"(A)") "!!                   Log file for C2-Ray run                       !!" 
      write(logf,"(A)") "!!                                                                 !!"
      write(logf,"(A)") "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      write(logf,*) 
      write(logf,*) "Beginning of mpi setup"
      write(logf,*) "Number of MPI processors used: ", nprocs
    endif



    ! Find number of OpenMP threads 
    !$omp parallel default(shared)
#ifdef MY_OPENMP
    nthreads = omp_get_num_threads()
#endif
    !$omp end parallel

    ! Report OpenMP usage
#ifdef MY_OPENMP
    write(logf,*) "Running in OpenMP mode"
    write(logf,*) 'Number of OpenMP threads per processor: ',nthreads
    if (rank .eq. 0) write(*,*) 'Number of OpenMP threads per processor: ',nthreads
    flush(logf)
#endif

    ! Figure out hostname
    ! NOTE: compiler dependent!!!
    ierror = hostnm(hostname)
    if (ierror .eq. 0) then
      write(logf,*) "Running on processor named ",trim(adjustl(hostname))
    else 
      write(logf,*) "Error establishing identity of processor."
    endif

    ! Let OpenMP threads report

    !$omp parallel default(private)
    thread_id = omp_get_thread_num()+1
    write(logf,*) 'Thread number ',thread_id,' of processor ',rank,' reporting'
    !$omp end parallel

    MPI_COMM_NEW = MPI_COMM_WORLD

#ifdef MPI
    write(logf,*) "End of mpi setup"
    write(logf,*)
    flush(logf)
#endif

  end subroutine mpi_setup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_end

    integer :: ierror
    logical :: openlog

    ! Find out if log file is open
    inquire(unit=logf,opened=openlog)

    ! Close log file
    if (openlog) close(logf)

    ! Close MPI
    call MPI_FINALIZE(ierror)

  end subroutine mpi_end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module my_mpi
