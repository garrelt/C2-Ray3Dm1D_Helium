!>
!! \brief This module contains routines for calculating the ionization and temperature evolution of the entire grid (3D).
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date:
!!
!! \b Version: 3D, MPI & OpenMP

module evolve

  ! This module contains routines having to do with the calculation of
  ! the ionization evolution of the entire grid (3D).
     
  ! This version has been adapted for efficiency in order to be able
  ! to calculate large meshes.
    
  ! - evolve3D : step through grid
  ! - evolve2D : step through z-plane
  ! - evolve0D : take care of one grid point
  ! - evolve0D_global: update entire grid

  ! Needs:
  ! doric : ionization calculation for one point + photo-ionization rates
  ! tped : temperature,pressure,electron density calculation

  use precision, only: dp
  use my_mpi ! supplies all the MPI and OpenMP definitions and variables
  use file_admin, only: logf,iterdump
  use sizes, only: Ndim, mesh
  use grid, only: x,y,z,vol,dr
  use material, only: ndens, xh,xhe, get_temperature_point, &
      set_temperature_point,set_final_temperature_point,isothermal
  use sourceprops, only: SrcSeries, NumSrc, srcpos, NormFlux, NormFluxPL
  use radiation, only: NumFreqBnd
  use photonstatistics, only: state_before, calculate_photon_statistics, &
       photon_loss, LLS_loss, report_photonstatistics, state_after, total_rates, &
       total_ionizations, update_grandtotal_photonstatistics
  use c2ray_parameters, only: convergence_fraction, subboxsize, max_subbox, S_star_nominal
  use abundances, only: abu_he
  use cgsconstants, only: ini_rec_colion_factors
  use thermalevolution, only: thermal
  
  implicit none

  private

  public :: evolve3D, phih_grid, evolve_ini, phihe_grid

  !> Periodic boundary conditions, has to be true for this version
  logical,parameter :: periodic_bc = .true.

  !> Minimum number of MPI processes for using the master-slave setup 
  integer, parameter ::  min_numproc_master_slave=10 !10

  ! Grid variables

  !> H Photo-ionization rate on the entire grid
  real(kind=dp),dimension(:,:,:),allocatable :: phih_grid
  !> He Photo-ionization rate on the entire grid
  real(kind=dp),dimension(:,:,:,:),allocatable :: phihe_grid
  !> Heating  rate on the entire grid
  real(kind=dp),dimension(:,:,:),allocatable :: phiheat
  
  !> Time-averaged H ionization fraction
  real(kind=dp),dimension(:,:,:,:),allocatable :: xh_av
  !> Time-averaged He ionization fraction
  real(kind=dp),dimension(:,:,:,:),allocatable :: xhe_av
  !> Intermediate result for H ionization fraction
  real(kind=dp),dimension(:,:,:,:),allocatable :: xh_intermed
  !> Intermediate result for He ionization fraction
  real(kind=dp),dimension(:,:,:,:),allocatable :: xhe_intermed
  !> H0 Column density (outgoing)
  real(kind=dp),dimension(:,:,:),allocatable :: coldensh_out
  !> He Column density (outgoing)
  real(kind=dp),dimension(:,:,:,:),allocatable :: coldenshe_out
  !> Buffer for MPI communication
  real(kind=dp),dimension(:,:,:),allocatable :: buffer
  !> Photon loss from the grid
  real(kind=dp) :: photon_loss_all(1:NumFreqBnd)
  !> Photon loss from one source
  real(kind=dp),dimension(:,:),allocatable :: photon_loss_src_thread
  real(kind=dp) :: photon_loss_src(1:NumFreqBnd)

  ! mesh positions of end points for RT
  integer,dimension(Ndim) :: lastpos_l !< mesh position of left end point for RT
  integer,dimension(Ndim) :: lastpos_r !< mesh position of right end point for RT
  integer,dimension(Ndim) :: last_l !< mesh position of left end point for RT
  integer,dimension(Ndim) :: last_r !< mesh position of right end point for RT

  integer :: sum_nbox !< sum of all nboxes (on one processor)
  integer :: sum_nbox_all !< sum of all nboxes (on all processors)

  ! GM/121127: This variable should always be set. If not running OpenMP
  ! it should be equal to 1. We initialize it to 1 here.
  integer :: tn=1 !< thread number
  
contains

  ! =======================================================================

  !> Allocate the arrays needed for evolve
  subroutine evolve_ini ()
    
    allocate(phih_grid(mesh(1),mesh(2),mesh(3)))
    allocate(phihe_grid(mesh(1),mesh(2),mesh(3),0:1))
    allocate(phiheat(mesh(1),mesh(2),mesh(3)))
    
    allocate(xh_av(mesh(1),mesh(2),mesh(3),0:1))
    allocate(xhe_av(mesh(1),mesh(2),mesh(3),0:2))

    allocate(xh_intermed(mesh(1),mesh(2),mesh(3),0:1))
    allocate(xhe_intermed(mesh(1),mesh(2),mesh(3),0:2))

    allocate(coldensh_out(mesh(1),mesh(2),mesh(3)))
    allocate(coldenshe_out(mesh(1),mesh(2),mesh(3),0:1))

    allocate(buffer(mesh(1),mesh(2),mesh(3)))
    allocate(photon_loss_src_thread(1:NumFreqBnd,nthreads))
  !   allocate(photon_loss_src_thread(1:NumFreqBnd,1))   
  end subroutine evolve_ini

  ! ===========================================================================

  !> Evolve the entire grid over a time step dt
  subroutine evolve3D (dt,restart)

    ! Calculates the evolution of the hydrogen ionization state
     
    ! Author: Garrelt Mellema
     
    ! Date: 28-Feb-2008 (21-Aug-2006 (f77/OMP: 13-Jun-2005))

    ! Version: Multiple sources / Using average fractions to converge
    ! loop over sources
    
    ! History:
    ! 11-Jun-2004 (GM) : grid arrays now passed via common (in grid.h)
    !    and material arrays also (in material.h).
    ! 11-Jun-2004 (GM) : adapted for multiple sources.
    !  3-Jan-2005 (GM) : reintegrated with updated Ifront3D
    ! 20-May-2005 (GM) : split original eveolve0D into two routines
    ! 13-Jun-2005 (HM) : OpenMP version : Hugh Merz
    ! 21-Aug-2006 (GM) : MPI parallelization over the sources (static model).
    ! 28-Feb-2008 (GM) : Added master-slave model for distributing
    !                    over the processes. The program decides which
    !                    model to use.

    ! The time step
    real(kind=dp),intent(in) :: dt !< time step
    integer,intent(in) :: restart !< restart flag

    ! Loop variables
    integer :: niter  ! iteration counter

    ! Wall clock counting
    integer(kind=8) :: wallclock1
    integer(kind=8) :: wallclock2
    integer(kind=8) :: countspersec
    ! Flag variable (passed back from evolve0D_global)
    integer :: conv_flag
    ! Minimum number of cells which are allowed to be non-converged
    integer :: conv_criterion 

#ifdef MPI
    integer :: mympierror
#endif

    ! End of declarations

    ! Initialize wall clock counter (for dumps)
    call system_clock(wallclock1)

     ! Initial state (for photon statistics)
    call state_before (xh,xhe)

    ! initialize average and intermediate results to initial values
    if (restart == 0) then
       xh_av(:,:,:,:)=xh(:,:,:,:)
       xh_intermed(:,:,:,:)=xh(:,:,:,:)
       xhe_av(:,:,:,:)=xhe(:,:,:,:)
       xhe_intermed(:,:,:,:)=xhe(:,:,:,:)
       niter=0 ! iteration starts at zeroF
       conv_flag=mesh(1)*mesh(2)*mesh(3) ! initialize non-convergence 
    else
       ! Reload xh_av,xh_intermed,photon_loss,niter
       call start_from_dump(niter)
       call global_pass (conv_flag,dt)
    endif
    
    ! Set the conv_criterion, if there are few sources we should make
    ! sure that things are converged around these sources.

   ! conv_criterion=min(int(convergence_fraction*mesh(1)*mesh(2)*mesh(3)),(NumSrc-1)/3 )
    conv_criterion=min(int(convergence_fraction*mesh(1)*mesh(2)*mesh(3)),NumSrc)   
    
    ! Iterate to reach convergence for multiple sources
#ifdef MPILOG
    ! Report
    write(logf,*) 'About to start evolve '
    flush(logf)
#endif
    do
       ! Update xh if converged and exit
!       write(*,*) 'start of do'
       if (conv_flag < conv_criterion) then
          xh(:,:,:,:)=xh_intermed(:,:,:,:)
          xhe(:,:,:,:)=xhe_intermed(:,:,:,:)
          call set_final_temperature_point
          exit
       else
          if (niter > 500) then
             ! Complain about slow convergence
             if (rank == 0) write(logf,*) 'Multiple sources not converging'
             exit
          endif
       endif
 
!*** TEST
      conv_flag=0.0
 
       ! Iteration loop counter
       niter=niter+1
!       write(*,*) 'before pass_all_sources'
 !      write(*,*) niter,dt
       call pass_all_sources (niter,dt)
 !   write(*,*) 'after pass_all_sources'
        
       ! Report subbox statistics
       if (rank == 0) &
            write(logf,*) "Average number of subboxes: ",sum_nbox_all/NumSrc

       if (rank == 0) then
          call system_clock(wallclock2,countspersec)
          ! Write iteration dump if more than 15 minutes have passed.
          ! system_clock starts counting at 0 when it reaches
          ! a max value. To catch this, test also for negative
          ! values of wallclock2-wallclock1
          if (wallclock2-wallclock1 > 15*60*countspersec .or. &
               wallclock2-wallclock1 < 0 ) then
             call write_iteration_dump(niter)
             wallclock1=wallclock2
          endif
       endif
 !      write(*,*) 'before global pass'
       call global_pass (conv_flag,dt)  
 !      write(*,*) 'after global pass'
    enddo 
    ! Calculate photon statistics
    call calculate_photon_statistics (dt,xh,xh_av,xhe,xhe_av) 
   !      write(*,*) 'after calc photon stats'  
    call report_photonstatistics (dt)
  !        write(*,*) 'after report photon stats'     
    call update_grandtotal_photonstatistics (dt)
  !           write(*,*) 'end evolve 3D'    
  end subroutine evolve3D

  ! ===========================================================================
  !!!*** add T
  subroutine write_iteration_dump (niter)   !*** This interdumfile has to include he, for now, just added in the end
    use material, only:temperature_grid
    integer,intent(in) :: niter  ! iteration counter

    integer :: ndump=0
    
    character(len=20) :: iterfile

    ndump=ndump+1
    if (mod(ndump,2) == 0) then
       iterfile="iterdump2.bin"
    else
       iterfile="iterdump1.bin"
    endif

    open(unit=iterdump,file=trim(adjustl(results_dir))//iterfile,form="unformatted", &
         status="unknown")

    write(iterdump) niter
    write(iterdump) photon_loss_all
    write(iterdump) phih_grid
    write(iterdump) xh_av
    write(iterdump) xh_intermed
    write(iterdump) phihe_grid
    write(iterdump) xhe_av
    write(iterdump) xhe_intermed
    if (.not.isothermal) then
       write(iterdump) phiheat
       write(iterdump) temperature_grid
    endif 
    close(iterdump)

  end subroutine write_iteration_dump

  ! ===========================================================================
  !!!*** add T
  subroutine start_from_dump(niter)  !***read the same way as written
    
    use material, only:temperature_grid
    integer,intent(out) :: niter  ! iteration counter

 
#ifdef MPI
    integer :: mympierror
#endif
   
    if (rank == 0) then
       open(unit=iterdump,file="iterdump.bin",form="unformatted", &
            status="old")
       
       read(iterdump) niter
       read(iterdump) photon_loss_all
       read(iterdump) phih_grid
       read(iterdump) xh_av
       read(iterdump) xh_intermed
       read(iterdump) phihe_grid
       read(iterdump) xhe_av
       read(iterdump) xhe_intermed
	 if (.not.isothermal) then
	     read(iterdump) phiheat
	     read(iterdump) temperature_grid
	 endif 
       
       close(iterdump)
       write(logf,*) "Read iteration ",niter," from dump file"
       write(logf,*) 'photon loss counter: ',photon_loss_all
       write(logf,*) "Intermediate result for mean ionization fraction: ", &
            sum(xh_intermed(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3)), &
            sum(xhe_intermed(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3)), &
            sum(xhe_intermed(:,:,:,2))/real(mesh(1)*mesh(2)*mesh(3))
    endif

#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(niter,1, &
         MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(photon_loss_all,NumFreqBnd, &
         MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(phih_grid,mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(phihe_grid,mesh(1)*mesh(2)*mesh(3)*2, &
         MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(phiheat,mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)        
    call MPI_BCAST(xh_av,mesh(1)*mesh(2)*mesh(3)*2, &
         MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(xhe_av,mesh(1)*mesh(2)*mesh(3)*3, &
         MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
    call MPI_BCAST(xh_intermed,mesh(1)*mesh(2)*mesh(3)*2, &
         MPI_DOUBLE_PRECISION,0,&
         MPI_COMM_NEW,mympierror)
    call MPI_BCAST(xhe_intermed,mesh(1)*mesh(2)*mesh(3)*3, &
         MPI_DOUBLE_PRECISION,0,&
         MPI_COMM_NEW,mympierror)
#endif
    
  end subroutine start_from_dump

  ! ===========================================================================

  subroutine pass_all_sources(niter,dt)
    
    ! For random permutation of sources:
    use  m_ctrper, only: ctrper
    use c2ray_parameters, only: use_LLS

    integer,intent(in) :: niter  ! iteration counter
    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    
    real(kind=dp) :: LLS_loss_all
    
    integer :: i,j,k
#ifdef MPI
    integer :: mympierror
#endif

 !   write(*,*) 'in pass all sources'
    if (rank == 0) write(logf,*) 'Doing all sources '
    ! reset global rates to zero for this iteration
    phih_grid(:,:,:)=0.0
    phihe_grid(:,:,:,:)=0.0    
    phiheat(:,:,:)=0.0
    ! reset photon loss counter
    photon_loss(:)=0.0
    LLS_loss = 0.0 ! make this a NumFreqBnd vector if needed later (GM/101129)
    
    ! Reset sum of subboxes counter
    sum_nbox=0

    ! Make a randomized list of sources :: call in serial
    if ( rank == 0 ) call ctrper (SrcSeries(1:NumSrc),1.0)
 
#ifdef MPI
    ! Distribute the source list to the other nodes
    call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
    
#ifdef MPILOG
    ! Report
    write(logf,*) 'after broadcasting sources'
    flush(logf)
#endif
    
  !  write(*,*) 'npr', npr
    ! Ray trace the whole grid for all sources.
    ! We can do this in two ways, depending on
    ! the number of processors. For many processors
    ! the master-slave setup should be more efficient.
    if (npr > min_numproc_master_slave) then
#ifdef MPILOG
    ! Report
    write(logf,*) 'Doing master-slave parallelization '
    flush(logf)
#endif
       call do_grid_master_slave (dt,niter)

    else
#ifdef MPILOG
    ! Report
    write(logf,*) 'Doing static parallelization '
    flush(logf)
#endif

    !   write(*,*) 'before do_grid_static'
       call do_grid_static (dt,niter)
    !   write(*,*) 'after do_grid_static'
    endif

#ifdef MPI
#ifdef MPILOG
    write(logf,*) "Start accumulation"
#endif
    ! accumulate (sum) the MPI distributed photon losses
    call MPI_ALLREDUCE(photon_loss, photon_loss_all, NumFreqBnd, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
         
    ! accumulate (sum) the MPI distributed photon losses
    if (use_LLS)then
    	call MPI_ALLREDUCE(LLS_loss, LLS_loss_all, 1, &
           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
    	! Put LLS_loss_all back in the LLS variable
    	LLS_loss = LLS_loss_all
    endif
#ifdef MPILOG
    write(logf,*) "Photon loss"
    flush(logf)
#endif
    ! accumulate (sum) the MPI distributed phih_grid
    call MPI_ALLREDUCE(phih_grid, buffer, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
    ! Overwrite the processor local values with the accumulated value
    phih_grid(:,:,:)=buffer(:,:,:)

#ifdef MPILOG
    write(logf,*) "Accumulation 2"
    flush(logf)
#endif
     call MPI_ALLREDUCE(phihe_grid(:,:,:,0), buffer, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)    
    ! Overwrite the processor local values with the accumulated value
    phihe_grid(:,:,:,0)=buffer(:,:,:)

#ifdef MPILOG
    write(logf,*) "Accumulation 3"
    flush(logf)
#endif
    call MPI_ALLREDUCE(phihe_grid(:,:,:,1), buffer, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)    
    ! Overwrite the processor local values with the accumulated value
    phihe_grid(:,:,:,1)=buffer(:,:,:)
    
    
    call MPI_ALLREDUCE(phiheat, buffer, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
    ! Overwrite the processor local values with the accumulated value
    phiheat(:,:,:)=buffer(:,:,:)    
        

#ifdef MPILOG
    write(logf,*) "Accumulation 4"
    flush(logf)
#endif
    ! accumulate (sum) the MPI distributed sum of number of boxes
    call MPI_ALLREDUCE(sum_nbox, sum_nbox_all, 1, &
         MPI_INTEGER, MPI_SUM, MPI_COMM_NEW, mympierror)
#ifdef MPILOG
    write(logf,*) "Accumulation 5"
    flush(logf)
#endif

    ! Only on the first iteration does evolve2D (evolve0D) change the
    ! ionization fractions
!***
    if (niter == 1) then
       ! accumulate (max) MPI distributed xh_av
       call MPI_ALLREDUCE(xh_av(:,:,:,1), buffer, mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_NEW, mympierror)
       ! Overwrite the processor local values with the accumulated value
       xh_av(:,:,:,1) = buffer(:,:,:)

#ifdef MPILOG
    write(logf,*) "Accumulation 6"
    flush(logf)
#endif
        call MPI_ALLREDUCE(xhe_av(:,:,:,0), buffer, mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_NEW, mympierror)      
       ! Overwrite the processor local values with the accumulated value
       xhe_av(:,:,:,0) = buffer(:,:,:)

#ifdef MPILOG
    write(logf,*) "Accumulation 7"
    flush(logf)
#endif
        call MPI_ALLREDUCE(xhe_av(:,:,:,2), buffer, mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_NEW, mympierror)   
       ! Overwrite the processor local values with the accumulated value
       xhe_av(:,:,:,2) = buffer(:,:,:)
#ifdef MPILOG
    write(logf,*) "Accumulation 8"
    flush(logf)
#endif
       !xh_av(:,:,:,0) = max(0.0_dp,min(1.0_dp,1.0-xh_av(:,:,:,1)))
       !xhe_av(:,:,:,1) = max(0.0_dp,min(1.0_dp,1.0-xhe_av(:,:,:,0)-xhe_av(:,:,:,2)))    
       xh_av(:,:,:,0) = 1.0-xh_av(:,:,:,1)
       xhe_av(:,:,:,1) = 1.0-xhe_av(:,:,:,0)-xhe_av(:,:,:,2)
       do k=1,mesh(3)
          do j=1,mesh(2)
             do i=1,mesh(1)
                if (xh_av(i,j,k,0) < 0.0d0) then
#ifdef MPILOG
                   write(logf,*) 'xh_av(0) < 0 at ',i,j,k,xh_av(i,j,k,0)
#endif
                   xh_av(i,j,k,0)=0.0d0
                endif
                if (xhe_av(i,j,k,0) < 0.0d0) then
#ifdef MPILOG
                   write(logf,*) 'xhe_av(0) < 0 at ',i,j,k,xhe_av(i,j,k,0)
#endif
                   xhe_av(i,j,k,0)=0.0d0
                endif
                if (xh_av(i,j,k,0) > 1.0d0) then
#ifdef MPILOG
                   write(logf,*) 'xh_av(0) > 1 at ',i,j,k,xh_av(i,j,k,0)
#endif
                   xh_av(i,j,k,0)=1.0d0
                endif
                if (xhe_av(i,j,k,0) > 1.0d0) then
#ifdef MPILOG
                   write(logf,*) 'xhe_av(0) > 1 at ',i,j,k,xhe_av(i,j,k,0)
#endif
                   xhe_av(i,j,k,0)=1.0d0
                endif
             enddo
          enddo
       enddo
  
       ! accumulate (max) MPI distributed xh_intermed
       call MPI_ALLREDUCE(xh_intermed(:,:,:,1), buffer, &
            mesh(1)*mesh(2)*mesh(3), MPI_DOUBLE_PRECISION, MPI_MAX, &
            MPI_COMM_NEW, mympierror)
       ! Overwrite the processor local values with the accumulated value
       xh_intermed(:,:,:,1)=buffer(:,:,:)

#ifdef MPILOG
    write(logf,*) "Accumulation 9"
    flush(logf)
#endif
       call MPI_ALLREDUCE(xhe_intermed(:,:,:,0), buffer, &
            mesh(1)*mesh(2)*mesh(3), MPI_DOUBLE_PRECISION, MPI_MIN, &
            MPI_COMM_NEW, mympierror)
       ! Overwrite the processor local values with the accumulated value
       xhe_intermed(:,:,:,0)=buffer(:,:,:)

#ifdef MPILOG
    write(logf,*) "Accumulation 10"
    flush(logf)
#endif
       call MPI_ALLREDUCE(xhe_intermed(:,:,:,2), buffer, &
            mesh(1)*mesh(2)*mesh(3), MPI_DOUBLE_PRECISION, MPI_MAX, &
            MPI_COMM_NEW, mympierror)
       ! Overwrite the processor local values with the accumulated value
       xhe_intermed(:,:,:,2)=buffer(:,:,:)

#ifdef MPILOG
    write(logf,*) "Accumulation 11"
    flush(logf)
#endif
       !xh_intermed(:,:,:,0)=max(0.0_dp,min(1.0_dp,1.0-xh_intermed(:,:,:,1)))
       !xhe_intermed(:,:,:,1)=max(0.0_dp,min(1.0_dp,1.0-xhe_intermed(:,:,:,0)-xhe_intermed(:,:,:,2)))
       xh_intermed(:,:,:,0)=1.0-xh_intermed(:,:,:,1)
       xhe_intermed(:,:,:,1)=1.0-xhe_intermed(:,:,:,0)-xhe_intermed(:,:,:,2)
       do k=1,mesh(3)
          do j=1,mesh(2)
             do i=1,mesh(1)
                if (xh_intermed(i,j,k,0) < 0.0d0) then
#ifdef MPILOG
                   write(logf,*) 'xh_intermed(0) < 0 at ',i,j,k,xh_intermed(i,j,k,0)
#endif
                   xh_intermed(i,j,k,0)=0.0d0
                endif
                if (xhe_intermed(i,j,k,0) < 0.0d0) then
#ifdef MPILOG
                   write(logf,*) 'xhe_intermed(0) < 0 at ',i,j,k,xhe_intermed(i,j,k,0)
#endif
                   xhe_intermed(i,j,k,0)=0.0d0
                endif
                if (xh_intermed(i,j,k,0) > 1.0d0) then
#ifdef MPILOG
                   write(logf,*) 'xh_intermed(0) > 1 at ',i,j,k,xh_intermed(i,j,k,0)
#endif
                   xh_intermed(i,j,k,0)=1.0d0
                endif
                if (xhe_intermed(i,j,k,0) > 1.0d0) then
#ifdef MPILOG
                   write(logf,*) 'xhe_intermed(0) > 1 at ',i,j,k,xhe_intermed(i,j,k,0)
#endif
                   xhe_intermed(i,j,k,0)=1.0d0
                endif
             enddo
          enddo
       enddo
    endif
#else
    photon_loss_all(:)=photon_loss(:)
    sum_nbox_all=sum_nbox
#endif
    
#ifdef MPI
    call MPI_BARRIER(MPI_COMM_NEW,mympierror)
#endif

#ifdef MPILOG
    write(logf,*) "End accumulation"
#endif
    
  end subroutine pass_all_sources

  ! ===========================================================================

  subroutine global_pass (conv_flag,dt)

    ! Flag variable (passed back from evolve0D_global)
    integer,intent(out) :: conv_flag
    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D

    integer :: i,j,k  ! mesh position

    ! Mesh position of the cell being treated
    integer,dimension(Ndim) :: pos

#ifdef MPI
    integer :: mympierror
#endif

    ! Report photon losses over grid boundary 
   ! if (rank == 0) write(logf,*) 'photon loss counter: ',photon_loss_all(:) 
    
    ! Turn total photon loss into a mean per cell (used in evolve0d_global)
    ! GM/110225: Possibly this should be 
    ! photon_loss_all(:)/(real(mesh(1))*real(mesh(2))
    ! Since this is the correct answer for an optically thin volume
    ! Compromise: photon_loss_all(:)/(real(mesh(1))**3)*sum(xh_av(:,:,:,1))
    ! then if the box is fully ionized the optically thin 1/N^2 is used
    ! and if the box is more neutral something close to 1/N^3 is used.
    ! Not sure    
    photon_loss(:)=photon_loss_all(:)/(real(mesh(1))*real(mesh(2))*real(mesh(3)))
#ifdef MPILOG
    write(logf,*) photon_loss(:)
#endif

    ! Report minimum value of xh_av(0) to check for zeros
    if (rank == 0) then
       write(logf,*) "min xh_av: ",minval(xh_av(:,:,:,0))
       write(logf,*) "min xhe_av: ",minval(xhe_av(:,:,:,0)) !*** he+ and he++ start with 0 so no sense to check
       flush(logf)
    endif   
#ifdef MPILOG
       write(logf,*) "min xh_av: ",minval(xh_av(:,:,:,0))
       write(logf,*) "min xhe_av: ",minval(xhe_av(:,:,:,0)) !*** he+ and he++ start with 0 so no sense to check
       flush(logf)
#endif
    ! Apply total photo-ionization rates from all sources (phih_grid)
    conv_flag=0 ! will be used to check for convergence
    
#ifdef MPI
    call MPI_BARRIER(MPI_COMM_NEW,mympierror)
#endif
    ! Loop through the entire mesh

    if (rank == 0) write(logf,*) 'Doing global '
    do k=1,mesh(3)
       do j=1,mesh(2)
          do i=1,mesh(1)
             pos=(/ i,j,k /)
             call evolve0D_global(dt,pos,conv_flag) 
          enddo
       enddo
    enddo     

    ! Report on convergence and intermediate result
    if (rank == 0) then
       write(logf,*) "Number of non-converged points: ",conv_flag
       write(logf,*) "Intermediate result for mean H ionization fraction: ", &
            sum(xh_intermed(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3))
       write(logf,*) "Intermediate result for mean He(+,++) ionization fraction: ", &
            sum(xhe_intermed(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3)),sum(xhe_intermed(:,:,:,2))/real(mesh(1)*mesh(2)*mesh(3))
    endif   
    
!    write(*,*) 'almost at end of global pass'
        ! Report on photon conservation
    call calculate_photon_statistics (dt,xh_intermed,xh_av,xhe_intermed,xhe_av)  
    call report_photonstatistics (dt) 
!    write(*,*) 'at the end of global pass'
  end subroutine global_pass

  ! ===========================================================================
  
  !> Ray tracing the entire grid for all the sources using the
  !! master-slave model for distributing the sources over the
  !! MPI processes.
  subroutine do_grid_master_slave (dt,niter)

    ! Ray tracing the entire grid for all the sources using the
    ! master-slave model for distributing the sources over the
    ! MPI processes.

    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    integer,intent(in) :: niter !< interation counter, passed on to evolve0D

    if (rank == 0) then
       call do_grid_master ()
    else
       call do_grid_slave (dt,niter)
    endif

  end subroutine do_grid_master_slave

  ! ===========================================================================

  !> The master task in the master-slave setup for distributing
  !! the ray-tracing over the sources over the MPI processes.
  subroutine do_grid_master ()

    ! The master task in the master-slave setup for distributing
    ! the ray-tracing over the sources over the MPI processes.

    integer :: ns1
    integer :: sources_done,whomax,who,answer
    ! counter for master-slave process
    integer,dimension(:),allocatable :: counts
#ifdef MPI
    integer :: mympierror
#endif

#ifdef MPI
    ! Source Loop - Master Slave with rank=0 as Master
    sources_done = 0
          
    ns1 = 0
    
    ! Allocate counter for master-slave process
    if (.not.(allocated(counts))) allocate(counts(0:npr-1))

    ! send tasks to slaves 
    
    whomax = min(NumSrc,npr-1)
    do who=1,whomax
       if (ns1 <= NumSrc) then
          ns1=ns1+1
          call MPI_Send (ns1, 1, MPI_INTEGER, who, 1,  &
               MPI_COMM_NEW, mympierror)
       endif
    enddo
    
    do while (sources_done < NumSrc)
       
       ! wait for an answer from a slave. 
       
       call MPI_Recv (answer,     & ! address of receive buffer
            1,		   & ! number of items to receive
            MPI_INTEGER,	   & ! type of data
            MPI_ANY_SOURCE,  & ! can receive from any other
            1,		   & ! tag
            MPI_COMM_NEW,	   & ! communicator
            mympi_status,	   & ! status
            mympierror)
       
       who = mympi_status(MPI_SOURCE) ! find out who sent us the answer
       sources_done=sources_done+1 ! and the number of sources done
       
       ! put the slave on work again,
       ! but only if not all tasks have been sent.
       ! we use the value of num to detect this */
       if (ns1 < NumSrc) then
          ns1=ns1+1
          call MPI_Send (ns1, 1, MPI_INTEGER, &
               who,		&	
               1,		&	
               MPI_COMM_NEW, &
               mympierror)
       endif
    enddo
    
    ! Now master sends a message to the slaves to signify that they 
    ! should end the calculations. We use a special tag for that:
    
    do who = 1,npr-1
       call MPI_Send (0, 1, MPI_INTEGER, &
            who,			  &
            2,			  & ! tag 
            MPI_COMM_NEW,	          &
            mympierror)
       
       ! the slave will send to master the number of calculations
       ! that have been performed. 
       ! We put this number in the counts array.
       
       call MPI_Recv (counts(who), & ! address of receive buffer
            1,                & ! number of items to receive
            MPI_INTEGER,      & ! type of data 
            who,              & ! receive from process who 
            7,                & ! tag 
            MPI_COMM_NEW,     & ! communicator 
            mympi_status,     & ! status
            mympierror)
    enddo
    
    write(logf,*) 'Mean number of sources per processor: ', &
         real(NumSrc)/real(npr-1)
    write(logf,*) 'Counted mean number of sources per processor: ', &
         real(sum(counts(1:npr-1)))/real(npr-1)
    write(logf,*) 'Minimum and maximum number of sources ', &
         'per processor: ', &
         minval(counts(1:npr-1)),maxval(counts(1:npr-1))
    flush(logf)

#endif

  end subroutine do_grid_master

  ! ===========================================================================

  !> The slave task in the master-slave setup for distributing
  !! the ray-tracing over the sources over the MPI processes.
  subroutine do_grid_slave(dt,niter)

    ! The slave task in the master-slave setup for distributing
    ! the ray-tracing over the sources over the MPI processes.

    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    integer,intent(in) :: niter !< interation counter, passed on to evolve0D

    integer :: local_count
    integer :: ns1
#ifdef MPI
    integer :: mympierror
#endif

#ifdef MPI
    local_count=0
    call MPI_Recv (ns1,  & ! address of receive buffer
         1,    & ! number of items to receive
         MPI_INTEGER,  & ! type of data
         0,		  & ! can receive from master only
         MPI_ANY_TAG,  & ! can expect two values, so
         ! we use the wildcard MPI_ANY_TAG 
         ! here
         MPI_COMM_NEW, & ! communicator
         mympi_status, & ! status
         mympierror)
    
    ! if tag equals 2, then skip the calculations
    
    if (mympi_status(MPI_TAG) /= 2) then
       do 
#ifdef MPILOG
          ! Report
          write(logf,*) 'Processor ',rank,' received: ',ns1
          write(logf,*) ' that is source ',SrcSeries(ns1)
          write(logf,*) ' at:',srcpos(:,ns1)
          flush(logf)
#endif
          ! Do the source at hand
          call do_source(dt,ns1,niter)
          
          ! Update local counter
          local_count=local_count+1
          
#ifdef MPILOG
          ! Report ionization fractions
          write(logf,*) sum(xh_av(:,:,:,1))/ & !this was _intermed
               real(mesh(1)*mesh(2)*mesh(3))
          write(logf,*) sum(xhe_av(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3))
          write(logf,*) sum(xhe_av(:,:,:,2))/real(mesh(1)*mesh(2)*mesh(3))
          write(logf,*) local_count
#endif
          ! Send 'answer'
          call MPI_Send (local_count, 1,  & ! sending one int 
               MPI_INTEGER, 0, & ! to master
               1,              & ! tag
               MPI_COMM_NEW,   & ! communicator
               mympierror)
          
          ! Receive new source number
          call MPI_Recv (ns1,     & ! address of receive buffer
               1,            & ! number of items to receive
               MPI_INTEGER,  & ! type of data
               0,            & ! can receive from master only
               MPI_ANY_TAG,  & !  can expect two values, so
               !  we use the wildcard MPI_ANY_TAG 
               !  here
               MPI_COMM_NEW, & ! communicator
               mympi_status, & ! status
               mympierror)
          
          ! leave this loop if tag equals 2
          if (mympi_status(MPI_TAG) == 2) then
#ifdef MPILOG
             write(logf,*) 'Stop received'
             flush(logf)
#endif
             exit 
          endif
       enddo
    endif
    
    ! this is the point that is reached when a task is received with
    ! tag = 2
    
    ! send the number of calculations to master and return
    
#ifdef MPILOG
    ! Report
    write(logf,*) 'Processor ',rank,' did ',local_count,' sources'
    flush(logf)
#endif
    call MPI_Send (local_count,  &
         1,           & 
         MPI_INTEGER, & ! sending one int
         0,           & ! to master
         7,           & ! tag
         MPI_COMM_NEW,& ! communicator
         mympierror)
#endif

  end subroutine do_grid_slave

  ! ===========================================================================

  !> Does the ray-tracing over the sources by distributing
  !! the sources evenly over the available MPI processes-
  subroutine do_grid_static (dt,niter)

    ! Does the ray-tracing over the sources by distributing
    ! the sources evenly over the available MPI processes-
    
    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    integer,intent(in) :: niter !< interation counter, passed on to evolve0D

    integer :: ns1
#ifdef MPILOG
    write(logf,*) rank,NumSrc,npr
    flush(logf)
#endif
    ! Source Loop - distributed for the MPI nodes
    do ns1=1+rank,NumSrc,npr
#ifdef MPILOG
       ! Report
       write(logf,*) 'Processor ',rank,' received: ',ns1
       write(logf,*) ' that is source ',SrcSeries(ns1)
       write(logf,*) ' at:',srcpos(:,ns1)
       flush(logf)
#endif
    !   if (rank == 0) write(logf,*) 'doing source', ns1 
       call do_source(dt,ns1,niter)
    enddo

  end subroutine do_grid_static

  ! ===========================================================================
  
  !> Does the ray-tracing over the entire 3D grid for one source.
  !! The number of this source in the current list is ns1.
  subroutine do_source(dt,ns1,niter)

    ! Does the ray-tracing over the entire 3D grid for one source.
    ! The number of this source in the current list is ns1.

    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    integer, intent(in) :: ns1 !< number of the source being done
    integer,intent(in) :: niter !< interation counter, passed on to evolve0D

    integer :: naxis,nplane,nquadrant
    integer :: ns
    integer :: k
    integer :: nbox
    integer ::nnt
    
    ! Mesh position of the cell being treated
    integer,dimension(Ndim) :: rtpos
      
    ! Pick up source number from the source list
    ns=SrcSeries(ns1)
#ifdef MPILOG
    write(logf,*) 'initialize things'
#endif
    ! reset column densities for new source point
    ! coldensh_out is unique for each source point
    coldensh_out(:,:,:)=0.0
    coldenshe_out(:,:,:,:)=0.0  
    ! Find the mesh position for the end points of the loop
    ! We trace until we reach max_subbox (set in c2ray_parameters)
    ! or the end of the grid. In the periodic case the end of the
    ! grid is always mesh/2 away from the source. If the grid is
    ! even-sized we trave mesh/2 cells to the left and mesh/2-1
    ! cell to the right. If it is odd, it is mesh/2 in either direction.
    ! The mod(mesh,2) takes care of handling this.    
    if (periodic_bc) then
       lastpos_r(:)=srcpos(:,ns)+min(max_subbox,mesh(:)/2-1+mod(mesh(:),2))
       lastpos_l(:)=srcpos(:,ns)-min(max_subbox,mesh(:)/2)
    else
       lastpos_r(:)=min(srcpos(:,ns)+max_subbox,mesh(:))
       lastpos_l(:)=max(srcpos(:,ns)-max_subbox,1)
    endif

    ! Loop through grid in the order required by 
    ! short characteristics
    
    ! Transfer is done in a set of cubes of increasing size.
    ! If the HII region is small we do not waste time calculating
    ! column densities of parts of the grid where no radiation
    ! penetrates. To test whether the current subbox is large
    ! enough we use the photon_loss_src. If this is non-zero,
    ! photons are leaving this subbox and we need to do another
    ! one. We also stop once we have done the whole grid.
    nbox=0 ! subbox counter
    photon_loss_src(:)=NormFlux(ns)*S_star_nominal+NormFluxPL(ns) !-1.0 ! to pass the first while test
    last_r(:)=srcpos(:,ns)!(/3,3,3/) ! to pass the first while test !
    last_l(:)=srcpos(:,ns)!(/3,3,3/) ! to pass the first while test srcpos(:,ns)

    ! Loop through boxes of increasing size
    ! NOTE: make this limit on the photon_loss a fraction of
    ! a source flux loss_fraction*NormFlux(ns)*S_star_nominal)
#ifdef MPILOG
    write(logf,*) 'before do loop'
#endif
!    write(*,*) 'I am just before the do while loop'
 !   write(*,*) '#,Flux, s*', ns,NormFlux(ns),S_star_nominal
    
    do while (all(photon_loss_src(:) > 1e-6*NormFlux(ns)*S_star_nominal) &
    !do while (all(photon_loss_src(:) /= 0.0) &
         .and. last_r(3) .le. lastpos_r(3) &
         .and. last_l(3) .ge. lastpos_l(3).and.nbox.lt.1)
       nbox=nbox+1 ! increase subbox counter
       photon_loss_src(:) = 0.0 ! reset photon_loss_src to zero
       photon_loss_src_thread(:,:) = 0.0 ! reset photon_loss_src to zero       
       last_r(:)=min(srcpos(:,ns)+subboxsize*nbox,lastpos_r(:))
       last_l(:)=max(srcpos(:,ns)-subboxsize*nbox,lastpos_l(:))

       ! OpenMP: if we have multiple OpenMP threads (nthreads > 1) we 
       ! parallelize over the threads by doing independent parts of
       ! the mesh.

#ifdef MY_OPENMP 
       tn=omp_get_thread_num()+1 
#else     
       tn=1    
#endif 
       
 !      write(*,*) 'nthreads',nthreads
       if (nthreads > 1) then ! OpenMP parallelization

          ! First do source point (on first pass)
          if (nbox == 1) then
             rtpos(:)=srcpos(:,ns)
             call evolve0D(dt,rtpos,ns,niter)
          endif
#ifdef MPILOG
    		write(logf,*) 'after source point'
#endif
          ! do independent areas of the mesh in parallel using OpenMP
          !$omp parallel default(shared) private(tn)
          !!!reduction(+:photon_loss_src)

          ! Find out your thread number
#ifdef MY_OPENMP
          tn=omp_get_thread_num()+1
#else
          tn=1
#endif
          
          ! Then do the the axes
          !$omp do schedule(dynamic,1)
          do naxis=1,6
             call evolve1D_axis(dt,ns,niter,naxis)
          enddo
          !$omp end do
#ifdef MPILOG
    		write(logf,*) 'after axis'
#endif
          ! Then the source planes
          !$omp do schedule (dynamic,1)
          do nplane=1,12
             call evolve2D_plane(dt,ns,niter,nplane)
          end do
          !$omp end do
#ifdef MPILOG
!    		write(logf,*) 'after planes'
#endif          
          ! Then the quadrants
          !$omp do schedule (dynamic,1)
          do nquadrant=1,8
             call evolve3D_quadrant(dt,ns,niter,nquadrant)
          end do
          !$omp end do
#ifdef MPILOG
    		write(logf,*) 'after quadrants'
#endif          
          !$omp end parallel
          ! Collect photon losses for each thread
          do nnt=1,nthreads
             photon_loss_src(:)=photon_loss_src(:) + &
                  photon_loss_src_thread(:,nnt)
          enddo

       else ! No OpenMP parallelization
          !write(logf,*) 'no omp'
          ! 1. transfer in the upper part of the grid 
          !    (srcpos(3)-plane and above)
!          write(*,*) 'in do source before evolve2D'

#ifdef MPILOG
    		write(logf,*) 'before evolve 2d'
#endif            

          do k=srcpos(3,ns),last_r(3)
             rtpos(3)=k
             call evolve2D(dt,rtpos,ns,niter)
          end do
          
          ! 2. transfer in the lower part of the grid (below srcpos(3))
          do k=srcpos(3,ns)-1,last_l(3),-1
             rtpos(3)=k
             call evolve2D(dt,rtpos,ns,niter)
          end do
          
#ifdef MPILOG
    		write(logf,*) 'after evolve 2d'
#endif          
          ! No OpenMP threads so we use position 1
          photon_loss_src(:)=photon_loss_src_thread(:,1) 
                   
       endif
    enddo

    ! Record the final photon loss, this is the photon loss that leaves
    ! the grid.
    photon_loss(:)=photon_loss(:) + photon_loss_src(:)
#ifdef MPILOG
    !		write(logf,*) 'after do loop'
#endif  
    ! Sum the total number of subboxes used for reporting later
    sum_nbox=sum_nbox+nbox

  end subroutine do_source

  ! ===========================================================================

  !> Traverse a z-plane (z=rtpos(3)) by sweeping in the x and y
  !! directions.
  subroutine evolve2D(dt,rtpos,ns,niter)

    ! Traverse a z-plane (z=rtpos(3)) by sweeping in the x and y
    ! directions.
    
    real(kind=dp),intent(in) :: dt      !! passed on to evolve0D
    integer,dimension(Ndim),intent(inout) :: rtpos !< mesh position, pos(3) is
                                                 !! intent(in)
    integer,intent(in) :: ns           !< current source
    integer,intent(in) :: niter        !< passed on to evolve0D

    integer :: i,j ! mesh positions
 !   write(logf,*) 'in evolve2D'
    ! sweep in `positive' j direction
    
    
    do j=srcpos(2,ns),last_r(2)
       rtpos(2)=j
       do i=srcpos(1,ns),last_r(1)
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) ! `positive' i
       end do
       do i=srcpos(1,ns)-1,last_l(1),-1
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) ! `negative' i
          
       end do
    end do
    
    ! sweep in `negative' j direction
    do j=srcpos(2,ns)-1,last_l(2),-1
       rtpos(2)=j
       do i=srcpos(1,ns),last_r(1)
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) ! `positive' i
       end do
       do i=srcpos(1,ns)-1,last_l(1),-1
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) ! `negative' i
       end do
    end do

  end subroutine evolve2D

  ! ===========================================================================

  ! Ray tracing for the axes going through the source point
  ! should be called after having done the source point
  subroutine evolve1D_axis(dt,ns,niter,naxis)

    real(kind=dp),intent(in) :: dt      ! passed on to evolve0D
    integer,intent(in) :: ns           ! current source
    integer,intent(in) :: niter        ! passed on to evolve0D
    integer,intent(in) :: naxis        ! axis to do

    integer :: i,j,k
    integer,dimension(Ndim) :: rtpos ! mesh position

    select case (naxis)
    case(1)
       ! sweep in +i direction
       rtpos(2:3)=srcpos(2:3,ns)
       do i=srcpos(1,ns)+1,last_r(1)
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) !# `positive' i
       enddo
    case(2)
       ! sweep in -i direction
       rtpos(2:3)=srcpos(2:3,ns)
       do i=srcpos(1,ns)-1,last_l(1),-1
          rtpos(1)=i
          call evolve0D(dt,rtpos,ns,niter) !# `negative' i
       end do
    case(3)
       ! sweep in +j direction
       rtpos(1)=srcpos(1,ns)
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)+1,last_r(2)
          rtpos(2)=j
          call evolve0D(dt,rtpos,ns,niter) !# `positive' j
       end do
    case(4)
       ! sweep in -j direction
       rtpos(1)=srcpos(1,ns)
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)-1,last_l(2),-1
          rtpos(2)=j
          call evolve0D(dt,rtpos,ns,niter) !# `negative' j
       end do
    case(5)
       ! sweep in +k direction
       rtpos(1:2)=srcpos(1:2,ns)
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          call evolve0D(dt,rtpos,ns,niter) !# `positive' k
       end do
    case(6)
       ! sweep in -k direction
       rtpos(1:2)=srcpos(1:2,ns)
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          call evolve0D(dt,rtpos,ns,niter) !# `negative' k
       end do
    end select
    
  end subroutine evolve1D_axis

  ! ===========================================================================

  !> Ray tracing for planes containing the source point
  !! should be called after evolve1D_axis
  subroutine evolve2D_plane(dt,ns,niter,nplane)

    ! find column density for the axes going through the source point
    ! should be called after having done the source point
    
    real(kind=dp),intent(in) :: dt      ! passed on to evolve0D
    integer,intent(in) :: ns           ! current source
    integer,intent(in) :: niter        ! passed on to evolve0D
    integer,intent(in) :: nplane        ! plane to do

    integer :: i,j,k
    integer,dimension(Ndim) :: rtpos ! mesh position

    select case (nplane)
    case(1)
       ! sweep in +i,+j direction
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)+1,last_r(2)
          rtpos(2)=j
          do i=srcpos(1,ns)+1,last_r(1)
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(2)
       ! sweep in +i,-j direction
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)-1,last_l(2),-1
          rtpos(2)=j
          do i=srcpos(1,ns)+1,last_r(1)
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(3)
       ! sweep in -i,+j direction
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)+1,last_r(2)
          rtpos(2)=j
          do i=srcpos(1,ns)-1,last_l(1),-1
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(4)
       ! sweep in -i,-j direction
       rtpos(3)=srcpos(3,ns)
       do j=srcpos(2,ns)-1,last_l(2),-1
          rtpos(2)=j
          do i=srcpos(1,ns)-1,last_l(1),-1
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(5)
       ! sweep in +i,+k direction
       rtpos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do i=srcpos(1,ns)+1,last_r(1)
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(6)
       ! sweep in -i,+k direction
       rtpos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do i=srcpos(1,ns)-1,last_l(1),-1
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(7)
       ! sweep in -i,-k direction
       rtpos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do i=srcpos(1,ns)-1,last_l(1),-1
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(8)
       ! sweep in +i,-k direction
       rtpos(2)=srcpos(2,ns)
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do i=srcpos(1,ns)+1,last_r(1)
             rtpos(1)=i
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(9) 
       ! sweep in +j,+k direction
       rtpos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(10) 
       ! sweep in -j,+k direction
       rtpos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(11) 
       ! sweep in +j,-k direction
       rtpos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
    case(12) 
       ! sweep in -j,-k direction
       rtpos(1)=srcpos(1,ns)
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             call evolve0D(dt,rtpos,ns,niter)
          enddo
       enddo
       
    end select
    
  end subroutine evolve2D_plane

  ! ===========================================================================

  !> Ray tracing for the 8 octants 
  !! should be called after evolve2D_plane
  subroutine evolve3D_quadrant(dt,ns,niter,nquadrant)

    ! find column density for a z-plane srcpos(3) by sweeping in x and y
    ! directions
    
    real(kind=dp),intent(in) :: dt     ! passed on to evolve0D
    integer,intent(in) :: ns           ! current source
    integer,intent(in) :: niter        ! passed on to evolve0D
    integer,intent(in) :: nquadrant    ! which quadrant to do    

    integer :: i,j,k
    integer,dimension(Ndim) :: rtpos ! mesh position

    select case (nquadrant)
    case (1)
       ! sweep in +i,+j,+k direction
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,ns)+1,last_r(1)
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter)
             end do
          enddo
       enddo
    case (2)
       ! sweep in -i,+j,+k direction
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,ns)-1,last_l(1),-1
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case (3)
       ! sweep in +i,-j,+k direction
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,ns)+1,last_r(1)
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case(4)
       ! sweep in -i,-j,+k direction
       do k=srcpos(3,ns)+1,last_r(3)
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,ns)-1,last_l(1),-1
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case (5)
       ! sweep in +i,+j,-k direction
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,ns)+1,last_r(1)
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `positive' i
             end do
          enddo
       enddo
    case (6)
       ! sweep in -i,+j,-k direction
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)+1,last_r(2)
             rtpos(2)=j
             do i=srcpos(1,ns)-1,last_l(1),-1
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case (7)
       ! sweep in +i,-j,-k direction
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,ns)+1,last_r(1)
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    case(8)
       ! sweep in -i,-j,-k direction
       do k=srcpos(3,ns)-1,last_l(3),-1
          rtpos(3)=k
          do j=srcpos(2,ns)-1,last_l(2),-1
             rtpos(2)=j
             do i=srcpos(1,ns)-1,last_l(1),-1
                rtpos(1)=i
                call evolve0D(dt,rtpos,ns,niter) !# `negative' i
             end do
          end do
       enddo
    end select

  end subroutine evolve3D_quadrant

  !=======================================================================

  !> Calculates the photo-ionization rate for one cell due to one source
  !! and adds this contribution to the collective rate.
  subroutine evolve0D(dt,rtpos,ns,niter)
    
    ! Calculates the photo-ionization rate for one cell due to one source
    ! and adds this contribution to the collective rate.
    
    ! Author: Garrelt Mellema
    
    ! Date: 01-Feb-2008 (21-Aug-2006, 20-May-2005, 5-Jan-2005, 02 Jun 2004)
    
    ! Version: multiple sources, fixed temperature
    
    ! Multiple sources
    ! We call this routine for every grid point and for every source (ns).
    ! The photo-ionization rates for each grid point are found and added
    ! to phih_grid, but the ionization fractions are not updated.
    ! For the first pass (niter = 1) it makes sense to DO update the
    ! ionization fractions since this will increase convergence speed
    ! in the case of isolated sources.
    use cgsconstants
    use tped, only: electrondens
    use doric_module, only: doric, coldens
    use radiation, only: photoion, photrates
    use material, only: clumping_point, ionstates
    use c2ray_parameters, only: epsilon,convergence1,convergence2, &
         type_of_clumping, convergence_frac,use_LLS,type_of_LLS
    use mathconstants, only: pi
    use cgsphotoconstants, only: s_H_heth, s_H_heLya,s_He_heLya,sighe0, s_He_he2,sighe1,s_He_he2,s_H_he2
    use material, only: coldensh_LLS, LLS_point
    use photonstatistics, only: total_LLS_loss    
    ! column density for stopping chemistry !***how should this criterion be for including he and more than one freq bands?
    ! for the moment, leave it as it is, it's probably ok. 
    real(kind=dp),parameter :: max_coldensh=2e29!2.0e22_dp!2e19_dp 
    
    logical :: falsedummy ! always false, for tests
    parameter(falsedummy=.false.)

    real(kind=dp),intent(in) :: dt ! time step
    integer,dimension(Ndim),intent(in) :: rtpos ! cell position (for RT)
    integer,intent(in)      :: ns ! source number 
    integer,intent(in)      :: niter ! global iteration number
    
    integer :: nx,nd,idim ! loop counters
    integer,dimension(Ndim) :: pos
    integer,dimension(Ndim) :: srcpos1
    real(kind=dp) :: dist2,path,vol_ph
    real(kind=dp) :: xs,ys,zs
!    real(kind=si) ::temper
    real(kind=dp) :: coldensh_in
    real(kind=dp),dimension(0:1) :: coldenshe_in,coldenshe_out_temp
    real(kind=dp) :: ndens_p
!    real(kind=dp),dimension(0:1) :: yh,yh_av,yh0
    real(kind=dp) :: convergence
    
    type(photrates) :: phi, dummiphi
    type(ionstates) :: ion

    ! set convergence tolerance
    convergence=convergence1
 !   write(logf,*) 'in evolve0D'
    ! Map pos to mesh pos, assuming a periodic mesh
    do idim=1,Ndim
       pos(idim)=modulo(rtpos(idim)-1,mesh(idim))+1
    enddo

    ! If coldensh_out is zero, we have not done this point
    ! yet, so do it. Otherwise do nothing. (grid is set to 0 for every source)
    if (coldensh_out(pos(1),pos(2),pos(3)) == 0.0) then
       ! Initialize local ionization states to the global ones
       do nx=0,1
          ion%h_av(nx)=max(xh_av(pos(1),pos(2),pos(3),nx),epsilon)
          ion%h(nx)=max(xh_intermed(pos(1),pos(2),pos(3),nx),epsilon)
          ion%h_old(nx)=max(xh(pos(1),pos(2),pos(3),nx),epsilon)         
       enddo

       do nx=0,2
          ion%he_av(nx)=max(xhe_av(pos(1),pos(2),pos(3),nx),epsilon)
          ion%he(nx)=max(xhe_intermed(pos(1),pos(2),pos(3),nx),epsilon)
          ion%he_old(nx)=max(xhe(pos(1),pos(2),pos(3),nx),epsilon)         
       enddo
      
       ! Initialize local density and temperature
       ndens_p=ndens(pos(1),pos(2),pos(3))
!       call get_temperature_point(pos(1),pos(2),pos(3),temper)
       ! Find the column density at the entrance point of the cell (short
       ! characteristics)
       
       if ( all( rtpos(:) == srcpos(:,ns) ) ) then
          ! Do not call cinterp for the source point.
          ! Set coldensh and path by hand
          coldensh_in=0.0
          coldenshe_in(:)=0.0
          path=0.5*dr(1)

          ! Find the distance to the source (average?)
          !dist2=0.5*dr(1) !NOT NEEDED         ! this makes vol=dx*dy*dz
          !vol_ph=4.0/3.0*pi*dist2**3
          vol_ph=dr(1)*dr(2)*dr(3)
          
       else
          
          ! For all other points call cinterp to find the column density
          !do idim=1,Ndim
          !   srcpos1(idim)=srcpos(idim,ns)
          !enddo
          call cinterp(rtpos,srcpos(:,ns),coldensh_in,coldenshe_in(0),coldenshe_in(1),path)

          path=path*dr(1)
          
          ! Find the distance to the source
          xs=dr(1)*real(rtpos(1)-srcpos(1,ns))
          ys=dr(2)*real(rtpos(2)-srcpos(2,ns))
          zs=dr(3)*real(rtpos(3)-srcpos(3,ns))
          dist2=xs*xs+ys*ys+zs*zs
          
          ! Find the volume of the shell this cell is part of 
          ! (dilution factor).
          vol_ph=4.0*pi*dist2*path

          ! Add LLS opacity
          ! GM/110224: previously we added this to coldensh_out,
          ! but this gives funny results for the photo-ionization
          ! rate which is based on the difference between the
          ! in and out column density. To mimick a LLS fog, it
          ! may be better to add it here.
          ! Initialize local LLS (if type of LLS is appropriate)
          if (use_LLS) then
             if (type_of_LLS == 2) call LLS_point (pos(1),pos(2),pos(3))
             coldensh_in = coldensh_in + coldensh_LLS * path/dr(1)
          endif          
       endif
          !write(*,*) 'dr(1),path,vol'
          !write(*,*) dr(1),path,vol_ph
          !pause


       ! Only ray trace and exit. Do not touch the ionization
       ! fractions. They are updated using phih_grid in evolve0d_global
!write(logf,*) 'before niter check'
       ! Only do chemistry if this is the first pass over the sources,
       ! and if column density is below the maximum.
       ! On the first global iteration pass it may be beneficial to assume 
       ! isolated sources, but on later passes the effects of multiple sources 
       ! has to be taken into account. 
       ! Therefore no changes to xh, xh_av, etc. should happen on later passes!

       if (niter == -1 .and. coldensh_in < max_coldensh) then !*** ==1
          dummiphi%h=0.0_dp

          call do_chemistry (dt, ndens_p, ion, dummiphi, &
               coldensh_in,coldenshe_in, path, vol_ph, pos, ns, local=.true.)

          ! Copy ion fractions to global arrays.
          ! This will speed up convergence if
          ! the sources are isolated and only ionizing up.
          ! In other cases it does not make a difference.
          xh_intermed(pos(1),pos(2),pos(3),1)=max(ion%h(1), &
               xh_intermed(pos(1),pos(2),pos(3),1))
          xh_intermed(pos(1),pos(2),pos(3),0)=max(epsilon,1.0- &
               xh_intermed(pos(1),pos(2),pos(3),1))
          xh_av(pos(1),pos(2),pos(3),1)=max(ion%h_av(1), &
               xh_av(pos(1),pos(2),pos(3),1))
          xh_av(pos(1),pos(2),pos(3),0)=max(epsilon,1.0_dp- &
               xh_av(pos(1),pos(2),pos(3),1))

          xhe_intermed(pos(1),pos(2),pos(3),0)=min(ion%he(0), &
               xhe_intermed(pos(1),pos(2),pos(3),0))
          xhe_intermed(pos(1),pos(2),pos(3),2)=max(ion%he(2), &
               xhe_intermed(pos(1),pos(2),pos(3),2))
          xhe_intermed(pos(1),pos(2),pos(3),0)=max(epsilon,1.0_dp- &
               xhe_intermed(pos(1),pos(2),pos(3),0) - &
               xhe_intermed(pos(1),pos(2),pos(3),2))

          xhe_av(pos(1),pos(2),pos(3),0)=min(ion%he_av(0), &
               xhe_av(pos(1),pos(2),pos(3),0))
          xhe_av(pos(1),pos(2),pos(3),2)=max(ion%he_av(2), &
               xhe_av(pos(1),pos(2),pos(3),2))

          xhe_av(pos(1),pos(2),pos(3),1)=max(epsilon,1.0_dp- &
               xhe_av(pos(1),pos(2),pos(3),0) - &
               xhe_av(pos(1),pos(2),pos(3),2))          
       endif  !if niter==!
!    write(logf,*) 'after niter check'      
       ! Add the (time averaged) column density of this cell
       ! to the total column density (for this source)
       ! and add the LLS column density to this.
       ! GM/110224: No! This messes up phi since phi is based
       !  upon the difference between the in and out column density.
       !  Instead add the LLS to coldensh_in, see above
       coldensh_out(pos(1),pos(2),pos(3))=coldensh_in + &
            coldens(path,ion%h_av(0),ndens_p,(1.0_dp-abu_he))

       coldenshe_out(pos(1),pos(2),pos(3),0)=coldenshe_in(0) + &
            coldens(path,ion%he_av(0),ndens_p,abu_he)

       coldenshe_out(pos(1),pos(2),pos(3),1)=coldenshe_in(1) + &
            coldens(path,ion%he_av(1),ndens_p,abu_he)
       
       ! Calculate (photon-conserving) photo-ionization rate
!***
!if (1.0.gt.2.0) then  !*** TEST
       if (coldensh_in < max_coldensh) then

          coldenshe_out_temp=coldenshe_out(pos(1),pos(2),pos(3),:)
          call photoion(phi,coldensh_in,coldensh_out(pos(1),pos(2),pos(3)), &
			 coldenshe_in,coldenshe_out_temp, &
			 vol_ph,ns,ion%h_av(1))

          phi%h=phi%h/(ion%h_av(0)*ndens_p*(1.0_dp-abu_he))
          phi%he(0)=phi%he(0)/(ion%he_av(0)*ndens_p*abu_he)
          phi%he(1)=phi%he(1)/(ion%he_av(1)*ndens_p*abu_he)
          
          !if ( ( rtpos(1) == srcpos(1,ns) ).and. rtpos(2) == srcpos(2,ns) ) then
          !write(*,*) phi%h, ndens_p
          !pause
          !endif
          
          ! GM/110224: Add the factor vol/vol_ph to phi, just as we do for
          ! the photon losses.
          ! GM/110302: Use phi%h_in (not out) here since this is where we draw the
          ! photons from, above.          
          if (use_LLS) call total_LLS_loss(phi%h_in*vol/vol_ph, coldensh_LLS * path/dr(1))          
       else
          phi%h=0.0_dp
          phi%h_out=0.0_dp
          phi%he(0)=0.0_dp
          phi%he(1)=0.0_dp
          phi%he_out(0)=0.0_dp
          phi%he_out(1)=0.0_dp
          !GM/130801: not sure what this is useful for
          !phi%int_out(:)=0.0_dp
       endif
       
      
      
       ! Add photo-ionization rate to the global array 
       ! (applied in evolve0D_global)
       phih_grid(pos(1),pos(2),pos(3))= &
            phih_grid(pos(1),pos(2),pos(3))+phi%h
       phihe_grid(pos(1),pos(2),pos(3),:)=&
             phihe_grid(pos(1),pos(2),pos(3),:)+phi%he(:)   
       if (.not. isothermal) phiheat(pos(1),pos(2),pos(3))=phiheat(pos(1),pos(2),pos(3))+phi%hv_h
                                                                                      
                 
!endif !*** TEST       
       ! Photon statistics: register number of photons leaving the grid

       if ( (any(rtpos(:) == last_l(:))) .or. &
            (any(rtpos(:) == last_r(:))) ) then
          photon_loss_src_thread(1,tn)=photon_loss_src_thread(1,tn) + &
               phi%h_out*vol/vol_ph
          !photon_loss_src(1,tn)=photon_loss_src(1,tn) + phi%h_out*vol/vol_ph
       endif

    endif ! end of coldens test

  end subroutine evolve0D

  ! =======================================================================

  !> Calculates the evolution of the hydrogen ionization state for
  !! one cell (mesh position pos) and multiple sources.
  subroutine evolve0D_global(dt,pos,conv_flag)

    ! Calculates the evolution of the hydrogen ionization state for
    ! one cell (pos) and multiple sources.

    ! Author: Garrelt Mellema

    ! Date: 11-Feb-2008 (20-May-2005, 5-Jan-2005, 02 Jun 2004)
    
    ! Version: Multiple sources (global update, no ray tracing)

    ! Multiple sources
    ! Global update: the collected rates are applied and the new ionization 
    ! fractions and temperatures are calculated.
    ! We check for convergence.
    use cgsconstants
    use radiation, only: photrates, scale_int2,scale_int3
    use c2ray_parameters, only: convergence1,convergence2,type_of_clumping
    use c2ray_parameters, only: convergence_frac, convergence_frac2,epsilon
    use material, only: ionstates
    real(kind=dp),intent(in) :: dt ! time step
    integer,dimension(Ndim),intent(in) :: pos ! position on mesh
    integer,intent(inout) :: conv_flag ! convergence counter

    integer :: nx,nit ! loop counters
    real(kind=dp) :: de ! electron density
!    real(kind=dp),dimension(0:1) :: yh,yh_av,yh0 ! ionization fractions
    real(kind=dp) :: yh0_av_old, yhe0_av_old, yhe1_av_old, yhe2_av_old,yh1_av_old 
    real(kind=dp) :: avg_temper, temper ! temperature
    real(kind=dp) :: ndens_p ! local number density
    real(kind=dp) :: temper_old   
    real(kind=dp) :: temper1   
    real(kind=dp) :: temp_av_old,temp_av_new,temper_inter
    

    real(kind=dp) :: phih ! local H photo-ionization rate (only non-zero when local=.false.!)
    real(kind=dp) :: phih_total ! local total photo-ionization rate (including
                                ! photon loss term)
    real(kind=dp),dimension(0:1) :: phihe ! local He photo-ionization rate (only non-zero when local=.false.!)
    real(kind=dp),dimension(0:1) :: phihe_total ! local total photo-ionization rate (including
                                ! photon loss term)
    real(kind=dp),dimension(0:1) :: temp !***
    real(kind=dp) :: convergence
    type(ionstates) :: ion    
    type(photrates) :: phi 
    ! Set convergence tolerance
    convergence=convergence2
!    write(*,*) 'in evolve0D_global'

#ifdef MPILOG
 !   write(logf,*) 'pos in global: ',pos
#endif
    ! Initialize local ionization states to global ones
    do nx=0,1
       ion%h(nx)=max(epsilon,xh_intermed(pos(1),pos(2),pos(3),nx))
       ion%h_old(nx)=max(epsilon,xh(pos(1),pos(2),pos(3),nx))
       ion%h_av(nx)=max(epsilon,xh_av(pos(1),pos(2),pos(3),nx)) ! use calculated xh_av
    enddo

    do nx=0,2
       ion%he(nx)=max(epsilon,xhe_intermed(pos(1),pos(2),pos(3),nx))
       ion%he_old(nx)=max(epsilon,xhe(pos(1),pos(2),pos(3),nx))
       ion%he_av(nx)=max(epsilon,xhe_av(pos(1),pos(2),pos(3),nx)) ! use calculated xhe_av
    enddo

    ! Initialize local scalars for density and temperature
    ndens_p=ndens(pos(1),pos(2),pos(3))

!    avg_temper=temper
 !   write(*,*) 'after initializing ion and dens'

    ! Use the collected photo-ionization rates
    phi%h=phih_grid(pos(1),pos(2),pos(3))
    phi%he(:)=phihe_grid(pos(1),pos(2),pos(3),:)
    if(.not.isothermal) phi%hv_h=phiheat(pos(1),pos(2),pos(3))
 !   write(*,*) 'after initializing rates'    
   ! if(pos(1)==1.and.pos(2)==1.and.pos(3)==1) write(*,*) 'global',phi%hv_h   

   call get_temperature_point (pos(1),pos(2),pos(3),temper_inter,temp_av_old,temper_old)
!    write(*,*) 'after temperature point'       
!    write(*,*) phi%he(1),phi%he(0),phi%h
   if (phi%he(1).le.1.0e-100_dp .and. phi%he(0).le.1.0e-100_dp .and. phi%h.le.1.0e-100_dp) then
   else
   
       


#ifdef MPILOG
    write(logf,*) 'pos in global: ',pos
    write(logf,*)    ndens_p, ion, phi
#endif
    
    ! I think instead of calling here twice get_temp, it is perhaps better to pass t_new
    ! and t_old as arguments from/to do_chemistry. (?)

    temp=(/0.0_dp,0.0_dp/) !this is an empty array for the column density of helium ***
    call do_chemistry (dt, ndens_p, ion, phi,0.0_dp, &
         temp,  1.0_dp, 0.0_dp, pos, 0 , local=.false.)
         
 !   write(logf,*) 'after do chemistry'
    
#ifdef MPILOG
    write(logf,*) 'pos in global: ',pos
    flush(logf)
#endif

     endif  ! if (phi%he(1).le.1.0e-170_dp.and.phi%he(0).le.1.0e-170_dp.and.phi%h.le.1.0e-170_dp) then
     
    ! Test for global convergence using the time-averaged neutral fraction.
    ! For low values of this number assume convergence
    yh0_av_old=xh_av(pos(1),pos(2),pos(3),0) ! use previously calculated xh_av
    yh1_av_old=xh_av(pos(1),pos(2),pos(3),1)
    yhe0_av_old=xhe_av(pos(1),pos(2),pos(3),0)
    yhe1_av_old=xhe_av(pos(1),pos(2),pos(3),1)
    yhe2_av_old=xhe_av(pos(1),pos(2),pos(3),2)
    call get_temperature_point (pos(1),pos(2),pos(3),temper_inter,temp_av_new,temper_old)
!    write(*,*) 'after temperature point 2'  
!    write(*,*) ion%h_av(0),ion%he_av(0),ion%he_av(2),temp_av_new
    if ( (abs((ion%h_av(0)-yh0_av_old)) > convergence2                .and. &
          abs((ion%h_av(0)-yh0_av_old)/ion%h_av(0)) > convergence2   .and. &
              (ion%h_av(0) > convergence_frac)  ).or.                       &
         (abs((ion%he_av(0)-yhe0_av_old)) > convergence2              .and. &
          abs((ion%he_av(0)-yhe0_av_old)/ion%he_av(0)) > convergence2 .and. &
              (ion%he_av(0) > convergence_frac)  ).or.                      &
    !     (abs((ion%he_av(1)-yhe1_av_old)) > convergence2              .and. &
    !      abs((ion%he_av(1)-yhe1_av_old)/ion%he_av(1)) > convergence2 .and. &
    !          (ion%he_av(1) > convergence_frac)  ).or.                      & 
         (abs((ion%he_av(2)-yhe2_av_old)) > convergence2              .and. &
          abs((ion%he_av(2)-yhe2_av_old)/ion%he_av(2)) > convergence2 .and. &
              (ion%he_av(2) > convergence_frac)  ).or.                & 
         !(abs((temper1-temper_old)/temper1) > 1.0e-1_dp).and.              &
         !(abs(temper1-temper_old) >     100.0_dp)                          &
         (abs((temp_av_old-temp_av_new)/temp_av_new) > 1.0e-1_dp).and.              &
         (abs(temp_av_new-temp_av_old) >     100.0_dp)                          &                  
                                                                      ) then
       conv_flag=conv_flag+1
       if ((abs((temp_av_old-temp_av_new)/temp_av_new) > 1.0e-1_dp).and.              &
          (abs(temp_av_new-temp_av_old) > 100.0_dp).and.pos(2).eq.150.and.pos(3).eq.150) then
       !	write(logf,*),pos
       !	write(logf,*) temp_av_new,abs(temp_av_new-temp_av_old)
       endif
    endif
!    write(*,*) 'after conv check' 
! if (1.0.gt.2.0) then !*** TEST
    ! Copy ion fractions to the global arrays.
    do nx=0,1
       xh_intermed(pos(1),pos(2),pos(3),nx)=ion%h(nx)
       xh_av(pos(1),pos(2),pos(3),nx)=ion%h_av(nx)
    enddo

    do nx=0,2
       xhe_intermed(pos(1),pos(2),pos(3),nx)=ion%he(nx)
       xhe_av(pos(1),pos(2),pos(3),nx)=ion%he_av(nx)
    enddo
    
 !  if (.not.isothermal) call set_temperature_point (pos(1),pos(2),pos(3),temper1,av)
    
!  endif !*** TEST

  end subroutine evolve0D_global

  ! ===========================================================================
  subroutine do_chemistry (dt, ndens_p, ion, &
       phi, coldensh_in,coldenshe_in, path, vol_ph, pos, ns, local)

    use c2ray_parameters, only: convergence1,convergence2,type_of_clumping, & 
         convergence_frac, add_photon_losses, epsilon, convergence_frac2
    use tped, only: electrondens
    use doric_module, only: doric, coldens
    use cgsphotoconstants, only: s_H_heth, s_H_heLya,s_He_heLya,sighe0, s_He_he2,sighe1,s_He_he2,s_H_he2
    use material, only: clumping_point, ionstates
    use radiation, only: photoion, photrates, scale_int2,scale_int3

    real(kind=dp),intent(in) :: dt !< time step
    real(kind=dp),intent(in) :: ndens_p
!    real(kind=dp),dimension(0:1),intent(out) :: yh
!    real(kind=dp),dimension(0:1),intent(inout) :: yh_av
!    real(kind=dp),dimension(0:1),intent(in) :: yh0    
!    real(kind=dp),intent(inout) :: phi !< local photo-ionization rate
    real(kind=dp),intent(in) :: coldensh_in
    real(kind=dp),dimension(0:1),intent(in) :: coldenshe_in
    real(kind=dp),intent(in) :: path
    real(kind=dp),intent(in) :: vol_ph
    integer,dimension(Ndim),intent(in) :: pos !< position on mesh
    integer,intent(in)      :: ns !< source number 
    logical,intent(in) :: local !< true if doing a non-global calculation.

    real(kind=dp) :: avg_temper, temper0, temper1,temper2,temper_inter
    real(kind=dp) :: phih_cell, phihv_cell
    real(kind=dp),dimension(0:1) :: phihe_cell
    real(kind=dp) :: yh0_av_old,oldhe1av,oldhe0av,oldhav
    real(kind=dp) :: yh1_av_old
    real(kind=dp) :: yhe0_av_old
    real(kind=dp) :: yhe1_av_old
    real(kind=dp) :: yhe2_av_old
    real(kind=dp) :: de
    real(kind=dp) :: coldensh_cell
    real(kind=dp), dimension(0:1) :: coldenshe_cell, coldensheout
    real(kind=dp) :: yfrac 
    real(kind=dp) :: zfrac
    real(kind=dp) :: y2afrac 
    real(kind=dp) :: y2bfrac 
    real(kind=dp) :: tau_H_heth
    real(kind=dp) :: tau_He_heth
    real(kind=dp) :: tau_H_heLya
    real(kind=dp) :: tau_He_heLya 
    real(kind=dp) :: tau_He2_he2th
    real(kind=dp) :: tau_He_he2th
    real(kind=dp) :: tau_H_he2th

    real(kind=dp) :: Nhe0,Nhe1,Nh0, ovtotneutral
    real(kind=dp) :: N_He0_ov_H0,N_He1_ov_H0,N_He0_ov_He1,N_H0_ov_He1,N_H0_ov_He0,N_He1_ov_He0,fH_2a,fH_2b,fHe0_2a,fHe0_2b
    real(kind=dp) :: fH_3a,fH_3b,fH_3c,fH_3d,fHe0_3a,fHe0_3b,fHe0_3c,fHe0_3d,fHe1_3a,fHe1_3b,fHe1_3c,fHe1_3d
    real(kind=dp) :: ionh0old,ionh1old,ionhe0old,ionhe1old,ionhe2old
    integer :: nx
    integer :: nit

    type(photrates) :: phi
    type(ionstates) :: ion

#ifdef MPILOG
!    if (.not.(local)) write(logf,*) 'pos in do_chem: ',pos
#endif
    ! Initialize yh to initial value (*** not needed, everything in ion) 
    !yh(:)=yh0(:)

    ! Initialize local temperature
    !

    call get_temperature_point (pos(1),pos(2),pos(3),temper_inter,avg_temper,temper1)
    !avg_temper=temper1
    temper0   =temper1
  
    !write(logf,*) 'in do chemistry'
    
    phih_cell=phi%h
    phihe_cell(:)=phi%he(:)
    phihv_cell=phi%hv_h

    ! Initialize local clumping (if type of clumping is appropriate)
    if (type_of_clumping == 5) call clumping_point (pos(1),pos(2),pos(3))
  
   
  
    nit=0
    do 
       nit=nit+1
       temper2   =temper1  ! This is the temperature from last iteration
       ! Save the values of yh_av found in the previous iteration

    !write(logf,*) 'in do chemistry iteration loop', nit

       yh0_av_old=ion%h_av(0)
       yh1_av_old=ion%h_av(1)
       yhe0_av_old=ion%he_av(0)
       yhe1_av_old=ion%he_av(1)
       yhe2_av_old=ion%he_av(2)

       ! Copy ionic abundances back to initial values (doric assumes
       ! that it contains this) !*** this is not longer needed because I have the
       ! *** old values in ion%
       !yh(:)=yh0(:)
              
       ! Calculate (mean) electron density
       de=electrondens(ndens_p,ion%h_av,ion%he_av)

       ! Find total photo-ionization rate 
       if (local) then
          
          ! Calculate (time averaged) column density of cell
          coldensh_cell=coldens(path,ion%h_av(0),ndens_p,1.0_dp-abu_he)
          do nx=0,1
             coldenshe_cell(nx)=coldens(path,ion%he_av(nx),ndens_p,abu_he)
          enddo  
          coldensheout(:)=coldenshe_in(:)+coldenshe_cell(:)
          ! Calculate (photon-conserving) photo-ionization rate
          
          call photoion(phi,coldensh_in,coldensh_in+coldensh_cell, &
               coldenshe_in, coldensheout, &
               vol_ph,ns,ion%h_av(1))

          phi%h=phi%h/(ion%h_av(0)*ndens_p*(1.0_dp-abu_he))

          do nx=0,1
             phi%he(nx)=phi%he(nx)/(ion%he_av(nx)*ndens_p*abu_he)
          enddo
       else

          ! (direct plus photon losses)
          ! DO THIS HERE, yh_av is changing
          ! (if the cell is ionized, add a fraction of the lost photons)
          !if (xh_intermed(pos(1),pos(2),pos(3),1) > 0.5)
          phi%h=phih_cell
          phi%he(:)=phihe_cell
          phi%hv_h=phihv_cell
          if (.not.isothermal) call ini_rec_colion_factors(avg_temper) !initialize the collisional ion
        
          if (add_photon_losses) then               !*** this has to be extended!   

             Nhe0= abu_he*ion%he(0)              ! *** probably use average fractions
             Nhe1= abu_he*ion%he(1)              ! *** probably use average fractions
       	     Nh0 = (1.0-abu_he)*ion%h(0)
             N_He0_ov_H0         = max(epsilon,Nhe0)/max(epsilon,Nh0)
             N_He0_ov_He1        = max(epsilon,Nhe0)/max(epsilon,Nhe1)
             N_H0_ov_He0         = max(epsilon,Nh0)  /max(epsilon,Nhe0)
             N_H0_ov_He1         = max(epsilon,Nh0)  /max(epsilon,Nhe1)
             N_He1_ov_He0        = max(epsilon,Nhe1)/max(epsilon,Nhe0)
             N_He1_ov_H0         = max(epsilon,Nhe1)/max(epsilon,Nh0)
 	     !call scale_int3(fH_3a,fHe0_3a,fHe1_3a, &
             !N_H0_ov_He0,N_He0_ov_H0,N_H0_ov_He1,N_He1_ov_H0,N_He0_ov_He1,N_He1_ov_He0, &
             !1,2,3)
 	     !call scale_int3(fH_3b,fHe0_3b,fHe1_3b, &
             !N_H0_ov_He0,N_He0_ov_H0,N_H0_ov_He1,N_He1_ov_H0,N_He0_ov_He1,N_He1_ov_He0, &
             !4,5,6)
 	     !call scale_int3(fH_3c,fHe0_3c,fHe1_3c, &
             !N_H0_ov_He0,N_He0_ov_H0,N_H0_ov_He1,N_He1_ov_H0,N_He0_ov_He1,N_He1_ov_He0, &
             !7,8,9)
      	     !call scale_int3(fH_3d,fHe0_3d,fHe1_3d, &
             !N_H0_ov_He0,N_He0_ov_H0,N_H0_ov_He1,N_He1_ov_H0,N_He0_ov_He1,N_He1_ov_He0, &
             !10,11,12)
 	     !call scale_int2(fH_2a,fHe0_2a,N_He0_ov_H0,1,2)
   	     !call scale_int2(fH_2b,fHe0_2b,N_He0_ov_H0,3,4)
!*** it is not so clear to me if I should first scale the photon loss in every interval with the amount of 
!*** "neutral" ions those photons could possibly ionize, or if I should first sclae them and ascribe them to 
!*** either H, He0 or He+ and then scale them with only the species in question. 
!*** for now, I'll just take the first option (scale and then devide)
            ovtotneutral= 1.0_dp/(vol*ion%h_av(0)*ndens_p*(1.0-abu_he))
   	    photon_loss(1)=photon_loss(1)*ovtotneutral
            ovtotneutral=1.0_dp/(vol*ndens_p*((1.0-abu_he)*ion%h_av(0)+abu_he*ion%he_av(0)))
   	    photon_loss(2)=photon_loss(2)*ovtotneutral
  	    photon_loss(3)=photon_loss(3)*ovtotneutral
            ovtotneutral= 1.0_dp/(vol*ndens_p*((1.0-abu_he)*ion%h_av(0)+abu_he*(ion%he_av(0)+ion%he_av(1))))
   	    photon_loss(4)=photon_loss(4)*ovtotneutral
   	    photon_loss(5)=photon_loss(5)*ovtotneutral
   	    photon_loss(6)=photon_loss(6)*ovtotneutral
   	    photon_loss(7)=photon_loss(7)*ovtotneutral
   	    phi%h  = phi%h   +(photon_loss(1) +photon_loss(2)* fH_2a  +photon_loss(4)*fH_3a+photon_loss(5)*fH_3b & 
                                              +photon_loss(3)* fH_2b  +photon_loss(6)*fH_3c+photon_loss(7)*fH_3d) 
   	    phi%he(0)= phi%he(0) +(photon_loss(2)*fHe0_2a  +photon_loss(4)* fHe0_3a +photon_loss(5)* fHe0_3b &
                        + photon_loss(3)*fHe0_2b  +photon_loss(6)* fHe0_3c +photon_loss(7)* fHe0_3d)
   	    phi%he(1)= phi%he(1) +(photon_loss(4)*fHe1_3a +photon_loss(5)*fHe1_3b &
                        + photon_loss(6)*fHe1_3c +photon_loss(7)*fHe1_3d)
       endif   ! add photon loss if
    endif     ! local if/else

 !      Calculate the new and mean ionization states
!*** DO THIS ONLY IF PHI IS NOT ==0
  ! if (phi%h.ne.0.0_dp) then
  

  

       coldensh_cell    =coldens(path,ion%h(0),ndens_p,(1.0_dp-abu_he))
       coldenshe_cell(0)=coldens(path,ion%he(0),ndens_p,abu_he)
       coldenshe_cell(1)=coldens(path,ion%he(1),ndens_p,abu_he)
#ifdef MPILOG
!    if (.not.(local)) write(logf,*) 'pos in do_chem: ',pos
#endif

!     write(logf,*) 'coldens' ,   coldensh_cell,coldenshe_cell   

       tau_H_heth  = coldensh_cell*s_H_heth     ! sigh*(frthe0/frth0)**(-sh0) opt depth of H at he0 ion threshold
       tau_He_heth = coldenshe_cell(0)*sighe0     ! sighe0                  ! opt depth of He0 at he0 ion threshold
       tau_H_heLya = coldensh_cell*s_H_heLya    ! sigh*(40.817eV/frth0)**(-sh0)       opt depth of H  at he+Lya 
       tau_He_heLya= coldenshe_cell(0)*s_He_heLya      ! sighe0*(40.817eV/frthe0)**(-she0)   opt depth of He at he+Lya  
       tau_He2_he2th = coldenshe_cell(1)*sighe1 
       tau_He_he2th = coldenshe_cell(0)*s_He_he2
       tau_H_he2th = coldensh_cell*s_H_he2

#ifdef MPILOG
!    if (.not.(local)) write(logf,*) 'pos in do_chem: ',tau_H_heth,tau_He_heth,tau_H_heLya,tau_He_heLya
    flush(logf)
#endif
       yfrac= tau_H_heth /(tau_H_heth +tau_He_heth)
       zfrac= tau_H_heLya/(tau_H_heLya+tau_He_heLya)
       y2afrac=  tau_He2_he2th /(tau_He2_he2th +tau_He_he2th+tau_H_he2th)
       y2bfrac=  tau_He_he2th /(tau_He2_he2th +tau_He_he2th+tau_H_he2th)

#ifdef MPILOG
!    if (.not.(local)) write(logf,*) 'pos in do_chem: ',dt,de,ndens_p,ion,phi,yfrac,zfrac
    flush(logf)
#endif
   
!      write(logf,*) 'before doric', ion,phi

       call doric(dt,de  ,ndens_p,ion,phi,yfrac,zfrac,y2afrac,y2bfrac)!,local)! 
       de=electrondens(ndens_p,ion%h_av,ion%he_av)

#ifdef MPILOG
  !  if (.not.(local)) write(logf,*) 'pos in do_chem: ',pos
#endif
     !  if (ion%he(0).le.1.0e-30.or.ion%h(0).le.1.0e-30.or.ion%he(1).le.1.0e-30) then
         coldensh_cell  =coldens(path,ion%h(0),ndens_p,(1.0_dp-abu_he))
         coldenshe_cell(0)=coldens(path,ion%he(0),ndens_p,abu_he)
         coldenshe_cell(1)=coldens(path,ion%he(1),ndens_p,abu_he)

       tau_H_heth  = coldensh_cell*s_H_heth     ! sigh*(frthe0/frth0)**(-sh0) opt depth of H at he0 ion threshold
       tau_He_heth = coldenshe_cell(0)*sighe0     ! sighe0                  ! opt depth of He0 at he0 ion threshold
       tau_H_heLya = coldensh_cell*s_H_heLya    ! sigh*(40.817eV/frth0)**(-sh0)       opt depth of H  at he+Lya 
       tau_He_heLya= coldenshe_cell(0)*s_He_heLya  ! sighe0*(40.817eV/frthe0)**(-she0)   opt depth of He at he+Lya  
       tau_He2_he2th = coldenshe_cell(1)*sighe1 
       tau_He_he2th = coldenshe_cell(0)*s_He_he2
       tau_H_he2th = coldensh_cell*s_H_he2

         yfrac= tau_H_heth /(tau_H_heth +tau_He_heth)
         zfrac= tau_H_heLya/(tau_H_heLya+tau_He_heLya)
         y2afrac=  tau_He2_he2th /(tau_He2_he2th +tau_He_he2th+tau_H_he2th)
         y2bfrac=  tau_He_he2th /(tau_He2_he2th +tau_He_he2th+tau_H_he2th)

         ionh0old=ion%h(0)
         ionh1old=ion%h(1)
         ionhe0old=ion%he(0)
         ionhe1old=ion%he(1)
         ionhe2old=ion%he(2)
         oldhav=ion%h_av(0)
         oldhe0av=ion%he_av(0)
         oldhe1av=ion%he_av(1)        
         


         call doric(dt,de  ,ndens_p,ion,phi,yfrac,zfrac,y2afrac,y2bfrac)!,local)! 
         ion%h(0)=(ion%h(0)+ionh0old)/2.0_dp
         ion%h(1)=(ion%h(1)+ionh1old)/2.0_dp
         ion%he(0)=(ion%he(0)+ionhe0old)/2.0_dp
         ion%he(1)=(ion%he(1)+ionhe1old)/2.0_dp
         ion%he(2)=(ion%he(2)+ionhe2old)/2.0_dp
         ion%h_av(0)=(ion%h_av(0)+oldhav)/2.0_dp
         ion%he_av(0)=(ion%he_av(0)+oldhe0av)/2.0_dp
         ion%he_av(1)=(ion%he_av(1)+oldhe1av)/2.0_dp
         
         de=electrondens(ndens_p,ion%h_av,ion%he_av)
     !  endif

         temper1=temper0 
               if (.not.isothermal) &
               call thermal(dt,temper1,avg_temper,de,ndens_p, &
               ion,phi)    
       
       ! Test for convergence on time-averaged neutral fraction
       ! For low values of this number assume convergence
       if ((abs((ion%h_av(0)-yh0_av_old)/ion%h_av(0)) < convergence2           &
             .or. (ion%h_av(0) < convergence_frac)).and.                       &

	     !(abs((ion%h_av(1)-yh1_av_old)/ion%h_av(1)) < convergence2         &
             !.or. (ion%h_av(1) < convergence_frac)).and.                       &

	     (abs((ion%he_av(0)-yhe0_av_old)/ion%he_av(0)) < convergence2      &
             .or. (ion%he_av(0) < convergence_frac)).and.                      &

	    ! (abs((ion%he_av(1)-yhe1_av_old)/ion%he_av(1)) < convergence2      &
            ! .or. (ion%he_av(1) < convergence_frac)).and.                      &

             (abs((ion%he_av(2)-yhe2_av_old)/ion%he_av(2)) < convergence2      &
             .or. (ion%he_av(2) < convergence_frac2)) .and.                    &
             
              (abs((temper1-temper2)/temper1) < convergence2) & 
            ) then  
              exit
       endif  

       ! Warn about non-convergence and terminate iteration
       if (nit > 400) then
          if (rank == 0) then   
             write(logf,*) 'Convergence failing (global) nit=', nit
             write(logf,*) 'x',ion%h_av(0), ion%he_av(0:2)
             write(logf,*) 'h',yh0_av_old, yhe0_av_old,yhe1_av_old, yhe2_av_old
             write(logf,*) abs(ion%h_av(0)-yh0_av_old),abs(ion%he_av(0)-yhe0_av_old),abs(ion%he_av(1))
          endif
          exit
       endif
      !else
      !	exit
      !endif   ! if phi.ne.0.0
    enddo

   if (.not. isothermal) call set_temperature_point (pos(1),pos(2),pos(3),temper1,avg_temper)

  end subroutine do_chemistry

  ! ===========================================================================

  !> Finds the column density at pos as seen from the source point srcpos
  !! through interpolation. The interpolation
  !! depends on the orientation of the ray. The ray crosses either
  !! a z-plane, a y-plane or an x-plane.
    subroutine cinterp (pos,srcpos,cdensi,cdensihe0,cdensihe1,path)
  
    ! Author: Garrelt Mellema
    
    ! Date: 21-Mar-2006 (06-Aug-2004)
    
    ! History:
    ! Original routine written by Alex Raga, Garrelt Mellema, Jane Arthur
    ! and Wolfgang Steffen in 1999.
    ! This version: Modified for use with a grid based approach.
    ! Better handling of the diagonals.
    ! Fortran90
    
    ! does the interpolation to find the column density at pos
    ! as seen from the source point srcpos. the interpolation
    ! depends on the orientation of the ray. The ray crosses either
    ! a z-plane, a y-plane or an x-plane.
    
    integer,dimension(Ndim),intent(in) :: pos !< cell position (mesh)
    integer,dimension(Ndim),intent(in) :: srcpos !< source position (mesh)
    real(kind=dp),intent(out) :: cdensi !< column density to cell
    real(kind=dp),dimension(0:1),intent(out) :: cdensihe0 !< column density to cell
    real(kind=dp),intent(out) :: cdensihe1 !< column density to cell
    real(kind=dp),intent(out) :: path !< path length over cell

    real(kind=dp),parameter :: sqrt3=sqrt(3.0)
    real(kind=dp),parameter :: sqrt2=sqrt(2.0)

    integer :: i,j,k,i0,j0,k0

    integer :: idel,jdel,kdel
    integer :: idela,jdela,kdela
    integer :: im,jm,km
    integer :: ip,imp,jp,jmp,kp,kmp
    integer :: sgni,sgnj,sgnk
    real(kind=dp) :: alam,xc,yc,zc,dx,dy,dz,s1,s2,s3,s4
    real(kind=dp) :: c1,c2,c3,c4,c1He0,c2He0,c3He0,c4He0,c1He1,c2He1,c3He1,c4He1
    real(kind=dp) :: dxp,dyp,dzp
    real(kind=dp) :: w1,w2,w3,w4,w1He0,w2He0,w3He0,w4He0,w1He1,w2He1,w3He1,w4He1
    real(kind=dp) :: di,dj,dk


    !DEC$ ATTRIBUTES FORCEINLINE :: weightf
    ! map to local variables (should be pointers ;)
    i=pos(1)
    j=pos(2)
    k=pos(3)
    i0=srcpos(1)
    j0=srcpos(2)
    k0=srcpos(3)
 
    ! calculate the distance between the source point (i0,j0,k0) and 
    ! the destination point (i,j,k)
    idel=i-i0
    jdel=j-j0
    kdel=k-k0
    idela=abs(idel)
    jdela=abs(jdel)
    kdela=abs(kdel)
     
    ! Find coordinates of points closer to source
    sgni=sign(1,idel)
!      if (idel == 0) sgni=0
    sgnj=sign(1,jdel)
!      if (jdel == 0) sgnj=0
    sgnk=sign(1,kdel)
!      if (kdel == 0) sgnk=0
    im=i-sgni
    jm=j-sgnj
    km=k-sgnk
    di=real(idel)
    dj=real(jdel)
    dk=real(kdel)
 
    ! Z plane (bottom and top face) crossing
    ! we find the central (c) point (xc,xy) where the ray crosses 
    ! the z-plane below or above the destination (d) point, find the 
    ! column density there through interpolation, and add the contribution
    ! of the neutral material between the c-point and the destination
    ! point.
    
    if (kdela >= jdela.and.kdela >= idela) then
       
       ! alam is the parameter which expresses distance along the line s to d
       ! add 0.5 to get to the interface of the d cell.
       alam=(real(km-k0)+sgnk*0.5)/dk
              
       xc=alam*di+real(i0) ! x of crossing point on z-plane 
       yc=alam*dj+real(j0) ! y of crossing point on z-plane
       
       dx=2.0*abs(xc-(real(im)+0.5*sgni)) ! distances from c-point to
       dy=2.0*abs(yc-(real(jm)+0.5*sgnj)) ! the corners.
       
       s1=(1.-dx)*(1.-dy)    ! interpolation weights of
       s2=(1.-dy)*dx         ! corner points to c-point
       s3=(1.-dx)*dy
       s4=dx*dy
        
       ip =modulo(i-1, mesh(1))+1
       imp=modulo(im-1,mesh(1))+1
       jp =modulo(j-1, mesh(2))+1
       jmp=modulo(jm-1,mesh(2))+1
       kmp=modulo(km-1,mesh(3))+1
       c1=     coldensh_out(imp,jmp,kmp)    !# column densities at the
       c2=     coldensh_out(ip,jmp,kmp)     !# four corners
       c3=     coldensh_out(imp,jp,kmp)
       c4=     coldensh_out(ip,jp,kmp)

       c1he0=coldenshe_out(imp,jmp,kmp,0)    !# column densities at the
       c2he0=coldenshe_out(ip,jmp,kmp,0)     !# four corners
       c3he0=coldenshe_out(imp,jp,kmp,0)
       c4he0=coldenshe_out(ip,jp,kmp,0)

       c1he1=coldenshe_out(imp,jmp,kmp,1)    !# column densities at the
       c2he1=coldenshe_out(ip,jmp,kmp,1)     !# four corners
       c3he1=coldenshe_out(imp,jp,kmp,1)
       c4he1=coldenshe_out(ip,jp,kmp,1)
       
       ! extra weights for better fit to analytical solution
       w1=   s1*weightf(c1,0)
       w2=   s2*weightf(c2,0)
       w3=   s3*weightf(c3,0)
       w4=   s4*weightf(c4,0)

       w1he0=s1*weightf(c1he0,1)
       w2he0=s2*weightf(c2he0,1)
       w3he0=s3*weightf(c3he0,1)
       w4he0=s4*weightf(c4he0,1)

       w1he1=s1*weightf(c1he1,2)
       w2he1=s2*weightf(c2he1,2)
       w3he1=s3*weightf(c3he1,2)
       w4he1=s4*weightf(c4he1,2)

       ! column density at the crossing point
       cdensi   =(c1   *w1   +c2   *w2   +c3   *w3   +c4   *w4   )/(w1+w2+w3+w4) 
       cdensihe0=(c1he0*w1he0+c2he0*w2he0+c3he0*w3he0+c4he0*w4he0)/(w1he0+w2he0+w3he0+w4he0) 
       cdensihe1=(c1he1*w1he1+c2he1*w2he1+c3he1*w3he1+c4he1*w4he1)/(w1he1+w2he1+w3he1+w4he1) 
 
       ! Take care of diagonals
       ! if (kdela == idela.or.kdela == jdela) then
       ! if (kdela == idela.and.kdela == jdela) then
       ! cdensi=sqrt3*cdensi
       !else
       !cdensi=sqrt2*cdensi
       !endif
       !endif

       if (kdela == 1.and.(idela == 1.or.jdela == 1)) then
          if (idela == 1.and.jdela == 1) then
             cdensi=   sqrt3*cdensi
             cdensihe0=sqrt3*cdensihe0
             cdensihe1=sqrt3*cdensihe1
          else
             cdensi=   sqrt2*cdensi
             cdensihe0=sqrt2*cdensihe0
             cdensihe1=sqrt2*cdensihe1
          endif
       endif
       ! if (kdela == 1) then
       ! if ((w3 == 1.0).or.(w2 == 1.0)) cdensi=sqrt(2.0)*cdensi
       ! if (w1 == 1.0) cdensi=sqrt(3.0)*cdensi
       ! write(logf,*) idela,jdela,kdela
       !endif

       ! Path length from c through d to other side cell.
       !dxp=di/dk
       !dyp=dj/dk
       path=sqrt((di*di+dj*dj)/(dk*dk)+1.0) ! pathlength from c to d point  
 

       ! y plane (left and right face) crossing
       ! (similar approach as for the z plane, see comments there)
    elseif (jdela >= idela.and.jdela >= kdela) then
           
       alam=(real(jm-j0)+sgnj*0.5)/dj
       zc=alam*dk+real(k0)
       xc=alam*di+real(i0)
       dz=2.0*abs(zc-(real(km)+0.5*sgnk))
       dx=2.0*abs(xc-(real(im)+0.5*sgni))
       s1=(1.-dx)*(1.-dz)
       s2=(1.-dz)*dx
       s3=(1.-dx)*dz
       s4=dx*dz
       ip=modulo(i-1,mesh(1))+1
       imp=modulo(im-1,mesh(1))+1
       jmp=modulo(jm-1,mesh(2))+1
       kp=modulo(k-1,mesh(3))+1
       kmp=modulo(km-1,mesh(3))+1

       c1=  coldensh_out(imp,jmp,kmp)
       c1he0=coldenshe_out(imp,jmp,kmp,0)
       c1he1=coldenshe_out(imp,jmp,kmp,1)

       c2=  coldensh_out(ip,jmp,kmp)
       c2he0=coldenshe_out(ip,jmp,kmp,0)
       c2he1=coldenshe_out(ip,jmp,kmp,1)

       c3=  coldensh_out(imp,jmp,kp)
       c3he0=coldenshe_out(imp,jmp,kp,0)
       c3he1=coldenshe_out(imp,jmp,kp,1)

       c4=  coldensh_out(ip,jmp,kp)
       c4he0=coldenshe_out(ip,jmp,kp,0)
       c4he1=coldenshe_out(ip,jmp,kp,1)

       ! extra weights for better fit to analytical solution
       w1=s1*weightf(c1,0)
       w2=s2*weightf(c2,0)
       w3=s3*weightf(c3,0)
       w4=s4*weightf(c4,0)

       w1he0=s1*weightf(c1he0,1)
       w2he0=s2*weightf(c2he0,1)
       w3he0=s3*weightf(c3he0,1)
       w4he0=s4*weightf(c4he0,1)

       w1he1=s1*weightf(c1he1,2)
       w2he1=s2*weightf(c2he1,2)
       w3he1=s3*weightf(c3he1,2)
       w4he1=s4*weightf(c4he1,2)  
     
       cdensi=   (c1   *w1   +c2   *w2   +c3   *w3   +c4   *w4   )/(w1+w2+w3+w4)
       cdensihe0=(c1he0*w1he0+c2he0*w2he0+c3he0*w3he0+c4he0*w4he0)/(w1he0+w2he0+w3he0+w4he0) 
       cdensihe1=(c1he1*w1he1+c2he1*w2he1+c3he1*w3he1+c4he1*w4he1)/(w1he1+w2he1+w3he1+w4he1)        
       ! Take care of diagonals
       if (jdela == 1.and.(idela == 1.or.kdela == 1)) then
          if (idela == 1.and.kdela == 1) then
             !write(logf,*) 'error',i,j,k
             cdensi=   sqrt3*cdensi
             cdensihe0=sqrt3*cdensihe0
             cdensihe1=sqrt3*cdensihe1
          else
             !write(logf,*) 'diagonal',i,j,k
             cdensi=   sqrt2*cdensi
             cdensihe0=sqrt2*cdensihe0
             cdensihe1=sqrt2*cdensihe1
          endif
       endif

       !dxp=di/dj
       !dzp=dk/dj
       !path=sqrt(dxp*dxp+1.0+dzp*dzp)
       path=sqrt((di*di+dk*dk)/(dj*dj)+1.0)
       

       ! x plane (front and back face) crossing
       ! (similar approach as with z plane, see comments there)

    elseif(idela >= jdela.and.idela >= kdela) then

       alam=(real(im-i0)+sgni*0.5)/di
       zc=alam*dk+real(k0)
       yc=alam*dj+real(j0)
       dz=2.0*abs(zc-(real(km)+0.5*sgnk))
       dy=2.0*abs(yc-(real(jm)+0.5*sgnj))
       s1=(1.-dz)*(1.-dy)
       s2=(1.-dz)*dy
       s3=(1.-dy)*dz
       s4=dy*dz
  
       imp=modulo(im-1,mesh(1))+1
       jp= modulo(j-1,mesh(2))+1
       jmp=modulo(jm-1,mesh(2))+1
       kp= modulo(k-1,mesh(3))+1
       kmp=modulo(km-1,mesh(3))+1
       c1=  coldensh_out(imp,jmp,kmp)
       c2=  coldensh_out(imp,jp,kmp)
       c3=  coldensh_out(imp,jmp,kp)
       c4=  coldensh_out(imp,jp,kp)

       c1he0=coldenshe_out(imp,jmp,kmp,0)
       c2he0=coldenshe_out(imp,jp,kmp,0)
       c3he0=coldenshe_out(imp,jmp,kp,0)
       c4he0=coldenshe_out(imp,jp,kp,0)

       c1he1=coldenshe_out(imp,jmp,kmp,1)
       c2he1=coldenshe_out(imp,jp,kmp,1)
       c3he1=coldenshe_out(imp,jmp,kp,1)
       c4he1=coldenshe_out(imp,jp,kp,1)

       ! extra weights for better fit to analytical solution
       w1   =s1*weightf(c1,0)
       w2   =s2*weightf(c2,0)
       w3   =s3*weightf(c3,0)
       w4   =s4*weightf(c4,0)

       w1he0=s1*weightf(c1he0,1)
       w2he0=s2*weightf(c2he0,1)
       w3he0=s3*weightf(c3he0,1)
       w4he0=s4*weightf(c4he0,1)

       w1he1=s1*weightf(c1he1,2)
       w2he1=s2*weightf(c2he1,2)
       w3he1=s3*weightf(c3he1,2)
       w4he1=s4*weightf(c4he1,2)      
 
       cdensi   =(c1   *w1   +c2   *w2   +c3   *w3   +c4   *w4   )/(w1+w2+w3+w4)
       cdensihe0=(c1he0*w1he0+c2he0*w2he0+c3he0*w3he0+c4he0*w4he0)/(w1he0+w2he0+w3he0+w4he0) 
       cdensihe1=(c1he1*w1he1+c2he1*w2he1+c3he1*w3he1+c4he1*w4he1)/(w1he1+w2he1+w3he1+w4he1)       
       if ( idela == 1 .and. ( jdela == 1 .or. kdela == 1 ) ) then
          if ( jdela == 1 .and. kdela == 1 ) then
             cdensi=   sqrt3*cdensi
             cdensihe0=sqrt3*cdensihe0            
             cdensihe1=sqrt3*cdensihe1
          else
             cdensi   =sqrt2*cdensi
             cdensihe0=sqrt2*cdensihe0
             cdensihe1=sqrt2*cdensihe1
          endif
       endif
        
       !dyp=dj/di
       !dzp=dk/di
       !path=sqrt(1.0+dyp*dyp+dzp*dzp)
       path=sqrt(1.0+(dj*dj+dk*dk)/(di*di))
       
    end if
  
  end subroutine cinterp


  ! =========================================================================

  !> Weight function for interpolation in cinterp
  !> Weight function for interpolation in cinterp
  real(kind=dp) function weightf (cd,id)

    use cgsphotoconstants, only: sigh,sighe0,sighe1
    real(kind=dp):: sig
    real(kind=dp),intent(in) :: cd
    integer,intent(in) :: id
    real(kind=dp),parameter :: minweight=1.0_dp/0.6_dp


    !weightf=1.0
    ! weightf=1.0/max(1.0d0,cd**0.54)
    ! weightf=exp(-min(700.0,cd*0.15*6.3d-18))
  if (id==0) then
     sig=sigh
  elseif (id==1) then
     sig=sighe0
  elseif (id==2) then
     sig=sighe1
  endif

    weightf=1.0/max(0.6_dp,cd*sig)

    ! weightf=1.0/log(max(e_ln,cd))

  end function weightf

end module evolve
