!>
!! \brief This module contains routines for calculating the ionization
!and temperature evolution of the entire grid (3D).
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
  ! doric : ionization calculation for one point + photo-ionization
  ! rates
  ! tped : temperature,pressure,electron density calculation

  use precision, only: dp
  use my_mpi ! supplies all the MPI and OpenMP definitions and variables
  use file_admin, only: logf, timefile, iterdump, results_dir, dump_dir
  use clocks, only: timestamp_wallclock
  use c2ray_parameters, only: convergence_fraction, use_LLS
  use sizes, only: Ndim, mesh

  use material, only: ndens
  use material, only: xh_hot
  use material, only: xh,xhe
  use material, only: isothermal
  use material, only: set_final_temperature_point
  use material, only: protect_ionization_fractions
  use material, only: special
  use sourceprops, only: NumSrc
#ifdef PL
  use sourceprops, only: NumPLSrc
#endif  
#ifdef QUASARS
  use sourceprops, only: NumQSrc
#endif  
  use photonstatistics, only: photon_loss, LLS_loss
  use photonstatistics, only: state_before
  use photonstatistics, only: calculate_photon_statistics
  use photonstatistics, only: report_photonstatistics
  use photonstatistics, only: update_grandtotal_photonstatistics
  use radiation_sizes, only: NumFreqBnd
 
  use evolve_data, only: bb_phih_grid, bb_phihe_grid, bb_phiheat_grid
#ifdef QUASARS
  use evolve_data, only: qpl_phih_grid, qpl_phihe_grid, qpl_phiheat_grid
#endif
#ifdef PL
  use evolve_data, only: pl_phih_grid, pl_phihe_grid, pl_phiheat_grid
#endif
  use evolve_data, only: xh_av, xhe_av, xh_intermed, xhe_intermed
  use evolve_data, only: xh_hot_av, xh_hot_intermed
  use evolve_data, only: photon_loss_all
#ifdef MPI
  use evolve_data, only: buffer
#endif
  use evolve_point, only: local_chemistry
  use evolve_source, only: sum_nbox,sum_nbox_all

  use evolve_point, only: evolve0d_global
  use master_slave_processing, only: do_grid

  implicit none

  save

  private

  ! Minimum number of cells which are allowed to be non-converged
  integer :: conv_criterion 
  ! Flag variable (passed back from evolve0D_global)
  integer :: conv_flag

  ! Wall clock counting
  ! 8 bytes to beat the maxcount
  integer(kind=8) :: wallclock1
  integer(kind=8) :: wallclock2
  integer(kind=8) :: countspersec


  public :: evolve3D

#ifdef MPI
    integer :: mympierror
#endif

contains

  ! ===========================================================================

  !> Evolve the entire grid over a time step dt
  subroutine evolve3D (time,dt,iter_restart)

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
    ! 21-Aug-2006 (GM) : MPI parallelization over the sources (static
    ! model).
    ! 28-Feb-2008 (GM) : Added master-slave model for distributing
    !                    over the processes. The program decides which
    !                    model to use.

    ! The time step
    real(kind=dp),intent(in) :: time !< time 
    real(kind=dp),intent(in) :: dt !< time step
    integer,intent(in) :: iter_restart !< iteration restart flag

    ! Loop variables
    integer :: niter  ! iteration counter
    integer :: i,j,k

    ! Phase which is being processed
    character(len=1) :: phase_type
    
    ! End of declarations

    ! Initialize wall clock counter (for dumps)
    call system_clock(wallclock1)

    ! Initial state (for photon statistics)
    call state_before (xh,xhe)

    ! Set the conv_criterion, if there are few sources we should make
    ! sure that things are converged around these sources.

   ! conv_criterion=min(int(convergence_fraction*mesh(1)*mesh(2)*mesh(3)),(NumSrc-1)/3
   ! )
    conv_criterion=min(int(convergence_fraction*mesh(1)*mesh(2)*mesh(3)),NumSrc)   
    
    ! Report time
    if (rank == 0) write(timefile,"(A,F8.1)") &
         "Time before starting iteration: ", timestamp_wallclock ()

    !if (iter_restart /= 0)
    !check_dump_for_sourcetype(iter_restart,niter,source_type)
    ! 
    ! Iterate to reach convergence for multiple sources
    ! Process different groups of sources separately
    ! initialize average and intermediate results to initial values
    if (iter_restart == 0) then
       do k=1,mesh(3)
          do j=1,mesh(2)
             do i=1,mesh(1)
                if (special(i,j,k)) then
                   xh_av(i,j,k,0)=(1.0_dp-xh_hot(i,j,k))*xh(i,j,k,0)
                   ! We assume that in the partially ionized cells, the
                   ! He2+ fraction is the same over the entire cell.
                   ! The He+ fraction is 1 in the hot part and taken
                   ! from the cold value in the cold part.
                   xhe_av(i,j,k,1)=(1.0-xhe(i,j,k,2))*(xh_hot(i,j,k)+ &
                        (1.0-xh_hot(i,j,k))*xhe(i,j,k,1))
                else
                   xh_av(i,j,k,0)=xh(i,j,k,0)
                   xhe_av(i,j,k,1)=xhe(i,j,k,1)
                endif
             enddo
          enddo
       enddo
       xh_hot_av(:,:,:)=xh_hot(:,:,:)
       xh_hot_intermed(:,:,:)=xh_hot(:,:,:)
       xh_av(:,:,:,1)=xh_hot(:,:,:)*xh(:,:,:,1)
       xh_intermed(:,:,:,:)=xh_av(:,:,:,:)
       ! The He2+ fraction is assumed to be the same in the hot and cold
       ! phaes
       xhe_av(:,:,:,2)=xhe(:,:,:,2)
       xh_av(:,:,:,0)=1.0_dp-(xhe_av(:,:,:,1)+xhe_av(:,:,:,2))
       xhe_intermed(:,:,:,:)=xhe_av(:,:,:,:)
       niter=0 ! iteration starts at zero
       conv_flag=mesh(1)*mesh(2)*mesh(3) ! initialize non-convergence
       phase_type="H"
    else
       ! Reload xh_av,xh_intermed,photon_loss,niter
       call start_from_dump(iter_restart,niter,phase_type)
       call global_pass (conv_flag,dt,phase_type)
    endif


    if (phase_type == "H") then
       call process_phase(niter,dt,"H")
       niter=0 ! iteration starts at zero
       conv_flag=mesh(1)*mesh(2)*mesh(3) ! initialize non-convergence 
       call process_phase(niter,dt,"C")
    else
       call process_phase(niter,dt,"C")
    endif

    ! Calculate photon statistics
    call calculate_photon_statistics (dt,xh,xh_av,xhe,xhe_av) 
    call report_photonstatistics (dt)
    call update_grandtotal_photonstatistics (dt)

  end subroutine evolve3D

  ! ===========================================================================

  subroutine process_phase(niter,dt,phase_type)

    integer,intent(inout) :: niter  ! iteration counter
    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    character(len=1),intent(in) :: phase_type !< type of phase ("hot" or"cold")
    
    ! Iterate to reach convergence for multiple sources
    do
       ! Update xh if converged and exit
       ! This should be < and NOT <= for the case of few sources:
       ! n sources will affect at least n cells on the first pass.
       ! We need to give them another iteration to allow them to
       ! work together.
       ! GM/130819: We additionally force to do at least two iterations
       ! by
       ! testing for niter.
       
       if (conv_flag < conv_criterion .and. niter > 1) then
          if (phase_type == "H" ) then
             xh_hot(:,:,:)=xh_hot_intermed(:,:,:)
          else
             xh(:,:,:,:)=xh_intermed(:,:,:,:)
             xhe(:,:,:,:)=xhe_intermed(:,:,:,:)
             call set_final_temperature_point ()
          endif

          ! Report
          if (rank == 0) then
             write(logf,"(2A)") "Multiple sources convergence reached for phase ", &
                  phase_type
             write(logf,*) "Test 1 values: ",conv_flag, conv_criterion
             !write(logf,*) "Test 2 values: ",rel_change_sum_xh, &
             !     convergence_fraction
          endif
          exit
       else
          if (niter > 500) then
             ! Complain about slow convergence
             if (rank == 0) write(logf,*) 'Multiple sources not converging'
             exit
          endif
       endif
       
       ! Iteration loop counter
       niter=niter+1
       
       ! Set all photo-ionization rates to zero
       call set_rates_to_zero
       
       ! Pass over all sources
       if (NumSrc > 0) then
          call pass_all_sources (niter,dt,phase_type)
          
          ! Report subbox statistics
          if (rank == 0) &
               write(logf,*) "Average number of subboxes: ", &
               real(sum_nbox_all)/real(NumSrc)
          
          if (rank == 0) then
             call system_clock(wallclock2,countspersec)
             ! Write iteration dump if more than 15 minutes have passed.
             ! system_clock starts counting at 0 when it reaches
             ! a max value. To catch this, test also for negative
             ! values of wallclock2-wallclock1
             write(logf,*) "Time and limit are: ", &
                  wallclock2-wallclock1, 15.0*60.0*countspersec
             if (wallclock2-wallclock1 > 15*60*countspersec .or. &
                  wallclock2-wallclock1 < 0 ) then
                call write_iteration_dump(niter,phase_type)
                wallclock1=wallclock2
             endif
          endif
          
       endif
       
       ! Do a global pass applying all the rates
       call global_pass (conv_flag,dt,phase_type)
       
       ! Report time
       if (rank == 0) write(timefile,"(A,I3,A,F8.1)") &
            "Time after iteration ",niter," : ", timestamp_wallclock ()
       
    enddo
    
  end subroutine process_phase

  ! ===========================================================================

  subroutine write_iteration_dump (niter,phase_type)

    use material, only:temperature_grid

    integer,intent(in) :: niter  ! iteration counter
    character(len=1),intent(in) :: phase_type

    integer :: ndump=0
    
    character(len=20) :: iterfile

    ! Report time
    write(timefile,"(A,F8.1)") &
         "Time before writing iterdump: ", timestamp_wallclock ()

    ndump=ndump+1
    if (mod(ndump,2) == 0) then
       iterfile="iterdump2.bin"
    else
       iterfile="iterdump1.bin"
    endif

    open(unit=iterdump,file=trim(adjustl(dump_dir))//iterfile,form="unformatted",&
    status="unknown")

    write(iterdump) niter
    write(iterdump) phase_type
    write(iterdump) photon_loss_all
    write(iterdump) bb_phih_grid
#ifdef QUASARS
    write(iterdump) qpl_phih_grid
#endif
#ifdef PL         
    write(iterdump) pl_phih_grid
#endif
    write(iterdump) xh_hot_av
    write(iterdump) xh_hot_intermed
    write(iterdump) xh_av
    write(iterdump) xh_intermed
    write(iterdump) bb_phihe_grid
#ifdef QUASARS
    write(iterdump) qpl_phihe_grid
#endif
#ifdef PL         
    write(iterdump) pl_phihe_grid
#endif
    write(iterdump) xhe_av
    write(iterdump) xhe_intermed
    if (.not.isothermal) then
       write(iterdump) bb_phiheat_grid
#ifdef QUASARS
       write(iterdump) qpl_phiheat_grid
#endif
#ifdef PL
       write(iterdump) pl_phiheat_grid
#endif
       write(iterdump) temperature_grid
    endif
    close(iterdump)

    ! Report time
    write(timefile,"(A,F8.1)") &
         "Time after writing iterdump: ", timestamp_wallclock ()

  end subroutine write_iteration_dump

  ! ===========================================================================

  subroutine start_from_dump(iter_restart,niter,phase_type)

    use material, only: temperature_grid

    integer,intent(in) :: iter_restart  ! restart flag
    integer,intent(out) :: niter  ! iteration counter
    character(len=1),intent(out) :: phase_type

    character(len=20) :: iterfile

    if (iter_restart == 0) then
       if (rank == 0) &
            write(logf,*) "Warning: start_from_dump called incorrectly"
    else
       if (rank == 0) then

          ! Report time
          write(timefile,"(A,F8.1)") &
               "Time before reading iterdump: ", timestamp_wallclock ()

          ! Set file to read (depending on restart flag)
          select case (iter_restart)
          case (1) 
             iterfile="iterdump1.bin"
          case (2) 
             iterfile="iterdump2.bin"
          case (3) 
             iterfile="iterdump.bin"
          end select

          open(unit=iterdump,file=trim(adjustl(dump_dir))//iterfile, &
               form="unformatted",status="old")

          read(iterdump) niter
          read(iterdump) phase_type
          read(iterdump) photon_loss_all
          read(iterdump) bb_phih_grid
#ifdef QUASARS
          read(iterdump) qpl_phih_grid
#endif
#ifdef PL         
          read(iterdump) pl_phih_grid
#endif
          read(iterdump) xh_hot_av
          read(iterdump) xh_hot_intermed
          read(iterdump) xh_av
          read(iterdump) xh_intermed
          read(iterdump) bb_phihe_grid
#ifdef QUASARS
          read(iterdump) qpl_phihe_grid
#endif
#ifdef PL         
          read(iterdump) pl_phihe_grid
#endif
          read(iterdump) xhe_av
          read(iterdump) xhe_intermed
          if (.not.isothermal) then
             read(iterdump) bb_phiheat_grid
#ifdef QUASARS
             read(iterdump) qpl_phiheat_grid
#endif
#ifdef PL
             read(iterdump) pl_phiheat_grid
#endif
             read(iterdump) temperature_grid
          endif

          close(iterdump)
          write(logf,*) "Read iteration ",niter," from dump file"
          write(logf,*) 'photon loss counter: ',photon_loss_all
          write(logf,*) "Intermediate result for mean ionization fraction: ", &
               sum(xh_intermed(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3)),&
               sum(xhe_intermed(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3)),&
               sum(xhe_intermed(:,:,:,2))/real(mesh(1)*mesh(2)*mesh(3))
       endif
       
#ifdef MPI       
       ! Distribute the input parameters to the other nodes
       call MPI_BCAST(niter,1, &
            MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(phase_type,1, &
            MPI_CHARACTER,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(photon_loss_all,NumFreqBnd, &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(bb_phih_grid,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#ifdef QUASARS
       call MPI_BCAST(qpl_phih_grid,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif
#ifdef PL
       call MPI_BCAST(pl_phih_grid,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif
       call MPI_BCAST(bb_phihe_grid,mesh(1)*mesh(2)*mesh(3)*2, &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#ifdef QUASARS
       call MPI_BCAST(qpl_phihe_grid,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif
#ifdef PL
       call MPI_BCAST(pl_phihe_grid,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif
       call MPI_BCAST(xh_av,mesh(1)*mesh(2)*mesh(3)*2, &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(xh_hot_av,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(xhe_av,mesh(1)*mesh(2)*mesh(3)*3, &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(xh_intermed,mesh(1)*mesh(2)*mesh(3)*2, &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(xh_hot_intermed,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       call MPI_BCAST(xhe_intermed,mesh(1)*mesh(2)*mesh(3)*3, &
            MPI_DOUBLE_PRECISION,0,&
            MPI_COMM_NEW,mympierror)
       if (.not.isothermal) then
          call MPI_BCAST(bb_phiheat_grid,mesh(1)*mesh(2)*mesh(3), &
               MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#ifdef QUASARS
       call MPI_BCAST(qpl_phiheat_grid,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif
#ifdef PL
       call MPI_BCAST(pl_phiheat_grid,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
#endif
          call MPI_BCAST(temperature_grid,mesh(1)*mesh(2)*mesh(3)*3, &
               MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
       endif
#endif

       ! Report time
       write(timefile,"(A,F8.1)") &
            "Time after reading iterdump: ", timestamp_wallclock ()

    endif

  end subroutine start_from_dump

  ! ===========================================================================
  
  subroutine set_rates_to_zero
    
    ! reset global rates to zero for this iteration
    bb_phih_grid(:,:,:)=0.0
    bb_phihe_grid(:,:,:,:)=0.0    
    bb_phiheat_grid(:,:,:)=0.0
#ifdef QUASARS
    qpl_phih_grid(:,:,:)=0.0
    qpl_phihe_grid(:,:,:,:)=0.0    
    qpl_phiheat_grid(:,:,:)=0.0
#endif
#ifdef PL
    pl_phih_grid(:,:,:)=0.0
    pl_phihe_grid(:,:,:,:)=0.0    
    pl_phiheat_grid(:,:,:)=0.0
#endif
    ! reset photon loss counters
    photon_loss(:)=0.0
    LLS_loss = 0.0 ! make this a NumFreqBnd vector if needed later (GM/101129)

  end subroutine set_rates_to_zero

  ! ===========================================================================

  subroutine pass_all_sources(niter,dt,phase_type)
    
    ! For random permutation of sources:
    use  m_ctrper, only: ctrper

    integer,intent(in) :: niter  ! iteration counter
    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    character(len=1),intent(in) :: phase_type !< type of phase

    if (rank == 0) write(logf,*) 'Doing all sources '

    ! Reset sum of subboxes counter
    sum_nbox=0

    ! Reset local_chemistry flag
    local_chemistry=.false.

    ! Make a randomized list of sources :: call in serial
    ! disabled / GM110512
    !if ( rank == 0 ) call ctrper (SrcSeries(1:NumSrc),1.0)

#ifdef MPI
    ! Distribute the source list to the other nodes
    ! disabled / GM110512
    !call
    !MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
    
    ! Ray trace the whole grid for all sources.
    call do_grid (dt,niter,phase_type)

#ifdef MPI
    ! Collect (sum) the rates from all the MPI processes
    call mpi_accumulate_grid_quantities (phase_type)

    ! Only if the do_chemistry routine was called with local option
    ! where the ionization fractions changed during the pass over
    ! all sources
    ! DISABLED in this version
    !!if (local_chemistry) call mpi_accumulate_ionization_fractions
#else
    photon_loss_all(:)=photon_loss(:)
    sum_nbox_all=sum_nbox
#endif
    
#ifdef MPI
    call MPI_BARRIER(MPI_COMM_NEW,mympierror)
#endif

  end subroutine pass_all_sources

  ! ===========================================================================

  subroutine global_pass (conv_flag,dt,phase_type)

    ! Flag variable (passed back from evolve0D_global)
    integer,intent(out) :: conv_flag
    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D
    character(len=1),intent(in) :: phase_type !< type of phase

    integer :: i,j,k  ! mesh position

    ! Mesh position of the cell being treated
    integer,dimension(Ndim) :: pos

    ! Report photon losses over grid boundary 
    ! if (rank == 0) write(logf,*) 'photon loss counter:
    ! ',photon_loss_all(:) 
    
    ! Turn total photon loss into a mean per cell (used in
    ! evolve0d_global)
    ! GM/110225: Possibly this should be 
    ! photon_loss_all(:)/(real(mesh(1))*real(mesh(2))
    ! Since this is the correct answer for an optically thin volume
    ! Compromise:
    ! photon_loss_all(:)/(real(mesh(1))**3)*sum(xh_av(:,:,:,1))
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
    endif   

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
             call evolve0D_global(dt,pos,conv_flag,phase_type)
          enddo
       enddo
    enddo

    ! Report on convergence and intermediate result
    if (rank == 0) then
       write(logf,*) "Number of non-converged points: ",conv_flag
       write(logf,*) "Intermediate result for mean H ionization fraction: ", &
            sum(xh_intermed(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3))
       write(logf,*) &
            "Intermediate result for mean He(+,++) ionization fraction:", &
            sum(xhe_intermed(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3)), &
            sum(xhe_intermed(:,:,:,2))/real(mesh(1)*mesh(2)*mesh(3))
    endif
    
    ! Report on photon conservation
    call calculate_photon_statistics(dt,xh_intermed,xh_av,xhe_intermed,xhe_av) 
    call report_photonstatistics (dt)

  end subroutine global_pass

  ! ===========================================================================

  subroutine mpi_accumulate_grid_quantities(phase_type)
    
    character(len=1),intent(in) :: phase_type !< type of phase

    real(kind=dp) :: LLS_loss_all
    
#ifdef MPI
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

    ! accumulate (sum) the MPI distributed sum of number of boxes
    call MPI_ALLREDUCE(sum_nbox, sum_nbox_all, 1, &
         MPI_INTEGER, MPI_SUM, MPI_COMM_NEW, mympierror)

    ! accumulate (sum) the MPI distributed phih_grid
    call MPI_ALLREDUCE(bb_phih_grid, buffer, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
    ! Overwrite the processor local values with the accumulated value
    bb_phih_grid(:,:,:)=buffer(:,:,:)

    if (phase_type /= "H") then
       call MPI_ALLREDUCE(bb_phihe_grid(:,:,:,0), buffer,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)    
       ! Overwrite the processor local values with the accumulated value
       bb_phihe_grid(:,:,:,0)=buffer(:,:,:)
       
       call MPI_ALLREDUCE(bb_phihe_grid(:,:,:,1), buffer,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)    
       ! Overwrite the processor local values with the accumulated value
       bb_phihe_grid(:,:,:,1)=buffer(:,:,:)
       
       call MPI_ALLREDUCE(bb_phiheat_grid, buffer,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
       ! Overwrite the processor local values with the accumulated value
       bb_phiheat_grid(:,:,:)=buffer(:,:,:)    
       
#ifdef QUASARS
       ! accumulate (sum) the MPI distributed phih_grid
       call MPI_ALLREDUCE(qpl_phih_grid, buffer,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
       ! Overwrite the processor local values with the accumulated value
       qpl_phih_grid(:,:,:)=buffer(:,:,:)
       
       call MPI_ALLREDUCE(qpl_phihe_grid(:,:,:,0), buffer,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)    
       ! Overwrite the processor local values with the accumulated value
       qpl_phihe_grid(:,:,:,0)=buffer(:,:,:)
       
       call MPI_ALLREDUCE(qpl_phihe_grid(:,:,:,1), buffer,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)    
       ! Overwrite the processor local values with the accumulated value
       qpl_phihe_grid(:,:,:,1)=buffer(:,:,:)
       
       call MPI_ALLREDUCE(qpl_phiheat_grid, buffer,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
       ! Overwrite the processor local values with the accumulated value
       qpl_phiheat_grid(:,:,:)=buffer(:,:,:)    
#endif
       
#ifdef PL
       ! accumulate (sum) the MPI distributed phih_grid
       call MPI_ALLREDUCE(pl_phih_grid, buffer, mesh(1)*mesh(2)*mesh(3),&
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
       ! Overwrite the processor local values with the accumulated value
       pl_phih_grid(:,:,:)=buffer(:,:,:)
       
       call MPI_ALLREDUCE(pl_phihe_grid(:,:,:,0), buffer,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)    
       ! Overwrite the processor local values with the accumulated value
       pl_phihe_grid(:,:,:,0)=buffer(:,:,:)
       
       call MPI_ALLREDUCE(pl_phihe_grid(:,:,:,1), buffer,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)    
       ! Overwrite the processor local values with the accumulated value
       pl_phihe_grid(:,:,:,1)=buffer(:,:,:)
       
       call MPI_ALLREDUCE(pl_phiheat_grid, buffer,mesh(1)*mesh(2)*mesh(3), &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
       ! Overwrite the processor local values with the accumulated value
       pl_phiheat_grid(:,:,:)=buffer(:,:,:)    
#endif
    endif

#endif

  end subroutine mpi_accumulate_grid_quantities

  ! ===========================================================================

  subroutine mpi_accumulate_ionization_fractions

#ifdef MPI
    ! accumulate (max) MPI distributed xh_av
    call MPI_ALLREDUCE(xh_av(:,:,:,1), buffer, mesh(1)*mesh(2)*mesh(3),&
         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_NEW, mympierror)
    ! Overwrite the processor local values with the accumulated value
    xh_av(:,:,:,1) = buffer(:,:,:)
    
    call MPI_ALLREDUCE(xhe_av(:,:,:,0), buffer, mesh(1)*mesh(2)*mesh(3),&
         MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_NEW, mympierror)      
    ! Overwrite the processor local values with the accumulated value
    xhe_av(:,:,:,0) = buffer(:,:,:)
    
    call MPI_ALLREDUCE(xhe_av(:,:,:,2), buffer, mesh(1)*mesh(2)*mesh(3),&
         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_NEW, mympierror)   
    ! Overwrite the processor local values with the accumulated value
    xhe_av(:,:,:,2) = buffer(:,:,:)
    
    !xh_av(:,:,:,0) = max(0.0_dp,min(1.0_dp,1.0-xh_av(:,:,:,1)))
    !xhe_av(:,:,:,1) =
    !max(0.0_dp,min(1.0_dp,1.0-xhe_av(:,:,:,0)-xhe_av(:,:,:,2)))    
    xh_av(:,:,:,0) = 1.0-xh_av(:,:,:,1)
    xhe_av(:,:,:,1) = 1.0-xhe_av(:,:,:,0)-xhe_av(:,:,:,2)
    
    ! Check the hydrogen and helium fractions for unphysical values
    call protect_ionization_fractions(xh_av,0,1,0,"xh_av(0)")
    call protect_ionization_fractions(xhe_av,0,2,0,"xhe_av(0)")
    
    ! accumulate (max) MPI distributed xh_intermed
    call MPI_ALLREDUCE(xh_intermed(:,:,:,1), buffer, &
         mesh(1)*mesh(2)*mesh(3), MPI_DOUBLE_PRECISION, MPI_MAX, &
         MPI_COMM_NEW, mympierror)
    ! Overwrite the processor local values with the accumulated value
    xh_intermed(:,:,:,1)=buffer(:,:,:)
    
    call MPI_ALLREDUCE(xhe_intermed(:,:,:,0), buffer, &
         mesh(1)*mesh(2)*mesh(3), MPI_DOUBLE_PRECISION, MPI_MIN, &
         MPI_COMM_NEW, mympierror)
    ! Overwrite the processor local values with the accumulated value
    xhe_intermed(:,:,:,0)=buffer(:,:,:)
    
    call MPI_ALLREDUCE(xhe_intermed(:,:,:,2), buffer, &
         mesh(1)*mesh(2)*mesh(3), MPI_DOUBLE_PRECISION, MPI_MAX, &
         MPI_COMM_NEW, mympierror)
       ! Overwrite the processor local values with the accumulated value
    xhe_intermed(:,:,:,2)=buffer(:,:,:)
    
    !xh_intermed(:,:,:,0)=max(0.0_dp,min(1.0_dp,1.0-xh_intermed(:,:,:,1)))
    !xhe_intermed(:,:,:,1)=max(0.0_dp,min(1.0_dp,1.0-xhe_intermed(:,:,:,0)-xhe_intermed(:,:,:,2)))
    xh_intermed(:,:,:,0)=1.0-xh_intermed(:,:,:,1)
    xhe_intermed(:,:,:,1)=1.0-xhe_intermed(:,:,:,0)-xhe_intermed(:,:,:,2)
    
    ! Check the hydrogen and helium fractions for unphysical values
    call protect_ionization_fractions(xh_intermed,0,1,0,"xh_intermed(0)")
    call protect_ionization_fractions(xhe_intermed,0,2,0,"xhe_intermed(0)")
#endif
    
  end subroutine mpi_accumulate_ionization_fractions

end module evolve
