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
  use file_admin, only: logf, timefile, iterdump, results_dir, dump_dir
  use clocks, only: timestamp_wallclock
  use c2ray_parameters, only: convergence_fraction, use_LLS
  use sizes, only: Ndim, mesh

  use material, only: ndens
  use material, only: xh,xhe
  use material, only: isothermal
  use material, only: set_final_temperature_point
  use material, only: protect_ionization_fractions
  use sourceprops, only: NumSrc
  use photonstatistics, only: photon_loss, LLS_loss
  use photonstatistics, only: state_before
  use photonstatistics, only: calculate_photon_statistics
  use photonstatistics, only: report_photonstatistics
  use photonstatistics, only: update_grandtotal_photonstatistics
  use radiation_sizes, only: NumFreqBnd
 
  use evolve_data, only: phih_grid, phihe_grid, phiheat
  use evolve_data, only: xh_av, xhe_av, xh_intermed, xhe_intermed
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

  public :: evolve3D

#ifdef MPI
    integer :: mympierror
#endif

contains

  ! ===========================================================================

  !> Evolve the entire grid over a time step dt
  subroutine evolve3D (time,dt,restart)

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
    real(kind=dp),intent(in) :: time !< time 
    real(kind=dp),intent(in) :: dt !< time step
    integer,intent(in) :: restart !< restart flag

    ! Loop variables
    integer :: niter  ! iteration counter

    ! Wall clock counting
    ! 8 bytes to beat the maxcount
    integer(kind=8) :: wallclock1
    integer(kind=8) :: wallclock2
    integer(kind=8) :: countspersec

    ! Flag variable (passed back from evolve0D_global)
    integer :: conv_flag

    ! Minimum number of cells which are allowed to be non-converged
    integer :: conv_criterion 

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
       niter=0 ! iteration starts at zero
       conv_flag=mesh(1)*mesh(2)*mesh(3) ! initialize non-convergence 
    else
       ! Reload xh_av,xh_intermed,photon_loss,niter
       call start_from_dump(restart,niter)
       call global_pass (conv_flag,dt)
    endif

    ! Set the conv_criterion, if there are few sources we should make
    ! sure that things are converged around these sources.

   ! conv_criterion=min(int(convergence_fraction*mesh(1)*mesh(2)*mesh(3)),(NumSrc-1)/3 )
    conv_criterion=min(int(convergence_fraction*mesh(1)*mesh(2)*mesh(3)),NumSrc)   
    
    ! Report time
    if (rank == 0) write(timefile,"(A,F8.1)") &
         "Time before starting iteration: ", timestamp_wallclock ()

    ! Iterate to reach convergence for multiple sources
    do
       ! Update xh if converged and exit
       ! This should be < and NOT <= for the case of few sources:
       ! n sources will affect at least n cells on the first pass.
       ! We need to give them another iteration to allow them to
       ! work together.
       ! GM/130819: We additionally force to do at least two iterations by
       ! testing for niter.

       if (conv_flag < conv_criterion .and. niter > 1) then
          xh(:,:,:,:)=xh_intermed(:,:,:,:)
          xhe(:,:,:,:)=xhe_intermed(:,:,:,:)
          call set_final_temperature_point

          ! Report
          if (rank == 0) then
             write(logf,*) "Multiple sources convergence reached"
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
          call pass_all_sources (niter,dt)

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
                call write_iteration_dump(niter)
                wallclock1=wallclock2
             endif
          endif
          
       endif

       ! Do a global pass applying all the rates
       call global_pass (conv_flag,dt)

       ! Report time
       if (rank == 0) write(timefile,"(A,I3,A,F8.1)") &
            "Time after iteration ",niter," : ", timestamp_wallclock ()
    enddo

    ! Calculate photon statistics
    call calculate_photon_statistics (dt,xh,xh_av,xhe,xhe_av) 
    call report_photonstatistics (dt)
    call update_grandtotal_photonstatistics (dt)

  end subroutine evolve3D

  ! ===========================================================================

  subroutine write_iteration_dump (niter)

    use material, only:temperature_grid

    integer,intent(in) :: niter  ! iteration counter

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

    open(unit=iterdump,file=trim(adjustl(dump_dir))//iterfile,form="unformatted", &
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

    ! Report time
    write(timefile,"(A,F8.1)") &
         "Time after writing iterdump: ", timestamp_wallclock ()

  end subroutine write_iteration_dump

  ! ===========================================================================

  subroutine start_from_dump(restart,niter)

    use material, only: temperature_grid

    integer,intent(in) :: restart  ! restart flag
    integer,intent(out) :: niter  ! iteration counter

    character(len=20) :: iterfile

    if (restart == 0) then
       if (rank == 0) &
            write(logf,*) "Warning: start_from_dump called incorrectly"
    else
       if (rank == 0) then

          ! Report time
          write(timefile,"(A,F8.1)") &
               "Time before reading iterdump: ", timestamp_wallclock ()

          ! Set file to read (depending on restart flag)
          select case (restart)
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
       if (.not.isothermal) then
          call MPI_BCAST(phiheat,mesh(1)*mesh(2)*mesh(3), &
               MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW,mympierror)
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
    phih_grid(:,:,:)=0.0
    phihe_grid(:,:,:,:)=0.0    
    phiheat(:,:,:)=0.0
    ! reset photon loss counters
    photon_loss(:)=0.0
    LLS_loss = 0.0 ! make this a NumFreqBnd vector if needed later (GM/101129)

  end subroutine set_rates_to_zero

  ! ===========================================================================

  subroutine pass_all_sources(niter,dt)
    
    ! For random permutation of sources:
    use  m_ctrper, only: ctrper

    integer,intent(in) :: niter  ! iteration counter
    real(kind=dp),intent(in) :: dt  !< time step, passed on to evolve0D

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
    !call MPI_BCAST(SrcSeries,NumSrc,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif
    
    ! Ray trace the whole grid for all sources.
    call do_grid (dt,niter)

#ifdef MPI
    ! Collect (sum) the rates from all the MPI processes
    call mpi_accumulate_grid_quantities

    ! Only if the do_chemistry routine was called with local option
    ! where the ionization fractions changed during the pass over
    ! all sources
    if (local_chemistry) call mpi_accumulate_ionization_fractions
#else
    photon_loss_all(:)=photon_loss(:)
    sum_nbox_all=sum_nbox
#endif
    
#ifdef MPI
    call MPI_BARRIER(MPI_COMM_NEW,mympierror)
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
             call evolve0D_global(dt,pos,conv_flag)
          enddo
       enddo
    enddo

    ! Report on convergence and intermediate result
    if (rank == 0) then
       write(logf,*) "Number of non-converged points: ",conv_flag
       write(logf,*) "Intermediate result for mean H ionization fraction: ", &
            sum(xh_intermed(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3))
       write(logf,*) &
            "Intermediate result for mean He(+,++) ionization fraction: ", &
            sum(xhe_intermed(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3)), &
            sum(xhe_intermed(:,:,:,2))/real(mesh(1)*mesh(2)*mesh(3))
    endif
    
    ! Report on photon conservation
    call calculate_photon_statistics (dt,xh_intermed,xh_av,xhe_intermed,xhe_av) 
    call report_photonstatistics (dt)

  end subroutine global_pass

  ! ===========================================================================

  subroutine mpi_accumulate_grid_quantities
    
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

    ! accumulate (sum) the MPI distributed phih_grid
    call MPI_ALLREDUCE(phih_grid, buffer, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
    ! Overwrite the processor local values with the accumulated value
    phih_grid(:,:,:)=buffer(:,:,:)

     call MPI_ALLREDUCE(phihe_grid(:,:,:,0), buffer, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)    
    ! Overwrite the processor local values with the accumulated value
    phihe_grid(:,:,:,0)=buffer(:,:,:)

    call MPI_ALLREDUCE(phihe_grid(:,:,:,1), buffer, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)    
    ! Overwrite the processor local values with the accumulated value
    phihe_grid(:,:,:,1)=buffer(:,:,:)
    
    call MPI_ALLREDUCE(phiheat, buffer, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_NEW, mympierror)
    ! Overwrite the processor local values with the accumulated value
    phiheat(:,:,:)=buffer(:,:,:)    
        
    ! accumulate (sum) the MPI distributed sum of number of boxes
    call MPI_ALLREDUCE(sum_nbox, sum_nbox_all, 1, &
         MPI_INTEGER, MPI_SUM, MPI_COMM_NEW, mympierror)
#endif

  end subroutine mpi_accumulate_grid_quantities

  ! ===========================================================================

  subroutine mpi_accumulate_ionization_fractions

#ifdef MPI
    ! accumulate (max) MPI distributed xh_av
    call MPI_ALLREDUCE(xh_av(:,:,:,1), buffer, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_NEW, mympierror)
    ! Overwrite the processor local values with the accumulated value
    xh_av(:,:,:,1) = buffer(:,:,:)
    
    call MPI_ALLREDUCE(xhe_av(:,:,:,0), buffer, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_NEW, mympierror)      
    ! Overwrite the processor local values with the accumulated value
    xhe_av(:,:,:,0) = buffer(:,:,:)
    
    call MPI_ALLREDUCE(xhe_av(:,:,:,2), buffer, mesh(1)*mesh(2)*mesh(3), &
         MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_NEW, mympierror)   
    ! Overwrite the processor local values with the accumulated value
    xhe_av(:,:,:,2) = buffer(:,:,:)
    
    !xh_av(:,:,:,0) = max(0.0_dp,min(1.0_dp,1.0-xh_av(:,:,:,1)))
    !xhe_av(:,:,:,1) = max(0.0_dp,min(1.0_dp,1.0-xhe_av(:,:,:,0)-xhe_av(:,:,:,2)))    
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
