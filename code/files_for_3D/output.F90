!>
!! \brief This module contains routines for file output
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 2013-02-22 (but older)
!!
!! \b Version: 3D, hydrogen + helium

module output_module
  
  ! This file contains routines having to do with the output
  ! of the C2-Ray program.
  
  ! setup_out : open files
  ! close_down : close files
  ! output : write output

  use precision, only: dp
  use my_mpi
  use file_admin, only: stdinput, results_dir, file_input, logf
  use material, only: isothermal
  use sizes, only: mesh
  use grid, only: x, vol
  use material, only: xh, temperature_grid, ndens, xhe
  use evolve, only: phih_grid, phiheat
  use sourceprops, only: srcpos, NormFlux, NormFluxPL, NumSrc
  use photonstatistics, only: do_photonstatistics, total_ion, totrec
  use photonstatistics, only: totcollisions, dh0, dhe0, dhe2, grtotal_ion
  use photonstatistics, only: photon_loss, grtotal_src
  use photonstatistics, only: initialize_photonstatistics
  use radiation, only: T_eff,R_star,L_star,S_star, pl_S_star


  implicit none
  
  private

  ! To controle what kind of output is produced we use an array of flags.
  ! There can be at most max_output_streams types of output.
  integer,parameter :: max_output_streams=5 !< maximum number of output streams
  integer,dimension(max_output_streams) :: streams !< flag array for output streams

  public :: setup_output, output, close_down

contains
  !----------------------------------------------------------------------------

  !> Initializes output streams
  subroutine setup_output ()
    
    ! Sets up output stream
    
    ! Version: Five streams

    ! Stream1:
    ! Ifront1.out contains a line of constant y and z going through the
    ! centre for all timesteps. (formatted)
    
    ! Stream2: 
    ! Ionization fractions for the full data cube (unformatted)
    ! xfrac3d_",f5.3,".bin"
    ! xfrac3dHe1_",f5.3,".bin"
    ! xfrac3dHe2_",f5.3,".bin"
    ! and if non isothermal also the temperature for the full data cube
    ! Temper3d_",f5.3,".bin"

    ! Stream3: 
    ! Ionization rate for the full data cube (unformatted)
    ! IonRates3d_",f5.3,".bin"
    
    ! Stream 4:
    ! Ionization fractions in a plane for one time step
    ! Ifront2_xy_",f5.3,".bin"
    ! Ifront2_xz_",f5.3,".bin"
    ! Ifront2_yz_",f5.3,".bin"
    ! Densities in a plane for one time step
    ! ndens_xy_",f5.3,".bin"
    ! ndens_xz_",f5.3,".bin"
    ! ndens_yz_",f5.3,".bin"
    
    ! Stream 5:
    ! Not used at the moment
    ! photon statistics

    if (rank == 0) then
       if (.not.file_input) then
          write(*,*) "Which output streams do you want?"
          write(*,*) "Enter a mask for streams 1 through 5"
          write(*,*) "E.g. 1,0,1,0,0 means streams 1 and 3, but not 2, 4, 5"
       endif
       read(stdinput,*) streams(1),streams(2),streams(3),streams(4),streams(5)
       
       ! Open files
       if (do_photonstatistics) then
          ! Fortran 2003 standard
          open(unit=90,file=trim(adjustl(results_dir))//"PhotonCounts.out", &
               form="formatted",status="unknown",position="append")
          write(90,*) "This file is supposed to contain numbers related ", &
               "to photon conservation. ", & 
               "However, this is still being developed. "
          !write(90,*) "redshift, total ionizations, total photons, " &
          !     "photon conservation number, ", &
          !     "fraction new ionization, fraction recombinations, ", &
          !     "fraction photon losses, fraction collisional ionization, ", &
          !     "grand total photon conservation number"

          ! Fortran 2003 standard
          open(unit=95,file=trim(adjustl(results_dir))//"PhotonCounts2.out", &
               form="formatted",status="unknown",position="append")
          write(95,*) "Columns: redshift, total number of ions, ", &
               "grand total ionizing photons, ", &
               "mean ionization fraction (H2,He2,He3) ", &
               "(by volume and then mass)"
       endif
#ifdef MPILOG
       write(logf,*) "Making output streams according to: ", &
            streams(1),streams(2),streams(3),streams(4),streams(5)
#endif
    endif
    if (do_photonstatistics) call initialize_photonstatistics ()

#ifdef MPILOG
    write(logf,*) "End of setup output"
#endif

  end subroutine setup_output
  
  !-----------------------------------------------------------------------------

  !> Closes global output files which have been open the entire run

  subroutine close_down ()
    
    ! Closes down
    
    if (rank == 0) then
       close(unit=90)
       close(unit=95)
    endif

  end subroutine close_down
  
  !----------------------------------------------------------------------------

  !> Produce output for a time frame
  subroutine output(zred_now,time,dt,photcons_flag)

    ! Simple output routine.

    real(kind=dp),intent(in) :: zred_now !< current redshift
    real(kind=dp),intent(in) :: time !< current simulation time
    real(kind=dp),intent(in) :: dt !< time step taken
    integer,intent(out) :: photcons_flag

    ! Set photon conservation flag to zero on all processors
    photcons_flag=0

#ifdef MPILOG     
     write(logf,*) 'output 1'
#endif 
    if (streams(1) == 1) call write_stream1 (zred_now)
#ifdef MPILOG     
     write(logf,*) 'output 2'
#endif 
    if (streams(2) == 1) call write_stream2 (zred_now)
#ifdef MPILOG     
     write(logf,*) 'output 3'
#endif 
    if (streams(3) == 1) call write_stream3 (zred_now)
#ifdef MPILOG     
     write(logf,*) 'output 4'
#endif 
    if (streams(4) == 1) call write_stream4 (zred_now)
#ifdef MPILOG     
     write(logf,*) 'output 5'
#endif 
    if (streams(5) == 1) call write_stream5 (zred_now)

#ifdef MPILOG     
     write(logf,*) 'output 6'
#endif 
    call write_photonstatistics (zred_now,time,dt,photcons_flag)

#ifdef MPILOG     
     write(logf,*) 'output 7'
#endif 

  end subroutine output

  !----------------------------------------------------------------------------

  !> produces output for a time frame. See below for format
  subroutine write_stream1 (zred_now)

    real(kind=dp),intent(in) :: zred_now !< current redshift

    character(len=*),parameter :: basename="Ifront1_"
    character(len=*),parameter :: base_extension=".dat"
    character(len=80) :: file1
    integer :: i,j,k

    if (rank == 0) then
       ! Stream 1
       if (streams(1) == 1) then
          ! Construct file name
          write(file1,"(f6.3)") zred_now
          !file1=trim(adjustl(results_dir))//"Ifront1_"//trim(adjustl(file1))//".dat"
          file1=trim(adjustl(results_dir))//basename//trim(adjustl(file1)) &
               //base_extension

          ! Open file
          open(unit=51,file=file1,form="formatted",status="unknown")

          ! Write data
          if (.not.isothermal) then
             do i=1,mesh(1)
                write(51,"(7(es10.3,1x))")  &
                     xh(i,srcpos(2,1),srcpos(3,1),0), &
                     xh(i,srcpos(2,1),srcpos(3,1),1), &
                     ndens(i,srcpos(2,1),srcpos(3,1)), &
                     xhe(i,srcpos(2,1),srcpos(3,1),0), &
                     xhe(i,srcpos(2,1),srcpos(3,1),1), &                 
                     xhe(i,srcpos(2,1),srcpos(3,1),2), &                        
                     temperature_grid(i,srcpos(2,1),srcpos(3,1),0)
             enddo
          else
             do i=1,mesh(1)
                write(51,"(7(es10.3,1x))")  &
                     xh(i,srcpos(2,1),srcpos(3,1),0), &
                     xh(i,srcpos(2,1),srcpos(3,1),1), &
                     ndens(i,srcpos(2,1),srcpos(3,1)), &
                     xhe(i,srcpos(2,1),srcpos(3,1),0), &
                     xhe(i,srcpos(2,1),srcpos(3,1),1), &                 
                     xhe(i,srcpos(2,1),srcpos(3,1),2)
             enddo
          endif
          ! Close file
          close(51)
       else
          ! Report error
          write(logf,*) "Calling stream 1 output where we should not."
       endif

    endif
  end subroutine write_stream1

  !----------------------------------------------------------------------------

  !> produces output for a time frame. See below for format
  subroutine write_stream2 (zred_now)

    real(kind=dp),intent(in) :: zred_now !< current redshift

    !character(len=*),parameter :: basename="Ifront1_"
    character(len=*),parameter :: base_extension=".bin"
    character(len=80) :: file1
    integer :: i,j,k

    ! Stream 2
    if (rank == 0) then

       if (streams(2) == 1) then

          ! Construct file name
          write(file1,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))//"xfrac3d_"//trim(adjustl(file1))// &
               base_extension

          ! Open, write and close
          open(unit=52,file=file1,form="unformatted",status="unknown")
          write(52) mesh(1),mesh(2),mesh(3)
          write(52) (((xh(i,j,k,1),i=1,mesh(1)),j=1,mesh(2)), &
               k=1,mesh(3))
          close(52)

          ! Construct file name
          write(file1,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))//"xfrac3dHe1_"// &
               trim(adjustl(file1))//base_extension

          ! Open, write and close
          open(unit=62,file=file1,form="unformatted",status="unknown")
          write(62) mesh(1),mesh(2),mesh(3)
          write(62) (((xhe(i,j,k,1),i=1,mesh(1)),j=1,mesh(2)), &
               k=1,mesh(3))
          close(62)

          ! Construct file name
          write(file1,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))//"xfrac3dHe2_"// &
               trim(adjustl(file1))//base_extension

          ! Open, write and close
          open(unit=72,file=file1,form="unformatted",status="unknown")
          write(72) mesh(1),mesh(2),mesh(3)
          write(72) (((xhe(i,j,k,2),i=1,mesh(1)),j=1,mesh(2)), &
               k=1,mesh(3))
          close(72)

       else
          ! Report error
          write(logf,*) "Calling stream 2 output where we should not."
       endif
       
    endif
    
  end subroutine write_stream2

  !----------------------------------------------------------------------------

  !> produces output for a time frame. See below for format
  subroutine write_stream3 (zred_now)

    real(kind=dp),intent(in) :: zred_now !< current redshift

    !character(len=*),parameter :: basename="Ifront1_"
    character(len=*),parameter :: base_extension=".bin"
    character(len=80) :: file1
    integer :: i,j,k

    ! Stream 2
    if (rank == 0) then

       if (streams(3) == 1) then

          if (.not. isothermal) then

             write(file1,"(f6.3)") zred_now
             file1=trim(adjustl(results_dir))//"Temper3D_"// &
                  trim(adjustl(file1))//base_extension

             open(unit=153,file=file1,form="unformatted",status="unknown")
             write(153) mesh(1),mesh(2),mesh(3)
             write(153) (((real(temperature_grid(i,j,k,0)),i=1,mesh(1)), &
                  j=1,mesh(2)),k=1,mesh(3))
             close(153)

#ifdef MPILOG     
             write(logf,*) 'output 3: temper3d'
             flush(logf)
#endif 
          endif

#ifdef MPILOG
          write(logf,*) allocated(phih_grid)
          write(logf,*) 'shape phih_grid: ',shape(phih_grid)
          flush(logf)
#endif 
          write(file1,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))//"IonRates3D_"// &
               trim(adjustl(file1))//base_extension

          open(unit=53,file=file1,form="unformatted",status="unknown")
          write(53) mesh(1),mesh(2),mesh(3)
          write(53) (((real(phih_grid(i,j,k)),i=1,mesh(1)),j=1,mesh(2)), &
               k=1,mesh(3))
          close(53)

          write(file1,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))//"HeatRates3D_"// &
               trim(adjustl(file1))//base_extension

          open(unit=53,file=file1,form="unformatted",status="unknown")
          write(53) mesh(1),mesh(2),mesh(3)
          write(53) (((real(phiheat(i,j,k)),i=1,mesh(1)),j=1,mesh(2)), &
               k=1,mesh(3))
          close(53)

#ifdef MPILOG     
          write(logf,*) 'output 3: IonRates3D'
          flush(logf)
#endif 
       else
          ! Report error
          write(logf,*) "Calling stream 3 output where we should not."
       endif

    endif

  end subroutine write_stream3

  !----------------------------------------------------------------------------

  !> produces output for a time frame. See below for format
  subroutine write_stream4 (zred_now)

    real(kind=dp),intent(in) :: zred_now !< current redshift

    !character(len=*),parameter :: basename="Ifront1_"
    character(len=*),parameter :: base_extension=".bin"
    character(len=80) :: file1,file2,file3
    character(len=6) :: zred_str
    integer :: i,j,k

    ! Stream 2
    if (rank == 0) then

       ! Stream 4
       if (streams(4).eq.1) then

          write(zred_str,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))//"Ifront2d_xy_"// &
               trim(adjustl(zred_str))//".bin"
          file2=trim(adjustl(results_dir))//"Ifront2d_xz_"// &
               trim(adjustl(zred_str))//".bin"
          file3=trim(adjustl(results_dir))//"Ifront2d_yz_"// &
               trim(adjustl(zred_str))//".bin"


          ! xy cut through source 
          open(unit=54,file=file1,form="unformatted",status="unknown")
          write(54) mesh(1),mesh(2)
          write(54) ((real(xh(i,j,mesh(3)/2,1)),i=1,mesh(1)), &
               j=1,mesh(2))
          close(54)

          ! xz cut through source 
          open(unit=55,file=file2,form="unformatted",status="unknown")
          write(55) mesh(1),mesh(3)
          write(55) ((real(xh(i,mesh(2)/2,k,1)),i=1,mesh(1)), &
               k=1,mesh(3))
          close(55)

          ! yz cut through source 
          open(unit=56,file=file3,form="unformatted",status="unknown")
          write(56) mesh(2),mesh(3)
          write(56) ((real(xh(mesh(1)/2,j,k,1)),j=1,mesh(2)), &
               k=1,mesh(3))
          close(56)
       else
          ! Report error
          write(logf,*) "Calling stream 4 output where we should not."
       endif
       
    endif

  end subroutine write_stream4

  !----------------------------------------------------------------------------

  !> produces output for a time frame. See below for format
  subroutine write_stream5 (zred_now)

    real(kind=dp),intent(in) :: zred_now !< current redshift

    !character(len=*),parameter :: basename="Ifront1_"
    character(len=*),parameter :: base_extension=".bin"
    character(len=80) :: file1,file2,file3
    character(len=6) :: zred_str
    integer :: i,j,k

    if (rank == 0) then

       ! Stream 5
       if (streams(5) == 1) then
          write(zred_str,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))//"ndens_xy_"//trim(adjustl(zred_str))//".bin"
          file2=trim(adjustl(results_dir))//"ndens_xz_"//trim(adjustl(zred_str))//".bin"
          file3=trim(adjustl(results_dir))//"ndens_yz_"//trim(adjustl(zred_str))//".bin"

          ! xy cut through source 
          open(unit=57,file=file1,form="unformatted",status="unknown")
          write(57) mesh(1),mesh(2)
          write(57) ((real(ndens(i,j,mesh(3)/2)),i=1,mesh(1)),j=1,mesh(2))
          close(57)

          ! xz cut through source 
          open(unit=58,file=file2,form="unformatted",status="unknown")
          write(58) mesh(1),mesh(3)
          write(58) ((real(ndens(i,mesh(2)/2,k)),i=1,mesh(1)),k=1,mesh(3))
          close(58)

          ! yz cut through source 
          open(unit=59,file=file3,form="unformatted",status="unknown")
          write(59) mesh(2),mesh(3)
          write(59) ((real(ndens(mesh(1)/2,j,k)),j=1,mesh(2)),k=1,mesh(3))
          close(59)
       else
          ! Report error
          write(logf,*) "Calling stream 5 output where we should not."
       endif
       
    endif

  end subroutine write_stream5
      
  !----------------------------------------------------------------------------

  !> produces output for a time frame. See below for format
  subroutine write_photonstatistics (zred_now,time,dt,photcons_flag)

    real(kind=dp),intent(in) :: zred_now !< current redshift
    real(kind=dp),intent(in) :: time !< current simulation time
    real(kind=dp),intent(in) :: dt !< time step taken
    integer,intent(out) :: photcons_flag

    real(kind=dp) :: totalsrc,photcons,total_photon_loss
    real(kind=dp) :: totions,totphots,volfrac(0:2),massfrac(0:2)

#ifdef MPI
    integer :: mympierror
#endif

    if (rank == 0) then
       ! Check if we are tracking photon conservation
       if (do_photonstatistics) then
          ! Photon Statistics
          ! total_ion is the total number of new ionization, plus
          ! the total number of recombinations, and here is also
          ! added the number of photons lost from the grid. Since
          ! this number was divided by the number of cells, we
          ! multiply by this again.
          total_photon_loss=sum(photon_loss)*dt* &
               real(mesh(1))*real(mesh(2))*real(mesh(3))
          !total_ion=total_ion + total_photon_loss
          totalsrc=(sum(NormFlux(1:NumSrc))*S_star + &
               sum(NormFluxPL(1:NumSrc))*pl_S_star)*dt
          photcons=(total_ion-totcollisions)/totalsrc
          !PhotonCounts: time
          !              Number of (ionizations + recombinations) / photons 
          !                   during time step
          !              Number of ionizations /(ionizations + recombinations)
          !              Number of recombinations /(ionizations + recombinations)
          !              Number of (ionizations + recombinations) / photons 
          !              Number of (ionizations + recombinations) / photons 
          !                   since t=0
          if (time > 0.0) then
             !write(90,"(f6.3,8(es10.3))") &
             !     zred_now, &
             !     total_ion, totalsrc, &
             !     photcons, &
             !     (dh0+dhe0+dhe2)/total_ion, &
             !     totrec/total_ion, &
             !     total_photon_loss/totalsrc, &
             !     totcollisions/total_ion, &
             !     grtotal_ion/grtotal_src
          endif
          totions=sum(ndens(:,:,:)*(xhe(:,:,:,2)*2.0_dp+xhe(:,:,:,1)+xh(:,:,:,1)))*vol
          volfrac(0)=sum(xh(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3))
          volfrac(1)=sum(xhe(:,:,:,1))/real(mesh(1)*mesh(2)*mesh(3))
          volfrac(2)=sum(xhe(:,:,:,2))/real(mesh(1)*mesh(2)*mesh(3))
          massfrac(0)=sum(xh(:,:,:,1)*ndens(:,:,:) ) /sum(real(ndens,dp))
          massfrac(1)=sum(xhe(:,:,:,1)*ndens(:,:,:) ) /sum(real(ndens,dp))
          massfrac(2)=sum(xhe(:,:,:,2)*ndens(:,:,:) ) /sum(real(ndens,dp))
          write(95,"(f6.3,8(es10.3))") zred_now,totions,grtotal_src, &
               volfrac,massfrac

!*** for the moment, I turn that off, until I checked, how I calculate those quantities.
          photcons_flag=0
          !if (abs(1.0-photcons) > 0.15) then
             !if ((1.0-photcons) > 0.15 .and. &
              !    total_photon_loss/totalsrc < (1.0-photcons) ) then
              !  photcons_flag=1
                ! Report photon conservation
              !  write(logf,"(A,2(es10.3,x))") &
                  !   "Photon conservation problem: ", &
                    ! photcons, total_photon_loss/totalsrc

             !endif
          !endif
!***
       endif
    endif
    
#ifdef MPI
    call MPI_BCAST(photcons_flag,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

  end subroutine write_photonstatistics

end module output_module
