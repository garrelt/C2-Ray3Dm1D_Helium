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
  use parameter, only: mesh, do_LTE, do_constant
  use array, only: temperature_array
  use array, only: intermediate_end_temperature
  use array, only: intermediate_average_temperature
  use array, only: xHI_array,xHII_array
  use array, only: xHeI_array,xHeII_array,xHeIII_array
  use array, only: photoionization_HI_array
  use array, only: heating_array
  use sourceprops, only: srcpos
  use array, only: number_density_array
  use grid, only: x_coordinate
  use photonstatistics, only: grand_total_ionization_photons, &
                              grand_total_recombination_photons, &
                              grand_total_collisional_photons, &
                              grand_total_recombation_ionization_He, &
                              grand_total_source_photons, &
                              grand_total_escaped_photons
  ! photon statistics
  use parameter, only: do_photonstatistics, total_simulation_time, &
                       time_partition, number_timesteps

  implicit none
  
  private

  ! To controle what kind of output is produced we use an array of flags.
  ! There can be at most max_output_streams types of output.
  integer,parameter :: max_output_streams=5 !< maximum number of output streams
  integer,dimension(max_output_streams) :: streams !< flag array for output streams

  public :: setup_output, output, close_down

  ! applicable to both Test 1 and Test 4
  logical, parameter :: do_test4 = .true.
  !logical, parameter :: do_test4 = .false.

  
  ! applicable to both Test 1 
  !logical, parameter :: do_bb = .true.
  logical, parameter :: do_bb = .false.

  

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initializes output streams
  subroutine setup_output()
    
    implicit none

    if (rank .eq. 0) then
      write(logf,*) "Beginning of output setup" 
    endif

    streams(1) = 1
    streams(2) = 1
    streams(3) = 1
    streams(4) = 1
    streams(5) = 1

    ! Open PhotoCounts files
    if (do_photonstatistics) then

      open(unit=96,file=trim(adjustl(results_dir))//"myphoton.out", &
           form="formatted",status="unknown")
      write(96,*) " grand_total_ionization_photons ", &
                  " grand_total_recombination_photons ", &
                  " grand_total_collisional_photons ", &
                  " grand_total_recombation_ionization_He ", &
                  " grand_total_source_photons ", &
                  " grand_total_escaped_photons"

    endif

    if (rank .eq. 0) then
      write(logf,*) "stream(1) = ",streams(1)
      write(logf,*) "stream(2) = ",streams(2)
      write(logf,*) "stream(3) = ",streams(3)
      write(logf,*) "stream(4) = ",streams(4)
      write(logf,*) "stream(5) = ",streams(5)
      write(logf,*) "End of output setup"
      write(logf,*)      
      flush(logf) 
    endif

  end subroutine setup_output
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Closes global output files which have been open the entire run

  subroutine close_down ()
    
    ! Closes down
    
    if (rank == 0) then
       close(unit=90)
       close(unit=95)
    endif

  end subroutine close_down
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output(zred_now,sim_time,dt,photcons_flag)

    implicit none

    real(kind=dp),intent(in) :: zred_now !< current redshift
    real(kind=dp),intent(in) :: sim_time !< current simulation time
    real(kind=dp),intent(in) :: dt !< time step taken
    integer,intent(out) :: photcons_flag

    ! Set photon conservation flag to zero on all processors
    photcons_flag=0

#ifdef MPILOG     
     write(logf,*) 'output 1'
#endif 
    if (streams(1) == 1) call write_stream1 (zred_now,sim_time)

#ifdef MPILOG     
     write(logf,*) 'output 2'
#endif 
    if (streams(2) == 1) call write_stream2 (zred_now,sim_time)

#ifdef MPILOG     
     write(logf,*) 'output 3'
#endif 
    if (streams(3) == 1) call write_stream3 (zred_now,sim_time)

#ifdef MPILOG     
     write(logf,*) 'output 4'
#endif 
    if (streams(4) == 1) call write_stream4 (zred_now,sim_time)

#ifdef MPILOG     
     write(logf,*) 'output 5'
#endif 
    if (streams(5) == 1) call write_stream5 (zred_now,sim_time)


#ifdef MPILOG     
     write(logf,*) 'output 6'
#endif 
    call write_photonstatistics (zred_now,sim_time,dt,photcons_flag)

#ifdef MPILOG     
     write(logf,*) 'output 7'
#endif 

  end subroutine output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> produces output for a time frame. See below for format
  subroutine write_stream1 (zred_now,sim_time)

    implicit none

    real(kind=dp),intent(in) :: zred_now !< current redshift
    real(kind=dp),intent(in) :: sim_time !< current simulation time
    character(len=*),parameter :: basename="Ifront1_"
    character(len=*),parameter :: base_extension=".dat"
    character(len=80) :: file1
    integer :: i,j,k

    if (rank == 0) then
       ! Stream 1
       if (streams(1) == 1) then
          ! Construct file name
         if (do_test4) then
          write(file1,"(f6.3)") sim_time/total_simulation_time !Test 4
         else
          write(file1,"(f6.3)") zred_now
         endif
          file1=trim(adjustl(results_dir))//basename//trim(adjustl(file1))//base_extension

          ! Open file
          open(unit=51,file=file1,form="formatted",status="unknown")

          ! Write data
             do i=1,mesh
                write(51,"(7(es10.3,1x))")  &
                     x_coordinate(i), &
                     xHI_array(i,srcpos(2,1),srcpos(3,1)), &
                     xHII_array(i,srcpos(2,1),srcpos(3,1)), &
                     temperature_array(i,srcpos(2,1),srcpos(3,1)), &
                     number_density_array(i,srcpos(2,1),srcpos(3,1)), &
                     xHeI_array(i,srcpos(2,1),srcpos(3,1)), &
                     xHeII_array(i,srcpos(2,1),srcpos(3,1)), &                 
                     xHeIII_array(i,srcpos(2,1),srcpos(3,1))                        

             enddo

          ! Close file
          close(51)
       else
          ! Report error
          write(logf,*) "Calling stream 1 output where we should not."
       endif

    endif
  end subroutine write_stream1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> produces output for a time frame. See below for format
  subroutine write_stream2 (zred_now,sim_time)

    implicit none

    real(kind=dp),intent(in) :: zred_now !< current redshift
    real(kind=dp),intent(in) :: sim_time !< current simulation time
    !character(len=*),parameter :: basename="Ifront1_"
    character(len=*),parameter :: base_extension=".bin"
    character(len=80) :: file1, file_timestep, file_function, file_LTE, file_source,file_mesh
    integer :: i,j,k
    integer :: number_of_ionized_cells

if (do_constant) then
  file_timestep = "constant_"
  write(file_function,"(I10.1)") number_timesteps 
  file_function = "t"//trim(adjustl(file_function))//"_"
else
  file_timestep = "adaptive_"
  write(file_function,"(I10.1)") time_partition
  file_function = "f"//trim(adjustl(file_function))//"_"
endif

if (do_LTE) then
  file_LTE = "do_LTE_"
else
  file_LTE = "no_LTE_"
endif

if (do_bb) then
  file_source = "bb_"
else
  file_source = "pl_"
endif

write(file_mesh,"(I10.1)") mesh

    ! Stream 2
    if (rank == 0) then

       if (streams(2) == 1) then


          ! Construct file name
         if (do_test4) then
          write(file1,"(f6.3)") sim_time/total_simulation_time !Test 4
         else
          write(file1,"(f6.3)") zred_now
         endif
          file1=trim(adjustl(results_dir))//"mesh_"//trim(adjustl(file1))//base_extension
          ! Open, write and close
          open(unit=50,file=file1,form="formatted",status="unknown")
          do i=1,mesh
          write(50,*) x_coordinate(i)
          enddo
          close(50)



          ! Construct file name
         if (do_test4) then
          write(file1,"(f6.3)") sim_time/total_simulation_time !Test 4
          if (do_constant) then
            file1=trim(adjustl(results_dir))//"constant_density_"//trim(adjustl(file1))//base_extension
          else
            file1=trim(adjustl(results_dir))//"adaptive_density_"//trim(adjustl(file1))//base_extension
          endif
         else
          write(file1,"(f6.3)") zred_now
          file1=trim(adjustl(results_dir))//"density_"//trim(adjustl(file1))//base_extension
         endif

          ! Open, write and close
          open(unit=51,file=file1,form="unformatted",status="unknown")
          write(51) mesh,mesh,mesh
          write(51) (((number_density_array(i,j,k),i=1,mesh),j=1,mesh), &
               k=1,mesh)
          close(51)

          ! Construct file name
         if (do_test4) then
          write(file1,"(f6.3)") sim_time/total_simulation_time !Test 4
file1=trim(adjustl(results_dir))//trim(adjustl(file_timestep))//trim(adjustl(file_LTE))//trim(adjustl(file_function))// &
trim(adjustl(file_source))//"mesh_"//trim(adjustl(file_mesh))//"_Temper3d_"//trim(adjustl(file1))//base_extension

          if (do_constant) then

            if (do_LTE) then
              !file1=trim(adjustl(results_dir))//"constant_do_LTE_t10_T1e5_Temper3d_"//trim(adjustl(file1))//base_extension
            else
              !file1=trim(adjustl(results_dir))//"constant_no_LTE_t10_T1e5_Temper3d_"//trim(adjustl(file1))//base_extension
            endif

          else

            if (do_LTE) then
              !file1=trim(adjustl(results_dir))//"adaptive_do_LTE_T1e5_Temper3d_"//trim(adjustl(file1))//base_extension
            else
              !file1=trim(adjustl(results_dir))//"adaptive_no_LTE_T1e5_Temper3d_"//trim(adjustl(file1))//base_extension
            endif

          endif
         else

          write(file1,"(f6.3)") zred_now
file1=trim(adjustl(results_dir))//trim(adjustl(file_timestep))//trim(adjustl(file_LTE))//trim(adjustl(file_function))// &
trim(adjustl(file_source))//"mesh_"//trim(adjustl(file_mesh))//"_Temper3d_"//trim(adjustl(file1))//base_extension

          if (do_constant) then

            if (do_LTE) then      

              if (do_BB) then  
                !file1=trim(adjustl(results_dir))//"constant_do_LTE_bb_Temper3d_"//trim(adjustl(file1))//base_extension
              else
                !file1=trim(adjustl(results_dir))//"constant_do_LTE_pl_Temper3d_"//trim(adjustl(file1))//base_extension
              endif

            else

              if (do_BB) then 
                !file1=trim(adjustl(results_dir))//"constant_no_LTE_bb_Temper3d_"//trim(adjustl(file1))//base_extension
              else
                !file1=trim(adjustl(results_dir))//"constant_no_LTE_pl_Temper3d_"//trim(adjustl(file1))//base_extension
              endif

            endif

          else

            if (do_LTE) then

              if (do_BB) then 
                !file1=trim(adjustl(results_dir))//"adaptive_do_LTE_bb_Temper3d_"//trim(adjustl(file1))//base_extension
              else
                !file1=trim(adjustl(results_dir))//"adaptive_do_LTE_pl_Temper3d_"//trim(adjustl(file1))//base_extension
              endif

            else

              if (do_BB) then 
                !file1=trim(adjustl(results_dir))//"adaptive_no_LTE_bb_Temper3d_"//trim(adjustl(file1))//base_extension
              else
                !file1=trim(adjustl(results_dir))//"adaptive_no_LTE_pl_Temper3d_"//trim(adjustl(file1))//base_extension
              endif

            endif

          endif

         endif

          ! Open, write and close
          open(unit=51,file=file1,form="unformatted",status="unknown")
          write(51) mesh,mesh,mesh
          write(51) (((temperature_array(i,j,k),i=1,mesh),j=1,mesh), &
               k=1,mesh)
          close(51)



          ! Construct file name
         if (do_test4) then
          write(file1,"(f6.3)") sim_time/total_simulation_time !Test 4
file1=trim(adjustl(results_dir))//trim(adjustl(file_timestep))//trim(adjustl(file_LTE))//trim(adjustl(file_function))// &
trim(adjustl(file_source))//"mesh_"//trim(adjustl(file_mesh))//"_xfrac3d_"//trim(adjustl(file1))//base_extension

          if (do_constant) then

            if (do_LTE) then
              !file1=trim(adjustl(results_dir))//"constant_do_LTE_t10_T1e5_xfrac3d_"//trim(adjustl(file1))//base_extension
            else
              !file1=trim(adjustl(results_dir))//"constant_no_LTE_t10_T1e5_xfrac3d_"//trim(adjustl(file1))//base_extension
            endif

          else

            if (do_LTE) then
              !file1=trim(adjustl(results_dir))//"adaptive_do_LTE_T1e5_xfrac3d_"//trim(adjustl(file1))//base_extension
            else
              !file1=trim(adjustl(results_dir))//"adaptive_no_LTE_T1e5_xfrac3d_"//trim(adjustl(file1))//base_extension
            endif
          endif
         else
          write(file1,"(f6.3)") zred_now
file1=trim(adjustl(results_dir))//trim(adjustl(file_timestep))//trim(adjustl(file_LTE))//trim(adjustl(file_function))// &
trim(adjustl(file_source))//"mesh_"//trim(adjustl(file_mesh))//"_xfrac3d_"//trim(adjustl(file1))//base_extension

          if (do_constant) then

            if (do_LTE) then

              if (do_BB) then 
                !file1=trim(adjustl(results_dir))//"constant_do_LTE_bb_xfrac3d_"//trim(adjustl(file1))//base_extension
              else
                !file1=trim(adjustl(results_dir))//"constant_do_LTE_pl_xfrac3d_"//trim(adjustl(file1))//base_extension
              endif

            else

              if (do_BB) then
                !file1=trim(adjustl(results_dir))//"constant_no_LTE_bb_xfrac3d_"//trim(adjustl(file1))//base_extension
              else
                !file1=trim(adjustl(results_dir))//"constant_no_LTE_pl_xfrac3d_"//trim(adjustl(file1))//base_extension
              endif

            endif

          else

            if (do_LTE) then

              if (do_BB) then
                !file1=trim(adjustl(results_dir))//"adaptive_do_LTE_bb_xfrac3d_"//trim(adjustl(file1))//base_extension
              else
                !file1=trim(adjustl(results_dir))//"adaptive_do_LTE_pl_xfrac3d_"//trim(adjustl(file1))//base_extension
              endif

            else

              if (do_BB) then
                !file1=trim(adjustl(results_dir))//"adaptive_no_LTE_bb_xfrac3d_"//trim(adjustl(file1))//base_extension
              else
                !file1=trim(adjustl(results_dir))//"adaptive_no_LTE_pl_xfrac3d_"//trim(adjustl(file1))//base_extension
              endif

            endif

          endif
         endif



          ! Open, write and close
          open(unit=52,file=file1,form="unformatted",status="unknown")
          write(52) mesh,mesh,mesh
          write(52) (((xHII_array(i,j,k),i=1,mesh),j=1,mesh), &
               k=1,mesh)
          close(52)

          number_of_ionized_cells = 0
          do i = 1,mesh
            do j = 1,mesh
              do k = 1,mesh
                if (xHI_array(i,j,k).le.1.0e-2) number_of_ionized_cells = number_of_ionized_cells + 1
              enddo
            enddo
          enddo

          ! Construct file name
         if (do_test4) then
          write(file1,"(f6.3)") sim_time/total_simulation_time !Test 4
file1=trim(adjustl(results_dir))//trim(adjustl(file_timestep))//trim(adjustl(file_LTE))//trim(adjustl(file_function))// &
trim(adjustl(file_source))//"mesh_"//trim(adjustl(file_mesh))//"_xfrac3dHe1_"//trim(adjustl(file1))//base_extension

          if (do_constant) then
            if (do_LTE) then
              !file1=trim(adjustl(results_dir))//"constant_do_LTE_t10_T1e5_xfrac3dHe1_"//trim(adjustl(file1))//base_extension
            else
              !file1=trim(adjustl(results_dir))//"constant_no_LTE_t10_T1e5_xfrac3dHe1_"//trim(adjustl(file1))//base_extension
            endif
          else
            if (do_LTE) then
             !file1=trim(adjustl(results_dir))//"adaptive_do_LTE_T1e5_xfrac3dHe1_"//trim(adjustl(file1))//base_extension
            else
             !file1=trim(adjustl(results_dir))//"adaptive_no_LTE_T1e5_xfrac3dHe1_"//trim(adjustl(file1))//base_extension
            endif
          endif
         else
          write(file1,"(f6.3)") zred_now
file1=trim(adjustl(results_dir))//trim(adjustl(file_timestep))//trim(adjustl(file_LTE))//trim(adjustl(file_function))// &
trim(adjustl(file_source))//"mesh_"//trim(adjustl(file_mesh))//"_xfrac3dHe1_"//trim(adjustl(file1))//base_extension

          if (do_constant) then

            if (do_LTE) then

              if (do_BB) then
                !file1=trim(adjustl(results_dir))//"constant_do_LTE_bb_xfrac3dHe1_"//trim(adjustl(file1))//base_extension
              else
                !file1=trim(adjustl(results_dir))//"constant_do_LTE_pl_xfrac3dHe1_"//trim(adjustl(file1))//base_extension
              endif

            else

              if (do_BB) then
                !file1=trim(adjustl(results_dir))//"constant_no_LTE_bb_xfrac3dHe1_"//trim(adjustl(file1))//base_extension
              else
                !file1=trim(adjustl(results_dir))//"constant_no_LTE_pl_xfrac3dHe1_"//trim(adjustl(file1))//base_extension
              endif

            endif

          else

            if (do_LTE) then

              if (do_BB) then
                !file1=trim(adjustl(results_dir))//"adaptive_do_LTE_bb_xfrac3dHe1_"//trim(adjustl(file1))//base_extension
              else
                !file1=trim(adjustl(results_dir))//"adaptive_do_LTE_pl_xfrac3dHe1_"//trim(adjustl(file1))//base_extension
              endif

            else

              if (do_BB) then
                !file1=trim(adjustl(results_dir))//"adaptive_no_LTE_bb_xfrac3dHe1_"//trim(adjustl(file1))//base_extension
              else
                !file1=trim(adjustl(results_dir))//"adaptive_no_LTE_pl_xfrac3dHe1_"//trim(adjustl(file1))//base_extension
              endif

            endif

          endif
         endif


          ! Open, write and close
          open(unit=62,file=file1,form="unformatted",status="unknown")
          write(62) mesh,mesh,mesh
          write(62) (((xHeII_array(i,j,k),i=1,mesh),j=1,mesh), &
               k=1,mesh)
          close(62)

          ! Construct file name
         if (do_test4) then
          write(file1,"(f6.3)") sim_time/total_simulation_time !Test 4
file1=trim(adjustl(results_dir))//trim(adjustl(file_timestep))//trim(adjustl(file_LTE))//trim(adjustl(file_function))// &
trim(adjustl(file_source))//"mesh_"//trim(adjustl(file_mesh))//"_xfrac3dHe2_"//trim(adjustl(file1))//base_extension

          if (do_constant) then
            if (do_LTE) then
              !file1=trim(adjustl(results_dir))//"constant_do_LTE_t10_T1e5_xfrac3dHe2_"//trim(adjustl(file1))//base_extension
            else
              !file1=trim(adjustl(results_dir))//"constant_no_LTE_t10_T1e5_xfrac3dHe2_"//trim(adjustl(file1))//base_extension
            endif
          else
            if (do_LTE) then
             !file1=trim(adjustl(results_dir))//"adaptive_do_LTE_T1e5_xfrac3dHe2_"//trim(adjustl(file1))//base_extension
            else
             !file1=trim(adjustl(results_dir))//"adaptive_no_LTE_T1e5_xfrac3dHe2_"//trim(adjustl(file1))//base_extension
            endif
          endif
         else
          write(file1,"(f6.3)") zred_now
file1=trim(adjustl(results_dir))//trim(adjustl(file_timestep))//trim(adjustl(file_LTE))//trim(adjustl(file_function))// &
trim(adjustl(file_source))//"mesh_"//trim(adjustl(file_mesh))//"_xfrac3dHe2_"//trim(adjustl(file1))//base_extension

          if (do_constant) then

            if (do_LTE) then

              if (do_BB) then
                !file1=trim(adjustl(results_dir))//"constant_do_LTE_bb_xfrac3dHe2_"//trim(adjustl(file1))//base_extension
              else
                !file1=trim(adjustl(results_dir))//"constant_do_LTE_pl_xfrac3dHe2_"//trim(adjustl(file1))//base_extension
              endif

            else

              if (do_BB) then
                !file1=trim(adjustl(results_dir))//"constant_no_LTE_bb_xfrac3dHe2_"//trim(adjustl(file1))//base_extension
              else
                !file1=trim(adjustl(results_dir))//"constant_no_LTE_pl_xfrac3dHe2_"//trim(adjustl(file1))//base_extension
              endif

            endif

          else

            if (do_LTE) then

              if (do_BB) then
               !file1=trim(adjustl(results_dir))//"adaptive_do_LTE_bb_xfrac3dHe2_"//trim(adjustl(file1))//base_extension
              else
               !file1=trim(adjustl(results_dir))//"adaptive_do_LTE_pl_xfrac3dHe2_"//trim(adjustl(file1))//base_extension
              endif

            else

              if (do_BB) then
               !file1=trim(adjustl(results_dir))//"adaptive_no_LTE_bb_xfrac3dHe2_"//trim(adjustl(file1))//base_extension
              else
               !file1=trim(adjustl(results_dir))//"adaptive_no_LTE_pl_xfrac3dHe2_"//trim(adjustl(file1))//base_extension
              endif

            endif

          endif
         endif


          ! Open, write and close
          open(unit=72,file=file1,form="unformatted",status="unknown")
          write(72) mesh,mesh,mesh
          write(72) (((xHeIII_array(i,j,k),i=1,mesh),j=1,mesh), &
               k=1,mesh)
          close(72)

       else
          ! Report error
          write(logf,*) "Calling stream 2 output where we should not."
       endif
       
    endif
    
  end subroutine write_stream2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> produces output for a time frame. See below for format
  subroutine write_stream3 (zred_now,sim_time)

    implicit none

    real(kind=dp),intent(in) :: zred_now !< current redshift
    real(kind=dp),intent(in) :: sim_time !< current simulation time
    !character(len=*),parameter :: basename="Ifront1_"
    character(len=*),parameter :: base_extension=".bin"
    character(len=80) :: file1
    integer :: i,j,k

    ! Stream 2
    if (rank == 0) then

       if (streams(3) == 1) then

         if (do_test4) then
          write(file1,"(f6.3)") sim_time/total_simulation_time !Test 4
         else
          write(file1,"(f6.3)") zred_now
         endif
          file1=trim(adjustl(results_dir))//"IonRates3D_"//trim(adjustl(file1))//base_extension

          open(unit=53,file=file1,form="unformatted",status="unknown")
          write(53) mesh,mesh,mesh
          write(53) (((real(photoionization_HI_array(i,j,k)),i=1,mesh),j=1,mesh), &
               k=1,mesh)
          close(53)

         if (do_test4) then
          write(file1,"(f6.3)") sim_time/total_simulation_time !Test 4
         else
          write(file1,"(f6.3)") zred_now
         endif
          file1=trim(adjustl(results_dir))//"HeatRates3D_"//trim(adjustl(file1))//base_extension

          open(unit=53,file=file1,form="unformatted",status="unknown")
          write(53) mesh,mesh,mesh
          write(53) (((real(heating_array(i,j,k)),i=1,mesh),j=1,mesh), &
               k=1,mesh)
          close(53)

      ! new thing
      write(96,*) zred_now, &
                  grand_total_ionization_photons, &
                  grand_total_recombination_photons, &
                  grand_total_collisional_photons, &
                  grand_total_recombation_ionization_He, &
                  grand_total_source_photons, &
                  grand_total_escaped_photons

       else
          ! Report error
          write(logf,*) "Calling stream 3 output where we should not."
       endif

    endif


  end subroutine write_stream3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> produces output for a time frame. See below for format
  subroutine write_stream4 (zred_now,sim_time)

    implicit none

    real(kind=dp),intent(in) :: zred_now !< current redshift
    real(kind=dp),intent(in) :: sim_time !< current simulation time
    !character(len=*),parameter :: basename="Ifront1_"
    character(len=*),parameter :: base_extension=".bin"
    character(len=80) :: file0,file1,file2,file3
    integer :: i,j,k

    ! Stream 2
    if (rank == 0) then

       ! Stream 4
       if (streams(4).eq.1) then

         if (do_test4) then
          write(file0,"(f6.3)") sim_time/total_simulation_time !Test 4
         else
          write(file0,"(f6.3)") zred_now
         endif
          file1=trim(adjustl(results_dir))//"Ifront2d_xy_"//trim(adjustl(file0))//".bin"
          file2=trim(adjustl(results_dir))//"Ifront2d_xz_"//trim(adjustl(file0))//".bin"
          file3=trim(adjustl(results_dir))//"Ifront2d_yz_"//trim(adjustl(file0))//".bin"


          ! xy cut through source 
          open(unit=54,file=file1,form="unformatted",status="unknown")
          write(54) mesh,mesh
          write(54) ((real(xHII_array(i,j,mesh/2)),i=1,mesh),j=1,mesh)
          close(54)

          ! xz cut through source 
          open(unit=55,file=file2,form="unformatted",status="unknown")
          write(55) mesh,mesh
          write(55) ((real(xHII_array(i,mesh/2,k)),i=1,mesh),k=1,mesh)
          close(55)

          ! yz cut through source 
          open(unit=56,file=file3,form="unformatted",status="unknown")
          write(56) mesh,mesh
          write(56) ((real(xHII_array(mesh/2,j,k)),j=1,mesh),k=1,mesh)
          close(56)
       else
          ! Report error
          write(logf,*) "Calling stream 4 output where we should not."
       endif
       
    endif

  end subroutine write_stream4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> produces output for a time frame. See below for format
  subroutine write_stream5 (zred_now,sim_time)

    implicit none

    real(kind=dp),intent(in) :: zred_now !< current redshift
    real(kind=dp),intent(in) :: sim_time !< current simulation time
    !character(len=*),parameter :: basename="Ifront1_"
    character(len=*),parameter :: base_extension=".bin"
    character(len=80) :: file0,file1,file2,file3
    character(len=6) :: zred_str
    integer :: i,j,k

    if (rank == 0) then

       ! Stream 5
       if (streams(5) == 1) then

         if (do_test4) then
          write(file0,"(f6.3)") sim_time/total_simulation_time !Test 4
         else
          write(file0,"(f6.3)") zred_now
         endif
          file1=trim(adjustl(results_dir))//"ndens_xy_"//trim(adjustl(file0))//".bin"
          file2=trim(adjustl(results_dir))//"ndens_xz_"//trim(adjustl(file0))//".bin"
          file3=trim(adjustl(results_dir))//"ndens_yz_"//trim(adjustl(file0))//".bin"

          ! xy cut through source 
          open(unit=57,file=file1,form="unformatted",status="unknown")
          write(57) mesh,mesh
          write(57) ((real(number_density_array(i,j,mesh/2)),i=1,mesh),j=1,mesh)
          close(57)

          ! xz cut through source 
          open(unit=58,file=file2,form="unformatted",status="unknown")
          write(58) mesh,mesh
          write(58) ((real(number_density_array(i,mesh/2,k)),i=1,mesh),k=1,mesh)
          close(58)

          ! yz cut through source 
          open(unit=59,file=file3,form="unformatted",status="unknown")
          write(59) mesh,mesh
          write(59) ((real(number_density_array(mesh/2,j,k)),j=1,mesh),k=1,mesh)
          close(59)
       else
          ! Report error
          write(logf,*) "Calling stream 5 output where we should not."
       endif
       
    endif

  end subroutine write_stream5
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> produces output for a time frame. See below for format
  subroutine write_photonstatistics (zred_now,sim_time,dt,photcons_flag)

    implicit none

    real(kind=dp),intent(in) :: zred_now !< current redshift
    real(kind=dp),intent(in) :: sim_time !< current simulation time
    real(kind=dp),intent(in) :: dt !< time step taken
    integer,intent(out) :: photcons_flag

    real(kind=dp) :: totalsrc,photcons,total_photon_loss
    real(kind=dp) :: totions,totphots,volfrac(0:2),massfrac(0:2)

#ifdef MPI
    integer :: mympierror
#endif

    if (rank == 0) then

    endif
    
#ifdef MPI
    call MPI_BCAST(photcons_flag,1,MPI_INTEGER,0,MPI_COMM_NEW,mympierror)
#endif

  end subroutine write_photonstatistics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module output_module
