!>
!! \brief Main to calculate properties of sources
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date: 02-July-2013
!<
Program C2Ray

  ! Goal:
  ! Calculate properties of sources in the same way as done in the C2-Ray
  ! code.

  ! Needs following `modules'
  ! c2ray_parameters : all the tunable parameters for the code
  ! my_mpi : sets up the MPI parallel environment
  ! output_module : output routines
  ! grid : sets up the grid
  ! radiation : radiation tools
  ! nbody : interface to N-body output
  ! cosmology : cosmological utilities
  ! material : material properties
  ! times : time and time step utilities
  ! sourceprops : source properties
  ! evolve : evolve grid in time

  use precision, only: dp
  use clocks, only: setup_clocks, update_clocks, report_clocks
  use file_admin, only: stdinput, logf, file_input, flag_for_file_input
  use my_mpi !, only: mpi_setup, mpi_end, rank
  use radiation, only: rad_ini
  !use radiation, only: hphotint,hphotintthin,photint,photintthin
  !use radiation, only: pl_hphotint,pl_hphotintthin,pl_photint,pl_photintthin
  use radiation, only: bb_photo_thick_table, bb_photo_thin_table
  use radiation, only: pl_photo_thick_table, pl_photo_thin_table
  use radiation, only: bb_heat_thick_table, bb_heat_thin_table
  use radiation, only: pl_heat_thick_table, pl_heat_thin_table
 
#ifdef XLF
  ! Place for xlf specific statements
#endif

  implicit none

#ifdef PGI
  include 'lib3f.h' ! for iargc, getargc
#endif

  integer :: ierror !< error flag
#ifdef MPI
  integer :: mympierror
#endif

  ! Input file
  character(len=512) :: inputfile !< name of input file
  character(len=1) :: answer !< y or n answer
  ! Initialize clocks (cpu and wall)
  call setup_clocks

  ! Set up MPI structure
  call mpi_setup()

  ! Set up input stream (either standard input or from file given
  ! by first argument)
  if (rank == 0) then
     flush(logf)
     if (COMMAND_ARGUMENT_COUNT () > 0) then
        call GET_COMMAND_ARGUMENT(1,inputfile)
        write(logf,*) "reading input from ",trim(adjustl(inputfile))
        open(unit=stdinput,file=inputfile)
        call flag_for_file_input(.true.)
     else
        write(logf,*) "reading input from command line"
     endif
     flush(logf)
  endif

#ifdef MPILOG
  write(logf,*) 'about to initialize rad'
#endif
  ! Initialize photo-ionization calculation
  call rad_ini ( )

  ! Write tables to file
  ! open(unit=101,file="photint.bin",form="unformatted",status="new")
  ! write(101) shape(photint)
  ! write(101) photint
  ! close(101)

  ! open(unit=101,file="hphotint.bin",form="unformatted",status="new")
  ! write(101) shape(hphotint)
  ! write(101) hphotint
  ! close(101)

  ! open(unit=101,file="photintthin.bin",form="unformatted",status="new")
  ! write(101) shape(photintthin)
  ! write(101) photintthin
  ! close(101)

  ! open(unit=101,file="hphotintthin.bin",form="unformatted",status="new")
  ! write(101) shape(hphotintthin)
  ! write(101) hphotintthin
  ! close(101)

  ! open(unit=101,file="pl_photint.bin",form="unformatted",status="new")
  ! write(101) shape(pl_photint)
  ! write(101) pl_photint
  ! close(101)

  ! open(unit=101,file="pl_hphotint.bin",form="unformatted",status="new")
  ! write(101) shape(pl_hphotint)
  ! write(101) pl_hphotint
  ! close(101)

  ! open(unit=101,file="pl_photintthin.bin",form="unformatted",status="new")
  ! write(101) shape(pl_photintthin)
  ! write(101) pl_photintthin
  ! close(101)

  ! open(unit=101,file="pl_hphotintthin.bin",form="unformatted",status="new")
  ! write(101) shape(pl_hphotintthin)
  ! write(101) pl_hphotintthin
  ! close(101)

  ! Write tables to file
  open(unit=101,file="bb_photo_thin_table.bin",form="unformatted",status="new")
  write(101) shape(bb_photo_thin_table)
  write(101) bb_photo_thin_table
  close(101)

  open(unit=101,file="bb_photo_thick_table.bin",form="unformatted",status="new")
  write(101) shape(bb_photo_thick_table)
  write(101) bb_photo_thick_table
  close(101)

  open(unit=101,file="bb_heat_thin_table.bin",form="unformatted",status="new")
  write(101) shape(bb_heat_thin_table)
  write(101) bb_heat_thin_table
  close(101)

  open(unit=101,file="bb_heat_thick_table.bin",form="unformatted",status="new")
  write(101) shape(bb_heat_thick_table)
  write(101) bb_heat_thick_table
  close(101)

  ! Report clocks (cpu and wall)
  call report_clocks ()

  ! End the run
  call mpi_end ()

end Program C2Ray
