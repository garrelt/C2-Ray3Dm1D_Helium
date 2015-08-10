# Makefile for C2Ray_3D.
#
# Author: Garrelt Mellema

# This Makefile can make different versions of C2Ray_3D.
# These versions differ in their parallelization and/or
# in their connection to specific N-body results.
#
# Note 1: Parallelization
# The parallelization intended is specified in the name
# of the executable: _omp means OpenMP (shared memory), 
# _mpi means MPI (distributed memory). Both can also be
# used at the same time (if your architecture supports
# it.
#
# Note 2: N-body module
# Different versions exist with different Nbody interfaces:
# pmfast - interface to older pmfast simulations
# cubep3m - interface to cubep3m simulations
# LG - interface to Local Group simulations (GADGET)
# Gadget - interface to LOFAR EoR GADGET simulations (not working)
#
# Note 3: Compiler & Flags
# The compiler is specified by the FC variable (MPIFC for the MPI
# compiler). We have only extensively used the Intel F90 compiler. 
# Support for other compilers will have to be added.
# Parts of the code need to know about the compiler, this is
# done through preprocessor statements. So when compiling with
# intel compiler, -DIFORT needs to be specified. Support for
# new compilers thus needs to be added in the code too.
#
# Note 4: Recompiling
# Some dependencies are through module parameters, and thus
# not recognized by make. Best practise is to run "make clean"
# before running "make".
#-------------------------------------------------------

# Compiler
#FC = gfortran # GNU compiler
#MPIFC = mpif90.openmpi # MPI compiler

# F90 options (ifort)
#GFORTFLAGS = -cpp -O0 -g -DGFORT 
#GFORTFLAGS = -O3 -DGFORT
# Processor dependent optimization
#F90FLAGS1 = $(GFORTFLAGS) 

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -I/usr/include/lam -DMPI # For LAM mpi (Stockholm)
#MPI_FLAGS = -DMPI # 
#MPI_FLAGS = -DMPI -DMPILOG # Add more (MPI node) diagnostic output
#OPENMP_FLAGS = -fopenmp -DMY_OPENMP # For Intel compiler

#-------------------------------------------------------
# Compiler
# Intel: best tested
FC = ifort # Intel compiler
MPIFC = mpif90 # MPI compiler

# F90 options (ifort)
#IFORTFLAGS = -DMPILOG -O0 -g -u -fpe0 -DIFORT -shared-intel -check all -traceback -fpp  #!<
#IFORTFLAGS = -fpp -O3 -u -fpe0 -ipo -DIFORT -shared-intel -check all -traceback -vec_report
IFORTFLAGS = -fpp -O3 -u -fpe0 -ipo -mcmodel=medium -shared-intel -DIFORT #-check all -traceback
# Processor dependent optimization
#F90FLAGS1 = $(IFORTFLAGS)  #!<
F90FLAGS1 = -xHost $(IFORTFLAGS)  #!<
#F90FLAGS1 = -xW $(IFORTFLAGS) 
#F90FLAGS1 = -xO $(IFORTFLAGS) 
#F90FLAGS1 = -xT $(IFORTFLAGS) # Laptop 
#F90FLAGS1 = -xB $(IFORTFLAGS)
#F90FLAGS1 = -xAVX $(IFORTFLAGS) # Curie thin nodes
#F90FLAGS1 = -xHost $(IFORTFLAGS) # Curie xlarge nodes 

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -I/usr/include/lam -DMPI # For LAM mpi (Stockholm)
MPI_FLAGS = -DMPI # 
#MPI_FLAGS = -DMPI -DMPILOG # Add more (MPI node) diagnostic output
OPENMP_FLAGS = -openmp -DMY_OPENMP # For Intel compiler

#-------------------------------------------------------

# Compiler
# Sun: problems with constant definition. Cannot have sqrt in constant
# definition.
#FC = f95 # Sun compiler
#MPIFC = mpif90 # MPI compiler

# F90 options (ifort)
#SUNFLAGS = -O3 -DSUN
# Processor dependent optimization
#F90FLAGS1 = $(SUNFLAGS) 
#F90FLAGS1 = -xW $(SUNFLAGS) 

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -I/usr/include/lam -DMPI # For LAM mpi (Stockholm)
#MPI_FLAGS = -DMPI # 
#MPI_FLAGS = $(MPI_FLAGS) -DMPILOG # Add more (MPI node) diagnostic output
#OPENMP_FLAGS = -openmp -DMY_OPENMP # For Sun compiler

#-------------------------------------------------------

# PGI compiler
#FC = pf90
#MPIFC = mpif77
#MPIFC = mpif90

# F90 options (pgi)
#PGIFLAGS = -O3 -fast -DPGI
#F90FLAGS1 = -tp barcelona-64  $(PGIFLAGS) # ranger processors

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -DMPI 
#MPI_FLAGS = $(MPI_FLAGS) -DMPILOG # Add more (MPI node) diagnostic output
#OPENMP_FLAGS = -mp -DMY_OPENMP

#-------------------------------------------------------

# Other F90 compilers
#FC = f90

#-------------------------------------------------------

#LDR     = $(F90)

OPTIONS = $(F90FLAGS)

LDFLAGS = $(OPTIONS) -L/afs/astro.su.se/pkg/intel/Compiler/11.1/056/lib/intel64/
LIBS = -lirc

#-------------------------------------------------------

UTILS1= ../precision.o ../mrgrnk.o ../ctrper.o ../romberg.o ../file_admin.o ../report_memory.o
UTILS2= precision.o mrgrnk.o ctrper.o romberg.o file_admin.o report_memory.o

CONSTANTS1 = ../mathconstants.o ../cgsconstants.o  ../cgsphotoconstants.o  ../cgsastroconstants.o ../c2ray_parameters.o ../cosmoparms.o ../abundances.o ../atomic.o

CONSTANTS2 = mathconstants.o cgsconstants.o  cgsphotoconstants.o  cgsastroconstants.o c2ray_parameters.o cosmoparms.o abundances.o atomic.o

CONSTANTS1_TEST4 = ../mathconstants.o ../cgsconstants.o  ../cgsphotoconstants.o  ../cgsastroconstants.o ../c2ray_parameters_TEST4.o ../cosmoparms.o ../abundances.o ../atomic.o

CONSTANTS2_TEST4 = mathconstants.o cgsconstants.o  cgsphotoconstants.o  cgsastroconstants.o c2ray_parameters_TEST4.o cosmoparms.o abundances.o atomic.o


#RADIATION1 = ../radiation.o
#RADIATION2 = radiation.o
RADIATION1 = ../radiation_sizes.o ../radiation_sed_parameters.o ../radiation_tables.o ../radiation_photoionrates.o
RADIATION2 = radiation_sizes.o radiation_sed_parameters.o radiation_tables.o radiation_photoionrates.o

#EVOLVE = evolve8.o
EVOLVE = evolve_data.o column_density.o evolve_point.o evolve_source.o master_slave.o evolve.o

GENERAL1=../cosmology.o  ../cooling_h.o ../radiation_including_pls_new.o   ../thermal.o ../time_ini.o ../doric.o photonstatistics.o $(EVOLVE) output.o

GENERAL2=cosmology.o  cooling_h.o radiation_including_pls_new.o   thermal.o time_ini.o doric.o photonstatistics.o $(EVOLVE) output.o

GENERAL1_KYL=../cosmology.o  ../cooling_h.o $(RADIATION1)   ../thermal.o ../time_ini.o ../doric.o photonstatistics.o $(EVOLVE) output.o

GENERAL2_KYL=cosmology.o  cooling_h.o $(RADIATION2) thermal.o time_ini.o doric.o photonstatistics.o $(EVOLVE) output.o

TEST= ../clocks.o  ../test.o grid.o ../tped.o mat_ini_test.o sourceprops_test.o

TEST_2= clocks.o  test.o grid.o tped.o mat_ini_test.o sourceprops_test.o

TEST1= ../clocks.o  ../test.o grid.o ../tped.o mat_ini_test.o sourceprops_test2.o

TEST2=  clocks.o  test.o grid.o tped.o mat_ini_test.o sourceprops_test2.o

TEST4_1= ../clocks.o  ../test4.o grid.o ../tped.o mat_ini_test4.o sourceprops_test4.o

TEST4_2=  clocks.o  test4.o grid.o tped.o mat_ini_test4.o sourceprops_test4.o

PMFAST1= ../clocks.o pmfast.o grid.o ../tped.o mat_ini_pmfast.o sourceprops_pmfast.o
PMFAST2= clocks.o pmfast.o grid.o tped.o mat_ini_pmfast.o sourceprops_pmfast.o

CUBEP3M1= ../clocks.o ../cubep3m.o grid.o ../tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o
CUBEP3M2= clocks.o cubep3m.o grid.o tped.o mat_ini_cubep3m.o sourceprops_cubep3m.o

LG1= ../clocks.o  LG.o grid.o ../tped.o mat_ini_LG.o sourceprops_LG.o
LG2= clocks.o  LG.o grid.o tped.o mat_ini_LG.o sourceprops_LG.o

GADGET1= ../clocks.o gadget.o grid.o ../tped.o mat_ini_Gadget.o sourceprops_gadget.o
GADGET2= clocks.o gadget.o grid.o tped.o mat_ini_Gadget.o sourceprops_gadget.o

#--------TEST----------------------------------------------------------------

C2Ray_3D_test_periodic: F90=$(FC)
C2Ray_3D_test_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_test_periodic: $(UTILS1) sizes.o $(CONSTANTS1) ../no_mpi.o $(TEST) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2) no_mpi.o $(TEST_2) $(GENERAL2) C2Ray.o

C2Ray_3D_test_kyl_periodic: F90=$(FC)
C2Ray_3D_test_kyl_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_test_kyl_periodic: $(UTILS1) sizes.o $(CONSTANTS1) ../no_mpi.o $(TEST) $(GENERAL1_KYL) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2) no_mpi.o $(TEST_2) $(GENERAL2_KYL) C2Ray.o

C2Ray_3D_test_kyl_periodic_omp: F90=$(FC)
C2Ray_3D_test_kyl_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS)
C2Ray_3D_test_kyl_periodic_omp: $(UTILS1) sizes.o $(CONSTANTS1) ../no_mpi.o $(TEST) $(GENERAL1_KYL) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2) no_mpi.o $(TEST_2) $(GENERAL2_KYL) C2Ray.o

C2Ray_3D_test_kyl_periodic_omp_mpi: F90=$(MPIFC)
C2Ray_3D_test_kyl_periodic_omp_mpi: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS)  $(MPI_FLAGS)
C2Ray_3D_test_kyl_periodic_omp_mpi: $(UTILS1) sizes.o $(CONSTANTS1) ../mpi.o $(TEST) $(GENERAL1_KYL) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2) mpi.o $(TEST_2) $(GENERAL2_KYL) C2Ray.o

C2Ray_3D_test_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_test_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_test_periodic_mpi: $(UTILS1) sizes.o $(CONSTANTS1) ../mpi.o $(TEST) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2) mpi.o $(TEST_2) $(GENERAL2) C2Ray.o

C2Ray_3D_test_kyl_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_test_kyl_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_test_kyl_periodic_mpi: $(UTILS1) sizes.o $(CONSTANTS1) ../mpi.o $(TEST) $(GENERAL1_KYL) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2) mpi.o $(TEST_2) $(GENERAL2_KYL) C2Ray.o

C2Ray_3D_test_periodic_omp: F90=$(FC)
C2Ray_3D_test_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS)
C2Ray_3D_test_periodic_omp: $(UTILS1) sizes.o $(CONSTANTS1) ../no_mpi.o $(TEST) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2) no_mpi.o $(TEST_2) $(GENERAL2) C2Ray.o

#--------TEST2----------------------------------------------------------------
C2Ray_3D_test2_periodic: F90=$(FC)
C2Ray_3D_test2_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_test2_periodic: $(UTILS1) sizes.o $(CONSTANTS1) ../no_mpi.o $(TEST1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@  $(UTILS2) sizes.o $(CONSTANTS2) no_mpi.o $(TEST2) $(GENERAL2) C2Ray.o

C2Ray_3D_test2_kyl_periodic: F90=$(FC)
C2Ray_3D_test2_kyl_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_test2_kyl_periodic: $(UTILS1) sizes.o $(CONSTANTS1) ../no_mpi.o $(TEST1) $(GENERAL1_KYL) C2Ray.o
	$(F90) $(OPTIONS) -o $@  $(UTILS2) sizes.o $(CONSTANTS2) no_mpi.o $(TEST2) $(GENERAL2_KYL) C2Ray.o

C2Ray_3D_test2_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_test2_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_test2_periodic_mpi:  $(UTILS1)  sizes.o $(CONSTANTS1)   mpi.o $(TEST1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   mpi.o $(TEST2) $(GENERAL2) C2Ray.o

C2Ray_3D_test2_kyl_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_test2_kyl_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_test2_kyl_periodic_mpi:  $(UTILS1)  sizes.o $(CONSTANTS1)   mpi.o $(TEST1) $(GENERAL1_KYL) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   mpi.o $(TEST2) $(GENERAL2_KYL) C2Ray.o

C2Ray_3D_test2_periodic_omp: F90=$(FC)
C2Ray_3D_test2_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS)
C2Ray_3D_test2_periodic_omp:  $(UTILS1)  sizes.o $(CONSTANTS1)  ../no_mpi.o $(TEST1) $(GENERAL1_KYL) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)  no_mpi.o $(TEST2) $(GENERAL2_KYL) C2Ray.o

C2Ray_3D_test2_kyl_periodic_omp: F90=$(FC)
C2Ray_3D_test2_kyl_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS)
C2Ray_3D_test2_kyl_periodic_omp:  $(UTILS1)  sizes.o $(CONSTANTS1)  ../no_mpi.o $(TEST1) $(GENERAL1_KYL) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)  no_mpi.o $(TEST2) $(GENERAL2_KYL) C2Ray.o

#--------TEST4----------------------------------------------------------------

C2Ray_3D_test4_periodic: F90=$(FC)
C2Ray_3D_test4_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_test4_periodic: $(UTILS1) sizes.o $(CONSTANTS1_TEST4)    ../no_mpi.o $(TEST4_1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@  $(UTILS2) sizes.o $(CONSTANTS2_TEST4)  no_mpi.o $(TEST4_2) $(GENERAL2) C2Ray.o

C2Ray_3D_test4_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_test4_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_test4_periodic_mpi:  $(UTILS1)  sizes.o $(CONSTANTS1_TEST4)   mpi.o $(TEST4_1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2_TEST4)   mpi.o $(TEST4_2) $(GENERAL2) C2Ray.o



C2Ray_3D_test4_periodic_omp: F90=$(FC)
C2Ray_3D_test4_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS)
C2Ray_3D_test4_periodic_omp:  $(UTILS1)  sizes.o $(CONSTANTS1_TEST4)  ../no_mpi.o $(TEST4_1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2_TEST4)  no_mpi.o $(TEST4_2) $(GENERAL2) C2Ray.o



#--------PMFAST---------------------------------------------------------------

C2Ray_3D_pmfast_periodic: F90=$(FC)
C2Ray_3D_pmfast_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_pmfast_periodic:  $(UTILS1)  sizes.o $(CONSTANTS1)   ../no_mpi.o $(PMFAST1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   no_mpi.o $(PMFAST2) $(GENERAL2) C2Ray.o

C2Ray_3D_pmfast_periodic_compr: F90=$(FC)
C2Ray_3D_pmfast_periodic_compr: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_pmfast_periodic_compr:  $(UTILS1)  sizes.o $(CONSTANTS1)   ../no_mpi.o clocks.o pmfast.o grid.o tped.o mat_ini_pmfast_compr.o sourceprops_pmfast_compr.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics_compr.o evolve4_periodic_compr.o output_compr.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   no_mpi.o clocks.o  pmfast.o grid.o tped.o mat_ini_pmfast_compr.o sourceprops_pmfast_compr.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics_compr.o evolve4_periodic_compr.o output_compr.o C2Ray.o

C2Ray_3D_pmfast_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_pmfast_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_pmfast_periodic_mpi:  $(UTILS1)  sizes.o $(CONSTANTS1)    mpi.o $(PMFAST1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   mpi.o $(PMFAST2) $(GENERAL2) C2Ray.o

C2Ray_3D_pmfast_periodic_omp_mpi: F90=$(MPIFC)
C2Ray_3D_pmfast_periodic_omp_mpi: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) $(MPI_FLAGS)
C2Ray_3D_pmfast_periodic_omp_mpi:  $(UTILS1)  sizes.o $(CONSTANTS1)  mpi.o $(PMFAST1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   mpi.o $(PMFAST2) $(GENERAL2) C2Ray.o

C2Ray_3D_pmfast_periodic_omp: F90=$(FC)
C2Ray_3D_pmfast_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS)
C2Ray_3D_pmfast_periodic_omp:  $(UTILS1)  sizes.o $(CONSTANTS1)    ../no_mpi.o $(PMFAST1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)  no_mpi.o $(PMFAST2) $(GENERAL2) C2Ray.o

#--------CUBEP3M-------------------------------------------------------------

C2Ray_3D_cubep3m_periodic: F90=$(FC)
C2Ray_3D_cubep3m_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_cubep3m_periodic:  $(UTILS1)  sizes.o $(CONSTANTS1)    ../no_mpi.o $(CUBEP3M1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   no_mpi.o $(CUBEP3M2) $(GENERAL2) C2Ray.o

C2Ray_3D_cubep3m_kyl_periodic: F90=$(FC)
C2Ray_3D_cubep3m_kyl_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_cubep3m_kyl_periodic:  $(UTILS1)  sizes.o $(CONSTANTS1)    ../no_mpi.o $(CUBEP3M1) $(GENERAL1_KYL) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   no_mpi.o $(CUBEP3M2) $(GENERAL2_KYL) C2Ray.o

C2Ray_3D_cubep3m_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_cubep3m_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_cubep3m_periodic_mpi:  $(UTILS1)  sizes.o $(CONSTANTS1)   ../mpi.o $(CUBEP3M1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   mpi.o $(CUBEP3M2) $(GENERAL2) C2Ray.o

C2Ray_3D_cubep3m_periodic_omp: F90=$(FC)
C2Ray_3D_cubep3m_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) 
C2Ray_3D_cubep3m_periodic_omp:  $(UTILS1)  sizes.o $(CONSTANTS1)   ../no_mpi.o $(CUBEP3M1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   no_mpi.o $(CUBEP3M2) $(GENERAL2) C2Ray.o

C2Ray_3D_cubep3m_periodic_omp_mpi: F90=$(MPIFC)
C2Ray_3D_cubep3m_periodic_omp_mpi: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) $(MPI_FLAGS)
C2Ray_3D_cubep3m_periodic_omp_mpi:  $(UTILS1)  sizes.o $(CONSTANTS1)  ../mpi.o $(CUBEP3M1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   mpi.o $(CUBEP3M2) $(GENERAL2) C2Ray.o

C2Ray_3D_cubep3m_kyl_periodic_omp_mpi: F90=$(MPIFC)
C2Ray_3D_cubep3m_kyl_periodic_omp_mpi: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) $(MPI_FLAGS)
C2Ray_3D_cubep3m_kyl_periodic_omp_mpi:  $(UTILS1)  sizes.o $(CONSTANTS1)  ../mpi.o $(CUBEP3M1) $(GENERAL1_KYL) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   mpi.o $(CUBEP3M2) $(GENERAL2_KYL) C2Ray.o

C2Ray_3D_cubep3m_periodic_compr_mpi: F90=$(MPIFC)
C2Ray_3D_cubep3m_periodic_compr_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_cubep3m_periodic_compr_mpi:  $(UTILS1)  sizes.o $(CONSTANTS1)    mpi.o clocks.o cubep3m.o grid.o tped.o mat_ini_cubep3m_compr.o sourceprops_cubep3m_compr.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics_compr.o evolve4_periodic_compr.o output_compr.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   mpi.o clocks.o  cubep3m.o grid.o tped.o mat_ini_cubep3m_compr.o sourceprops_cubep3m_compr.o cooling.o radiation.o cosmology.o time_ini.o doric.o photonstatistics_compr.o evolve4_periodic_compr.o output_compr.o C2Ray.o

#--------LG------------------------------------------------------------------
 
C2Ray_3D_LG_periodic: F90=$(FC)
C2Ray_3D_LG_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_LG_periodic:  $(UTILS1)  sizes.o $(CONSTANTS1)    ../no_mpi.o $(LG1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)  no_mpi.o $(LG2) $(GENERAL2) C2Ray.o

C2Ray_3D_LG_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_LG_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_LG_periodic_mpi:  $(UTILS1)  sizes.o $(CONSTANTS1)    mpi.o $(LG1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   mpi.o $(LG2) $(GENERAL2) C2Ray.o

C2Ray_3D_LG_periodic_omp_mpi: F90=$(MPIFC)
C2Ray_3D_LG_periodic_omp_mpi: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) $(MPI_FLAGS)
C2Ray_3D_LG_periodic_omp_mpi:  $(UTILS1)  sizes.o $(CONSTANTS1)    mpi.o $(LG1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   mpi.o $(LG2) $(GENERAL2) C2Ray.o

C2Ray_3D_LG_periodic_omp: F90=$(FC)
C2Ray_3D_LG_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) 
C2Ray_3D_LG_periodic_omp:  $(UTILS1)  sizes.o $(CONSTANTS1)    ../no_mpi.o $(LG1) $(GENERAL1) C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   no_mpi.o $(LG2) $(GENERAL2) C2Ray.o

#--------GADGET-------------------------------------------------------------------------------

C2Ray_3D_Gadget_periodic: F90=$(FC)
C2Ray_3D_Gadget_periodic: F90FLAGS = $(F90FLAGS1)
C2Ray_3D_Gadget_periodic:  $(UTILS1)  sizes.o $(CONSTANTS1)    ../no_mpi.o $(GADGET1) $(GENERAL1) C2Ray_GadgetTest.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   no_mpi.o $(GADGET2) $(GENERAL2) C2Ray_GadgetTest.o

C2Ray_3D_Gadget_periodic_omp: F90=$(FC)
C2Ray_3D_Gadget_periodic_omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS)
C2Ray_3D_Gadget_periodic_omp:  $(UTILS1)  sizes.o $(CONSTANTS1)    ../no_mpi.o $(GADGET1) $(GENERAL1) C2Ray_GadgetTest.o
	$(F90) $(OPTIONS) -o $@ $(UTILS2) sizes.o $(CONSTANTS2)   no_mpi.o $(GADGET2) $(GENERAL2) C2Ray_GadgetTest.o

clean : 
	rm -f *.o *.mod *.l *.il

.f.o:
	$(F90) -c $(OPTIONS) $<

.f90.o:
	$(F90) -c $(OPTIONS) $<

.F90.o:
	$(F90) -c $(OPTIONS) $<

f.mod:
	$(F90) -c $(OPTIONS) $<

.SUFFIXES: .f90 .F90 .mod .o


