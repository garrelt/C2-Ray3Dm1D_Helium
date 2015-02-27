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

##
FC = gfortran # GNU compiler
MPIFC = mpif90.openmpi # MPI compiler
##


# F90 options (ifort)

##
GFORTFLAGS = -cpp -O0 -g -DGFORT
##
 
#GFORTFLAGS = -cpp -O2 -finline-functions
#GFORTFLAGS = -cpp -O3 -DGFORT
# Processor dependent optimization

##
F90FLAGS1 = $(GFORTFLAGS) 
##

##
MY_OPTIONS = -ffpe-trap=invalid,zero,overflow -fcheck=bounds,mem
##

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -I/usr/include/lam -DMPI # For LAM mpi (Stockholm)

##
MPI_FLAGS = -DMPI # 
##

#MPI_FLAGS = -DMPI -DMPILOG # Add more (MPI node) diagnostic output

##
OPENMP_FLAGS = -fopenmp -DMY_OPENMP # For Intel compiler
##

#-------------------------------------------------------
# Compiler
# Intel: best tested
#FC = ifort # Intel compiler
#MPIFC = mpif90 # MPI compiler

# F90 options (ifort)
#IFORTFLAGS = -DMPILOG -O0 -g -u -fpe0 -DIFORT -shared-intel -check all -traceback -fpp  #!<
#IFORTFLAGS = -fpp -O3 -u -fpe0 -ipo -DIFORT -shared-intel -check all -traceback -vec_report
#IFORTFLAGS = -O3 -vec_report -u -fpe0 -ipo -mcmodel=medium -shared-intel -DIFORT #-check all -traceback
# Processor dependent optimization
#F90FLAGS1 = $(IFORTFLAGS)  #!<
#F90FLAGS1 = -xW $(IFORTFLAGS) 
#F90FLAGS1 = -xO $(IFORTFLAGS) 
#F90FLAGS1 = -xT $(IFORTFLAGS) # Laptop 
#F90FLAGS1 = -xB $(IFORTFLAGS)
#F90FLAGS1 = -xAVX $(IFORTFLAGS) # Curie thin nodes
#F90FLAGS1 = -xHost $(IFORTFLAGS) # Curie xlarge nodes 

# These flags should be added to the F90FLAGS1 depending on the executable
# made. Specify this below on a per executable basis.
#MPI_FLAGS = -I/usr/include/lam -DMPI # For LAM mpi (Stockholm)
#MPI_FLAGS = -DMPI # 
#MPI_FLAGS = -DMPI -DMPILOG # Add more (MPI node) diagnostic output
#OPENMP_FLAGS = -openmp -DMY_OPENMP # For Intel compiler

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

LDFLAGS = $(OPTIONS) #-L/afs/astro.su.se/pkg/intel/Compiler/11.1/056/lib/intel64/
LIBS = #-lirc

#-------------------------------------------------------

UTILS= precision.o mrgrnk.o ctrper.o file_admin.o

CONSTANTS = mathconstant.o cgsconstants.o  cgsphotoconstants.o  cgsastroconstants.o c2ray_parameters.o cosmoparms.o abundances.o atomic.o parameter.o type_definition.o 

EVOLVE = evolve_array.o 

GENERAL_KYL = cosmology.o cooling_h.o thermal.o time.o doric.o $(AT) short.o total_photo_rate.o ionization_equation.o local_evolution.o local_photo_rate.o photonstatistics.o output.o 

TEST= clocks.o input.o array.o grid.o tped.o material_array.o density.o temperature.o xfrac.o redshift.o source.o clumping.o scaling_factor.o nominal_source.o photo_table.o photo_rate.o 

AT = AT_array.o AT_ifront.o AT_transformation.o AT_octant_grid.o AT_triangular_grid.o AT_ray.o AT_photoionization.o AT_time_calculation.o

ADP = ADP_array.o ADP_convergence.o ADP_copy.o $(EVOLVE) ADP_evolution.o global_column_density.o global_photo_rate.o

LTE = LTE_array.o LTE_copy.o LTE_evolution.o LTE_field.o LTE_FOF.o LTE_convergence.o

column_density = column_density_0D.o column_density_1D.o column_density_2D.o column_density_3D.o

#--------TEST----------------------------------------------------------------

a: F90=$(FC)
a: F90FLAGS = $(F90FLAGS1)
a: $(UTILS) $(CONSTANTS) romberg.o no_mpi.o $(TEST) $(GENERAL_KYL) $(column_density) $(LTE) $(ADP) column_density_and_photo_rate.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS) $(CONSTANTS) romberg.o no_mpi.o $(TEST) $(GENERAL_KYL) $(column_density) $(LTE) $(ADP) column_density_and_photo_rate.o C2Ray.o

omp: F90=$(FC)
omp: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS)
omp: $(UTILS)  $(CONSTANTS) romberg.o no_mpi.o $(TEST) $(GENERAL_KYL) $(column_density) $(LTE) $(ADP) column_density_and_photo_rate.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS) $(CONSTANTS) romberg.o no_mpi.o $(TEST) $(GENERAL_KYL) $(column_density) $(LTE) $(ADP) column_density_and_photo_rate.o C2Ray.o

C2Ray_3D_test_kyl_periodic_mpi: F90=$(MPIFC)
C2Ray_3D_test_kyl_periodic_mpi: F90FLAGS = $(F90FLAGS1) $(MPI_FLAGS)
C2Ray_3D_test_kyl_periodic_mpi: $(UTILS) $(CONSTANTS) romberg.o mpi.o $(TEST) $(GENERAL_KYL) $(column_density) $(LTE) $(ADP) column_density_and_photo_rate.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS) $(CONSTANTS) romberg.o mpi.o $(TEST) $(GENERAL_KYL) $(column_density) $(LTE) $(ADP) column_density_and_photo_rate.o C2Ray.o

p: F90=$(MPIFC)
p: F90FLAGS = $(F90FLAGS1) $(OPENMP_FLAGS) $(MPI_FLAGS)
p: $(UTILS) $(CONSTANTS) romberg.o mpi.o $(TEST) $(GENERAL_KYL) $(column_density) $(LTE) $(ADP) column_density_and_photo_rate.o C2Ray.o
	$(F90) $(OPTIONS) -o $@ $(UTILS) $(CONSTANTS) romberg.o mpi.o $(TEST) $(GENERAL_KYL) $(column_density) $(LTE) $(ADP) column_density_and_photo_rate.o C2Ray.o



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


