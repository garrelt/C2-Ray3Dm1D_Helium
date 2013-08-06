# MAKROS----------------------------------------------------
# F90 compiler

F90 = ifort
#F90 = gfortran

# ifort options
IFORTFLAGS = -O3 -vec_report -u -fpe0 -ipo -DIFORT
#IFORTFLAGS = -g -C -traceback -fpe0 -DIFORT
#IFORTFLAGS = -g -C -DIFORT
#IFORTFLAGS = -pg -g  -O0 -u -fpe0 -check all -traceback -ipo -DIFORT

#F90FLAGS1 = $(IFORTFLAGS) 
#F90FLAGS1 = -axW $(IFORTFLAGS) 
F90FLAGS1 = -axO $(IFORTFLAGS) 

#GFORTFLAGS = -O
#F90FLAGS1 = $(GFORTFLAGS) 

#F90FLAGS = $(F90FLAGS1) -openmp -DMPI
#F90FLAGS = $(F90FLAGS1) -DMPI
#F90FLAGS = $(F90FLAGS1)
#F90FLAGS = $(F90FLAGS1) -openmp

F90FLAGS = $(F90FLAGS1)

#LDR     = $(F90)

OPTIONS = $(F90FLAGS)

#LDFLAGS = $(OPTIONS) -L/afs/astro.su.se/pkg/intel/Compiler/11.1/056/lib/intel64/
#LIBS = -lirc

UTILS1=../romberg.o ../string.o ../mrgrnk.o  ../ctrper.o 

UTILS2= romberg.o string.o mrgrnk.o  ctrper.o

CONSTANTS1 = ../mathconstants.o ../cgsconstants.o  ../cgsphotoconstants.o \
 ../cgsastroconstants.o ../c2ray_parameters.o ../cosmoparms.o ../abundances.o \
../atomic.o ../cosmoparms.o

CONSTANTS2 = mathconstants.o cgsconstants.o  cgsphotoconstants.o \
 cgsastroconstants.o c2ray_parameters.o cosmoparms.o abundances.o atomic.o \
cosmoparms.o

# there are two different versions of doric, one which solves the
# inhomogeneous system of 3 ode (doric2_july2010) and one which
# solves the homogeneous system of 5 ode (doric3b)

#DORIC= doric3b.o 
DORIC= ../doric.o # 

DORIC2= doric.o # 

# there are also some different versions of RADIATION.
# radiation_monocromatic takes different inputs, see input file
# inputs/isochromatic
RADIATION= ../radiation_including_pls_new.o 
#RADIATION= radiation_monocromatic.o 
#RADIATION= ../radiation.o
#RADIATION2= radiation_monocromatic.o
RADIATION2= radiation_including_pls_new.o
#RADIATION2= radiation.o
#-----------------------------------------------------------

# Building C2Ray_1D: 
# $@ means: name of the file to be made 
C2Ray_1D_test: ../precision.o sizes.o $(CONSTANTS1) $(UTILS1) ../file_admin.o \
 ../no_mpi.o   ../file_admin.o ../clocks.o grid.o ../tped.o  sourceprops_test_one_source.o  mat_ini.o ../cosmology.o\
../cooling_h.o $(RADIATION) ../thermal.o ../time.o $(DORIC) \
photonstatistics.o evolve_new.o output.o C2Ray.o


	$(F90) $(OPTIONS) -o $@ precision.o $(CONSTANTS2) $(UTILS2) \
	file_admin.o clocks.o sizes.o no_mpi.o   grid.o \
	tped.o  sourceprops_test_one_source.o  cosmology.o  mat_ini.o cooling_h.o \
	$(RADIATION2) thermal.o time.o $(DORIC2) photonstatistics.o \
	evolve_new.o output.o C2Ray.o

TableWrite: ../precision.o sizes.o $(CONSTANTS1) $(UTILS1) ../file_admin.o \
 ../no_mpi.o   ../file_admin.o ../clocks.o grid.o ../tped.o  sourceprops_test_one_source.o  mat_ini.o ../cosmology.o\
../cooling_h.o $(RADIATION) ../thermal.o ../time.o $(DORIC) \
photonstatistics.o TableWrite.o


	$(F90) $(OPTIONS) -o $@ precision.o $(CONSTANTS2) $(UTILS2) \
	file_admin.o clocks.o sizes.o no_mpi.o   grid.o \
	tped.o  sourceprops_test_one_source.o  cosmology.o  mat_ini.o cooling_h.o \
	$(RADIATION2) thermal.o time.o $(DORIC2) photonstatistics.o \
	TableWrite.o

clean : 
	rm -f *.o *.mod 

# Building object files: 
# $< means: name of related file that caused the action

%.o: %.mod

.f90.o:
	$(F90) -c $(OPTIONS) $<
.F90.o:
	$(F90) -c $(OPTIONS) $<
f.mod:
	$(F90) -c $(OPTIONS) $<
# Suffix rules: List of significant suffixes
.SUFFIXES: .f90 .F90 .mod .o
