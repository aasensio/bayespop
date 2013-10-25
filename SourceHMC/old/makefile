# Makefile for building the MCMC F90 code BayesME

# Select compiler
# Intel Fortran compiler
COMPILER_F90 = ifort
OPTIONS_F90 = -c -O3 
# -g -traceback

COMPILER_C = icc
CC_LANG_FLAGS=-std=c89
CC_OPT_FLAGS=-O3 -ip -parallel
CFLAGS=$(CC_LANG_FLAGS) $(CC_OPT_FLAGS)

COMPILER_C = gcc
CC_LANG_FLAGS=-std=c89
CC_OPT_FLAGS=-O3 -fomit-frame-pointer
CFLAGS=$(CC_LANG_FLAGS) $(CC_OPT_FLAGS)


HMC_LIB = hmc
GEODESIC_LIB = geodesic
BLAS_LIB = /scratch/LIBS/BLAS
LAPACK_LIB = /scratch/LIBS/lapack-3.2

# GFortran
#COMPILER = gfortran
#OPTIONS = -c -frecord-marker=4 -ffree-form

LIBS=-L$(HMC_LIB) -lhmc -L$(GEODESIC_LIB) -lgeodesic -L$(BLAS_LIB) -L$(LAPACK_LIB) -llapack_LINUX -lblas_LINUX

.SUFFIXES: .f .o

.f.o:
	$(COMPILER_F90) $(OPTIONS_F90) $< -o $@
	
.SUFFIXES: .f90 .o

.f90.o:
	$(COMPILER_F90) $(OPTIONS_F90) $< -o $@

.SUFFIXES: .c .o

.c.o:
	$(COMPILER_C) $(CFLAGS) -c $< -o $@

OBJECTS_HMC = hmc/guided_hmc.o hmc/mt19937.o hmc/hanson.o

OBJECTS_GEODESIC = geodesic/Accept.o  geodesic/converge.o  geodesic/destsv.o  geodesic/dgqt.o  geodesic/dpmpar.o  geodesic/fdavv.o \
	geodesic/fdjac.o  geodesic/lambda.o  geodesic/leastsq.o  geodesic/UpdateJac.o

OBJECTS = chk_g95.o chainAnalysis.o like.o maths.o vars.o init.o main.o sampling.o maxlike.o geodesiclm.o

mcmc: $(OBJECTS)
	$(COMPILER_F90) $(OBJECTS) -o mcmc $(LIBS)
	cp mcmc ../
	
all:
	make hmc
	make geodesic
	make mcmc
	
hmc: $(OBJECTS_HMC)
	ar cr hmc/libhmc.a $(OBJECTS_HMC)
	
geodesic: $(OBJECTS_GEODESIC)
	ar cr geodesic/libgeodesic.a $(OBJECTS_GEODESIC)

clean:
	find . -maxdepth 3 -name "*.o" -delete ; find . -maxdepth 3 -name "*.mod" -delete
	find . -maxdepth 1 -name "*.f90~" -delete ; find . -maxdepth 1 -name "mcmc" -delete
	find ../ -maxdepth 1 -name "mcmc" -delete ; find . -maxdepth 2 -name "*.a" -delete

maxlike.o: maxlike.f90 like.o vars.o geodesiclm.o
	$(COMPILER_F90) $(OPTIONS_F90) maxlike.f90
	
sampling.o: sampling.f90 like.o vars.o maths.o
	$(COMPILER_F90) $(OPTIONS_F90) sampling.f90

geodesiclm.o: geodesiclm.f90 like.o
	$(COMPILER_F90) $(OPTIONS_F90) geodesiclm.f90
	
chk_g95.o: chk_g95.f90
	$(COMPILER_F90) $(OPTIONS_F90) chk_g95.f90
	
maths.o: maths.f90
	$(COMPILER_F90) $(OPTIONS_F90) maths.f90

chainAnalysis.o: chainAnalysis.f90 vars.o maths.o
	$(COMPILER_F90) $(OPTIONS_F90) chainAnalysis.f90

vars.o: vars.f90
	$(COMPILER_F90) $(OPTIONS_F90) vars.f90

init.o: init.f90 maths.o vars.o
	$(COMPILER_F90) $(OPTIONS_F90) init.f90

like.o: like.f90 vars.o maths.o
	$(COMPILER_F90) $(OPTIONS_F90) like.f90
	
main.o: main.f90 init.o sampling.o maths.o vars.o chainAnalysis.o maxlike.o
	$(COMPILER_F90) $(OPTIONS_F90) main.f90