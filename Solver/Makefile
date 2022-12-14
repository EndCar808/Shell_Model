# Enda Carroll
# Sept 2022
# Makfile to compile programmes for solver for the Magnetohydrodynamic GOY shell Model

 
# ---------------------------------------------------------------------
#  Compiler flags
# ---------------------------------------------------------------------
# Get the operating system:
# Need to strip leading and/or trailing whitespaces otherwise if statement wont work
OS = $(strip $(shell lsb_release -si))

# Set the operating system dependent compiler variables
ifeq ($(OS), Ubuntu)
	# CHIRP
	GCC := gcc
	CCFLAG := -O3 -W -Wall -g
	LIBFLAG := -lm -lhdf5 -lhdf5_hl
	GSLFLAG := -lgsl -lgslcblas
	# LAPACKFLAG := -llapacke -llapack -lblas -lgfortran
	INCDIR = -I/usr/include/hdf5/serial/  -I/usr/include/gsl/ 
	LIBDIR = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ 
else ifeq ($(OS), RedHatEnterpriseServer)
	# SONIC
	GCC := gcc
	CCFLAG := -O3 -W -Wall -g
	LIBFLAG := -lm -lhdf5 -lhdf5_hl #-fopenmp -lpthread -lfftw3_omp
	INCDIR = 
	LIBDIR = 
else ifeq ($(OS), CentOS)
	# KAY
	GCC := gcc
	CCFLAG := -O3 -W -Wall -g
	LIBFLAG := -lm -lhdf5 -lhdf5_hl #-fopenmp -lpthread -lfftw3_omp
	INCDIR = 
	LIBDIR = 
else
	$(error Unknown OS found, please check OS and add appropriate flags in Makefile)
endif

# ---------------------------------------------------------------------
#  Executable flags
# ---------------------------------------------------------------------
# Flags to be passed to the solver for normal runs
HYDRO_SOLVER_IFRK4_FLAGS := -D__INT_FAC_RK4 -D__STATS #-D__HYPER -D__ADAPTIVE_STEP -D__CFL_STEP #-D__TRANSIENTS  -D__EKMN_DRAG
HYDRO_SOLVER_RK4_FLAGS := -D__RK4 -D__STATS #-D__HYPER -D__ADAPTIVE_STEP -D__CFL_STEP #-D__TRANSIENTS  -D__EKMN_DRAG
HYDRO_SOLVER_AB4CN_FLAGS := -D__AB4CN -D__STATS #-D__HYPER -D__ADAPTIVE_STEP -D__CFL_STEP #-D__TRANSIENTS  -D__EKMN_DRAG

MAG_HYDRO_SOLVER_FLAGS := -D__INT_FAC_RK4 -D__MAGNETO #-D__STATS -D__HYPER -D__ADAPTIVE_STEP -D__CFL_STEP #-D__TRANSIENTS  -D__EKMN_DRAG

# Phase only flags to be passed to solver
PHASE_ONLY_FLAGS := -D__PHASE_ONLY -D__RK4 #-D__MAGNETO -D__HYPER -D__TRANSIENTS

# Phase only flags to be passed to solver
PHASE_ONLY_DIRECT_FLAGS := -D__PHASE_ONLY_DIRECT -D__RK4 #-D__MAGNETO -D__HYPER -D__TRANSIENTS

# Dataset Flags to be passed to the solver
DSET_FLAGS := -D__SYS_MEASURES -D__STATS #-D__ENRG_SPECT -D__DISS_SPECT -D__ENRG_FLUX -D__FORCING
 
# Flags to be passed to the solver for testing
TEST_FLAGS    := -D__MAGNETO -D__TESTING -D__DEBUG #-D__ADAPTIVE_STEP -D__CFL_STEP -D__TRANSIENTS
# ---------------------------------------------------------------------
#  Directory creation
# ---------------------------------------------------------------------
OBJBIN := obj
OBJDIR := $(shell mkdir -p $(OBJBIN))

BIN    := bin
BINDIR := $(shell mkdir -p $(BIN))

OUT    := output
OUTDIR := $(shell mkdir -p $(OUT))

SRCDIR      	:= src
SRCFILE     	:= $(wildcard $(SRCDIR)/*.c)
SRCOBJ_H_IF    	:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_hydro_IF.o, $(SRCFILE))
SRCOBJ_H_RK    	:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_hydro_RK.o, $(SRCFILE))
SRCOBJ_H_AB    	:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_hydro_AB.o, $(SRCFILE))
SRCOBJ_MHD  	:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_mag_hydro.o, $(SRCFILE))
SRCOBJ_PO  		:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_phase_only.o, $(SRCFILE))
SRCOBJ_PO_D 	:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_phase_only_direct.o, $(SRCFILE))
SRCOBJ_Test 	:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_test.o, $(SRCFILE))

DEPS = $(wildcard *.h) 

# ---------------------------------------------------------------------
#  Builds
# ---------------------------------------------------------------------
############
.PHONY: all
############


all: $(BIN)/solver_hydro_IF $(BIN)/solver_hydro_RK $(BIN)/solver_hydro_AB $(BIN)/solver_mag_hydro $(BIN)/solver_phase_only $(BIN)/solver_phase_only_direct $(BIN)/solver_test


####--------
#  Solver
####--------
solver_hydro_IF: $(BIN)/solver_hydro_IF

$(BIN)/solver_hydro_IF: $(SRCOBJ_H_IF) $(DEPS)
	$(GCC) $(CCFLAG) $(HYDRO_SOLVER_IFRK4_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%_hydro_IF.o: $(SRCDIR)/%.c $(DEPS)
	$(GCC) $(CCFLAG) $(HYDRO_SOLVER_IFRK4_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

solver_hydro_RK: $(BIN)/solver_hydro_RK

$(BIN)/solver_hydro_RK: $(SRCOBJ_H_RK) $(DEPS)
	$(GCC) $(CCFLAG) $(HYDRO_SOLVER_RK4_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%_hydro_RK.o: $(SRCDIR)/%.c $(DEPS)
	$(GCC) $(CCFLAG) $(HYDRO_SOLVER_RK4_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

solver_hydro_AB: $(BIN)/solver_hydro_AB

$(BIN)/solver_hydro_AB: $(SRCOBJ_H_AB) $(DEPS)
	$(GCC) $(CCFLAG) $(HYDRO_SOLVER_AB4CN_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%_hydro_AB.o: $(SRCDIR)/%.c $(DEPS)
	$(GCC) $(CCFLAG) $(HYDRO_SOLVER_AB4CN_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

solver_mag_hydro: $(BIN)/solver_mag_hydro 

$(BIN)/solver_mag_hydro: $(SRCOBJ_MHD) $(DEPS)
	$(GCC) $(CCFLAG) $(MAG_HYDRO_SOLVER_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%_mag_hydro.o: $(SRCDIR)/%.c $(DEPS)
	$(GCC) $(CCFLAG) $(MAG_HYDRO_SOLVER_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

####-----------
#  Phase Only 
####-----------
solver_phase_only: $(BIN)/solver_phase_only 

$(BIN)/solver_phase_only: $(SRCOBJ_PO) $(DEPS) 
	$(GCC) $(CCFLAG) $(PHASE_ONLY_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%_phase_only.o: $(SRCDIR)/%.c $(DEPS) 
	$(GCC) $(CCFLAG) $(PHASE_ONLY_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

solver_phase_only_direct: $(BIN)/solver_phase_only_direct 

$(BIN)/solver_phase_only_direct: $(SRCOBJ_PO_D) $(DEPS) 
	$(GCC) $(CCFLAG) $(PHASE_ONLY_DIRECT_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%_phase_only_direct.o: $(SRCDIR)/%.c $(DEPS) 
	$(GCC) $(CCFLAG) $(PHASE_ONLY_DIRECT_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

####--------
#  Testing
####--------
test: $(BIN)/solver_test

test_run: $(BIN)/solver_test
	$(BIN)/solver_test -o "../Data/Test/" -n 19  -s 0.0 -e 12.0 -h 0.001 -i "N_SCALING" -v 0.0000001 -d 0.0000001 -a 1.5 -b 1.5 -t "TEST"


# Tetst Case
$(BIN)/solver_test: $(SRCOBJ_Test) $(DEPS)
	$(GCC) $(CCFLAG) $(TEST_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%_test.o: $(SRCDIR)/%.c $(DEPS)
	$(GCC) $(CCFLAG) $(TEST_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

# ---------------------------------------------------------------------
#  Clean up
# ---------------------------------------------------------------------
clean:
	rm -rf $(OBJBIN)
	rm -rf $(BIN)
	rm -rf $(OUT)
# ---------------------------------------------------------------------
#  End of File
# ---------------------------------------------------------------------