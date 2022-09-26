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
	LIBFLAG := -lfftw3 -lm -lhdf5 -lhdf5_hl #-fopenmp -lpthread -lfftw3_omp
	# GSLFLAG := -lgsl -lgslcblas
	# LAPACKFLAG := -llapacke -llapack -lblas -lgfortran
	INCDIR = -I/usr/include/hdf5/serial/  -I/usr/include/gsl/ 
	LIBDIR = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ 
else ifeq ($(OS), RedHatEnterpriseServer)
	# SONIC
	GCC := gcc
	CCFLAG := -O3 -W -Wall -g
	LIBFLAG := -lfftw3 -lfftw3 -lm -lhdf5 -lhdf5_hl #-fopenmp -lpthread -lfftw3_omp
	INCDIR = 
	LIBDIR = 
else ifeq ($(OS), CentOS)
	# KAY
	GCC := gcc
	CCFLAG := -O3 -W -Wall -g
	LIBFLAG := -lfftw3 -lfftw3 -lm -lhdf5 -lhdf5_hl #-fopenmp -lpthread -lfftw3_omp
	INCDIR = 
	LIBDIR = 
else
	$(error Unknown OS found, please check OS and add appropriate flags in Makefile)
endif

# ---------------------------------------------------------------------
#  Executable flags
# ---------------------------------------------------------------------
# Flags to be passed to the solver for normal runs
SOLVER_FLAGS := -D__MAGNETO -D__INT_FAC_RK4 #-D__HYPER -D__ADAPTIVE_STEP -D__CFL_STEP #-D__TRANSIENTS  -D__EKMN_DRAG

# Phase only flags to be passed to solver
PHASE_ONLY_FLAGS := -D__MAGNETO -D__PHASE_ONLY -D__RK4 # -D__HYPER -D__TRANSIENTS

# Phase only flags to be passed to solver
PHASE_ONLY_DIRECT_FLAGS := -D__MAGNETO -D__PHASE_ONLY_DIRECT -D__RK4 # -D__HYPER -D__TRANSIENTS

# Dataset Flags to be passed to the solver
DSET_FLAGS := -D__SYS_MEASURES -D__ENRG_FLUX
 
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

SRCDIR      := src
SRCFILE     := $(wildcard $(SRCDIR)/*.c)
SRCOBJ      := $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%.o, $(SRCFILE))
SRCOBJ_PO  	:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_phase_only.o, $(SRCFILE))
SRCOBJ_PO_D := $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_phase_only_direct.o, $(SRCFILE))
SRCOBJ_Test := $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_test.o, $(SRCFILE))

DEPS = $(wildcard *.h) 

# ---------------------------------------------------------------------
#  Builds
# ---------------------------------------------------------------------
############
.PHONY: all
############


all: $(BIN)/solver $(BIN)/solver_phase_only $(BIN)/solver_phase_only_direct $(BIN)/solver_test


####--------
#  Solver
####--------
solver: $(BIN)/solver 

$(BIN)/solver: $(SRCOBJ) $(DEPS)
	$(GCC) $(CCFLAG) $(SOLVER_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%.o: $(SRCDIR)/%.c $(DEPS)
	$(GCC) $(CCFLAG) $(SOLVER_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

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