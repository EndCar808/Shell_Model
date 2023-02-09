# Enda Carroll
# Sept 2022
# Makfile to compile programmes for solver for the Magnetohydrodynamic GOY shell Model

 
# ---------------------------------------------------------------------
#  Compiler flags
# ---------------------------------------------------------------------
# Get the operating system:
# Need to strip leading and/or trailing whitespaces otherwise if statement wont work
OS        = $(strip $(shell lsb_release -si))
MACH_NAME = $(strip $(shell hostname))

# Set the operating system dependent compiler variables
ifeq ($(OS), Ubuntu)
	ifeq ($(MACH_NAME), DESKTOP-SS58OJO)
		# LAPTOPL: WSL 2
		GCC := gcc
		CCFLAG := -O3 -W -Wall -g
		LIBFLAG := -lm -lhdf5 -lhdf5_hl
		GSLFLAG := -lgsl -lgslcblas
		# LAPACKFLAG := -llapacke -llapack -lblas -lgfortran
		INCDIR = -I/home/endacarroll/Downloads/hdf5-1.12.2/hdf5/include -I/usr/include/gsl/
		LIBDIR = -L/home/endacarroll/Downloads/hdf5-1.12.2/hdf5/lib
	else ifeq ($(MACH_NAME), turbo)
		# TURBO
		GCC := gcc
		CCFLAG := -O3 -W -Wall -g
		LIBFLAG := -lm -lhdf5 -lhdf5_hl
		GSLFLAG := -lgsl -lgslcblas
		# LAPACKFLAG := -llapacke -llapack -lblas -lgfortran
		INCDIR = -I/opt/software/hdf5-1.12.2-serial/include/ -I/usr/include/gsl/
		LIBDIR = -L/opt/software/hdf5-1.12.2-serial/lib/
	else
		# CHIRP
		GCC := gcc
		CCFLAG := -O3 -W -Wall -g
		LIBFLAG := -lm -lhdf5 -lhdf5_hl
		GSLFLAG := -lgsl -lgslcblas
		# LAPACKFLAG := -llapacke -llapack -lblas -lgfortran
		INCDIR = -I/usr/include/hdf5/serial/  -I/usr/include/gsl/ 
		LIBDIR = -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ 
	endif
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
# Flags to be passed to the hyrdo solver for normal runs
HYDRO_SOLVER_IFRK4_FLAGS := -D__INT_FAC_RK4 -D__VEL # -D__ADAPTIVE_STEP -D__CFL_STEP
HYDRO_SOLVER_RK4_FLAGS := -D__RK4 -D__VEL # -D__ADAPTIVE_STEP -D__CFL_STEP
HYDRO_SOLVER_AB4CN_FLAGS := -D__AB4CN -D__VEL # -D__ADAPTIVE_STEP -D__CFL_STEP

# Flags to be passed to the mag hyrdo solver for normal runs
MAG_HYDRO_SOLVER_IFRK4_FLAGS := -D__INT_FAC_RK4 -D__MAGNETO -D__VEL -D__MAG # -D__ADAPTIVE_STEP -D__CFL_STEP
MAG_HYDRO_SOLVER_RK4_FLAGS := -D__RK4 -D__MAGNETO -D__VEL -D__MAG # -D__ADAPTIVE_STEP -D__CFL_STEP
MAG_HYDRO_SOLVER_AB4CN_FLAGS := -D__AB4CN -D__MAGNETO -D__VEL -D__MAG # -D__ADAPTIVE_STEP -D__CFL_STEP

# Flags to be passed to the mag hyrdo solver for normal runs
ELSASSAR_MHD_SOLVER_IFRK4_FLAGS := -D__INT_FAC_RK4 -D__ELSASSAR_MHD -D__Z_PLUS -D__Z_MINUS # -D__ADAPTIVE_STEP -D__CFL_STEP
ELSASSAR_MHD_SOLVER_RK4_FLAGS := -D__RK4 -D__ELSASSAR_MHD -D__Z_PLUS -D__Z_MINUS # -D__ADAPTIVE_STEP -D__CFL_STEP
ELSASSAR_MHD_SOLVER_AB4CN_FLAGS := -D__AB4CN -D__ELSASSAR_MHD -D__Z_PLUS -D__Z_MINUS # -D__ADAPTIVE_STEP -D__CFL_STEP

# Phase only flags to be passed to solver
PHASE_ONLY_FLAGS := -D__PHASE_ONLY -D__RK4 -D__VEL_AMP -D__VEL_PHI #-D__MAGNETO -D__HYPER -D__TRANSIENTS

# Phase only flags to be passed to solver
PHASE_ONLY_FXD_AMP_FLAGS := -D__PHASE_ONLY_FXD_AMP -D__RK4 -D__VEL_AMP -D__VEL_PHI #-D__MAGNETO -D__HYPER -D__TRANSIENTS

# Dataset Flags to be passed to the solver
DSET_FLAGS = -D__SYS_MEASURES 
DSET_FLAGS += -D__CONSERVED_PHASES -D__PHASE_SYNC -D__PHASE_SYNC_STATS 
DSET_FLAGS += -D__STATS 
DSET_FLAGS += -D__ENRG_SPECT_AVG -D__DISS_SPECT_AVG -D__TOT_ENRG_FLUX -D__TOT_KIN_HEL_FLUX -D__ENRG_FLUX_AVG -D__VEL_AMP_AVG -D__MAG_AMP_AVG -D__KIN_HEL_FLUX_AVG -D__KIN_HEL_FLUX
DSET_FLAGS += -D__KIN_ENRG_SPECT_AVG -D__MAG_ENRG_SPECT_AVG -D__PSEUDO_ENRG_FLUX_AVG
DSET_FLAGS += -D__ENRG_SPECT -D__DISS_SPECT -D__ENRG_FLUX -D__KIN_ENRG_SPECT -D__MAG_ENRG_SPECT
DSET_FLAGS += -D__FORCING
DSET_FLAGS += -D__VEL 
DSET_FLAGS += -D__MAG
 
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
SRCOBJ_MHD_IF  	:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_mag_hydro_IF.o, $(SRCFILE))
SRCOBJ_MHD_RK  	:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_mag_hydro_RK.o, $(SRCFILE))
SRCOBJ_MHD_AB  	:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_mag_hydro_AB.o, $(SRCFILE))
SRCOBJ_ELS_IF  	:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_elsassar_mhd_IF.o, $(SRCFILE))
SRCOBJ_ELS_RK  	:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_elsassar_mhd_RK.o, $(SRCFILE))
SRCOBJ_ELS_AB  	:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_elsassar_mhd_AB.o, $(SRCFILE))
SRCOBJ_PO  		:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_phase_only.o, $(SRCFILE))
SRCOBJ_PO_D 	:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_phase_only_fxd_amp.o, $(SRCFILE))
SRCOBJ_Test 	:= $(patsubst $(SRCDIR)/%.c, $(OBJBIN)/%_test.o, $(SRCFILE))

DEPS = $(wildcard *.h) 

# ---------------------------------------------------------------------
#  Builds
# ---------------------------------------------------------------------
############
.PHONY: all
############


all: $(BIN)/solver_hydro_IF $(BIN)/solver_hydro_RK $(BIN)/solver_hydro_AB $(BIN)/solver_mag_hydro_IF $(BIN)/solver_mag_hydro_RK $(BIN)/solver_mag_hydro_AB $(BIN)/solver_elsassar_mhd_IF $(BIN)/solver_elsassar_mhd_RK $(BIN)/solver_elsassar_mhd_AB $(BIN)/solver_phase_only $(BIN)/solver_phase_only_fxd_amp $(BIN)/solver_test


####--------
#  Solvers
####--------
##----------------------- Hydro Solvers
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


##----------------------- Mag Hydro Solvers
solver_mag_hydro_IF: $(BIN)/solver_mag_hydro_IF 

$(BIN)/solver_mag_hydro_IF: $(SRCOBJ_MHD_IF) $(DEPS)
	$(GCC) $(CCFLAG) $(MAG_HYDRO_SOLVER_IFRK4_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%_mag_hydro_IF.o: $(SRCDIR)/%.c $(DEPS)
	$(GCC) $(CCFLAG) $(MAG_HYDRO_SOLVER_IFRK4_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

solver_mag_hydro_RK: $(BIN)/solver_mag_hydro_RK

$(BIN)/solver_mag_hydro_RK: $(SRCOBJ_MHD_RK) $(DEPS)
	$(GCC) $(CCFLAG) $(MAG_HYDRO_SOLVER_RK4_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%_mag_hydro_RK.o: $(SRCDIR)/%.c $(DEPS)
	$(GCC) $(CCFLAG) $(MAG_HYDRO_SOLVER_RK4_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

solver_mag_hydro_AB: $(BIN)/solver_mag_hydro_AB

$(BIN)/solver_mag_hydro_AB: $(SRCOBJ_MHD_AB) $(DEPS)
	$(GCC) $(CCFLAG) $(HYDRO_SOLVER_AB4CN_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%_mag_hydro_AB.o: $(SRCDIR)/%.c $(DEPS)
	$(GCC) $(CCFLAG) $(HYDRO_SOLVER_AB4CN_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)


##----------------------- Elsassar Solvers
solver_elsassar_mhd_IF: $(BIN)/solver_elsassar_mhd_IF 

$(BIN)/solver_elsassar_mhd_IF: $(SRCOBJ_ELS_IF) $(DEPS)
	$(GCC) $(CCFLAG) $(ELSASSAR_MHD_SOLVER_IFRK4_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%_elsassar_mhd_IF.o: $(SRCDIR)/%.c $(DEPS)
	$(GCC) $(CCFLAG) $(ELSASSAR_MHD_SOLVER_IFRK4_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

solver_elsassar_mhd_RK: $(BIN)/solver_elsassar_mhd_RK 

$(BIN)/solver_elsassar_mhd_RK: $(SRCOBJ_ELS_RK) $(DEPS)
	$(GCC) $(CCFLAG) $(ELSASSAR_MHD_SOLVER_RK4_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%_elsassar_mhd_RK.o: $(SRCDIR)/%.c $(DEPS)
	$(GCC) $(CCFLAG) $(ELSASSAR_MHD_SOLVER_RK4_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

solver_elsassar_mhd_AB: $(BIN)/solver_elsassar_mhd_AB 

$(BIN)/solver_elsassar_mhd_AB: $(SRCOBJ_ELS_AB) $(DEPS)
	$(GCC) $(CCFLAG) $(ELSASSAR_MHD_SOLVER_AB4CN_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%_elsassar_mhd_AB.o: $(SRCDIR)/%.c $(DEPS)
	$(GCC) $(CCFLAG) $(ELSASSAR_MHD_SOLVER_AB4CN_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)


####-----------
#  Phase Only 
####-----------
solver_phase_only: $(BIN)/solver_phase_only 

$(BIN)/solver_phase_only: $(SRCOBJ_PO) $(DEPS) 
	$(GCC) $(CCFLAG) $(PHASE_ONLY_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%_phase_only.o: $(SRCDIR)/%.c $(DEPS) 
	$(GCC) $(CCFLAG) $(PHASE_ONLY_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

solver_phase_only_fxd_amp: $(BIN)/solver_phase_only_fxd_amp 

$(BIN)/solver_phase_only_fxd_amp: $(SRCOBJ_PO_D) $(DEPS) 
	$(GCC) $(CCFLAG) $(PHASE_ONLY_FXD_AMP_FLAGS) $(DSET_FLAGS) -o $@ $^ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

$(OBJBIN)/%_phase_only_fxd_amp.o: $(SRCDIR)/%.c $(DEPS) 
	$(GCC) $(CCFLAG) $(PHASE_ONLY_FXD_AMP_FLAGS) $(DSET_FLAGS) -c $^ -o $@ ${INCDIR} ${LIBDIR} $(LIBFLAG) $(GSLFLAG)

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