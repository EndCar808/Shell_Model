/**
* @file data_types.h 
* @author Enda Carroll
* @date Sept 2022
* @brief file containing the main data types and global variables
*/
// ---------------------------------------------------------------------
//  Standard Libraries and Headers
// ---------------------------------------------------------------------
#ifndef __DATA_TYPES
#include <complex.h>
#ifndef __HDF5_HDR
#include <hdf5.h>
#include <hdf5_hl.h>
#define __HDF5_HDR
#endif
	
#include <gsl/gsl_histogram.h> 
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_math.h>

// ---------------------------------------------------------------------
//  Compile Time Macros and Definitions
// ---------------------------------------------------------------------
#define checkError(x) ({int __val = (x); __val == -1 ? \
	({fprintf(stderr, "ERROR ("__FILE__":%d) -- %s\n", __LINE__, strerror(errno)); \
	exit(-1);-1;}) : __val; })

// For coloured printing to screen
#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define BLUE    "\x1b[34m"
#define MAGENTA "\x1b[35m"
#define CYAN    "\x1b[36m"
#define RESET   "\x1b[0m"
// ---------------------------------------------------------------------
//  Integration Functionality
// ---------------------------------------------------------------------
// These definitions control the functionality of the solver. These definitions
// are also passed at compilation for more control of the solver. Definitions turned
// on in this file means that they WILL be turned on no matter what is passed at compilation

// Choose which system to solver 
// #define __MAGNETO
// Choose which integrator to use
#if defined(__INT_FAC_RK4)
#define INT_FAC_RK4
#endif
#if defined(__RK4)
#define RK4
#endif
#if defined(__AB4CN)
#define AB4CN
#endif
// #define __RK5
// #define __DPRK5
// Choose whether to print updates to screen
#define __PRINT_SCREEN
// Adaptive stepping Indicators 
#define ADAPTIVE_STEP 1          // Indicator for the adaptive stepping
#define CFL_STEP 1 				 // Indicator to use a CFL like condition with the adaptive stepping
// Solver Types 
#define HYPER_VISC 1			 // Turned on hyperviscosity if called for at compilation time
#define VISC_POW 2.0             // The power of the hyperviscosity -> 1.0 means no hyperviscosity
#define HYPO_DIFF 1   			 // Turn on Ekman drag if called for at compilation time
#define HYPO_DIFF_POW -2.0 		 // The power of the Eckman drag term -> 0.0 means no drag
// For allow transient dynamics
#define TRANSIENT_ITERS 1        // Indicator for transient iterations
#define TRANSIENT_FRAC 0.2       // Fraction of total iteration = transient iterations
// Allow for Phase Only mode
#if defined(__PHASE_ONLY_FXD_AMP) // Turn on phase only mode if called for at compilation
#define PHASE_ONLY_FXD_AMP
#endif
#if defined(__PHASE_ONLY) // Turn on phase only direct mode if called for at compilation
#define PHASE_ONLY
#endif
// Computing stats will be decided at compilation time
#if defined(__STATS)
#define STATS
#endif
// Testing the solver will be decided at compilation
#if defined(__TESTING)
#define TESTING
#endif
// If debugging / testing is called for at compile time
#if defined(__DEBUG)
#define DEBUG
#endif
// ---------------------------------------------------------------------
//  Datasets to Write to File
// ---------------------------------------------------------------------
// These definitions control which datasets are to be computed and written to file
// Turning these on in this file means that they WILL be on at compilation time

// Choose whether to save the Fourier Velocity
// #define __VEL
// Choose whether to save the Fourier Magnetic Field
// #define __MAG
// Choose whether to save the Fourier Velocity Phase and Amp
// #define __VEL_AMP
// #define __VEL_PHI
// Choose whether to save the Fourier Magnetic Field Phase and Amp
// #define __MAG_AMP
// #define __MAG_PSI
// Choose whether to save the Nonlinear term or RHS of equation of motion
// #define __RHS
// #define __NONLIN
// Choose whether to save the forcing or not
// define __FORCING
// Choose whether to compute system measures
// #define __SYS_MEASURES
// #define __ENRG_FLUX
// #define __TOT_ENRG_FLUX
// #define __ENRG_FLUX_AVG
// #define __VEL_AMP_AVG
// #define __MAG_AMP_AVG
// Choose whether to compute the spectra
// #define __ENRG_SPECT
// #define __DISS_SPECT
// #define __ENRG_SPECT_AVG
// #define __DISS_SPECT_AVG
// Choose whether to compute phase sync data
// #define __CONSERVED_PHASES
// #define __PHASE_SYNC
// #define __PHASE_SYNC_STATS
// Choose whether to save the time, collocation points and wavenumbers
#define __TIME
#define __WAVELIST
// Choose which stats data to record
#if defined(__STATS)
// #define __VEL_HIST
// #define __MAG_HIST
#define __STR_FUNC_VEL
#define __STR_FUNC_MAG
#define __STR_FUNC_VEL_FLUX
#define __STR_FUNC_MAG_FLUX
#endif
// Writing datasets to file
#define DSET_VEL 0
#define DSET_VEL_AMP 1
#define DSET_VEL_PHI 2
#define DSET_MAG 3
#define DSET_MAG_AMP 4
#define DSET_MAG_PSI 5
#define DSET_ENRG_FLUX 6
#define DSET_ENRG_DISS_VEL 7
#define DSET_ENRG_INPT_VEL 8
#define DSET_ENRG_DISS_MAG 9
#define DSET_ENRG_INPT_MAG 10
#define DSET_FORCING_U 11
#define DSET_FORCING_B 12
#define DSET_ENRG_SPECT 13
#define DSET_DISS_SPECT 14
#define DSET_VEL_TRIADS 15
#define DSET_MAG_TRIADS 16
#define DSET_VEL_PHASE_DIFF 17
#define DSET_MAG_PHASE_DIFF 18
#define DSET_VEL_TRIAD_ORDER 19
#define DSET_MAG_TRIAD_ORDER 20
#define DSET_VEL_PHASE_DIFF_ORDER 21
#define DSET_MAG_PHASE_DIFF_ORDER 22
#define NUM_DSETS 23
// ---------------------------------------------------------------------
//  Global Variables
// ---------------------------------------------------------------------
// These definitions define some of the solver parameters.
#define TRANS_FRAC 0.2          		// Fraction of time to ignore before saving to file
// Define the shell wavenumber parameters
#define K_0 1.0 						// The shell wavenumber prefactor
#define LAMBDA 2.0              		// The intershell ratio for the shell wavenumber 
// Define the equations of motion parameters
#define EPSILON 0.5 					// Interaction coefficient for the HD shell model
#define EPSILON_M 1.0/3.0 				// Interaction coefficient for the MHD shell model
// Define forcing parameters
#define FORC_STOC_SIGMA 0.1				// The ampitude of the Gaussian white noise 
// System checking parameters
#define MIN_STEP_SIZE 1e-10 			// The minimum allowed stepsize for the solver 
#define MAX_ITERS 1e+20					// The maximum iterations to perform
#define MAX_FIELD_LIM 1e+100    		// The maximum allowed velocity &/or magnetic
// Stats parameters
#define NUM_POW 6 						// The highest moment to compute for the structure function
#define NUM_RUN_STATS 7 				// The number of running stats moments to record
#define VEL_BIN_LIM	40					// The bin limits (in units of standard deviations) for the velocity histogram
#define VEL_NUM_BINS 1000				// The number of bins to use for the velocity histograms
// Initial Condtion parameters
#define NO_INPUT_FILE 0 				// Indicates if system is to staart by generating initial conditions
#define INPUT_FILE 1 					// Indicates if system is to staart from an input file
// Phase sync parameters
#define NUM_MAG_TRIAD_TYPES 3 			// The number of extra triad types under the magnetic field
#define NUM_PHASE_SYNC_HIST_BINS 1000	// The number of histogram bins for the phase sync stats data
// ---------------------------------------------------------------------
//  Global Struct Definitions
// ---------------------------------------------------------------------
// System variables struct
typedef struct system_vars_struct {
	char u0[64];						// String to indicate the initial condition to use
	char forcing[64];					// String to indicate what type of forcing is selected
	long int N;							// Array holding the no. of collocation shells
	long int num_t_steps;				// Number of iteration steps to perform
	long int num_print_steps;           // Number of times system was saved to file
	long int tot_iters;					// Records the total executed iterations
	long int tot_save_steps;			// Records the total saving iterations
	long int trans_iters;				// The number of transients iterations to perform before printing to file
	double t0;							// Intial time
	double T;							// Final time
	double t;							// Time variable
	double dt;							// Timestep
	double min_dt;						// Smallest timestep achieved when adaptive stepping is on
	double max_dt;						// Largest timestep achieved when adaptive stepping is on
	int print_every;                    // Records how many iterations are performed before printing to file
	long int SAVE_EVERY;				// For specifying how often to print
	double force_scale_var;				// The scaling variable for the forced modes
	int force_k; 						// The forcing wavenumber 
	int local_forcing_proc;				// Identifier used to indicate which process contains modes that need to be forced
	int num_forced_modes; 				// The number of modes to be forced on the local process
	int TRANS_ITERS_FLAG;				// Flag for indicating whether transient iterations are to be performed
	double TRANS_ITERS_FRAC;			// The fraction of total iterations that will be ignored
	int ADAPT_STEP_FLAG;			 	// Flag for indicating if adaptive stepping is to be used
	double CFL_CONST;					// The CFL constant for the adaptive step
	int CFL_COND_FLAG;					// Flag for indicating if the CFL like condition is to be used for the adaptive stepping
	double ALPHA;						// Slope of the velocity energy spectrum = 2 * ALPHA
	double BETA; 						// Slope of the magnetic energy spectrum = 2 * BETA
	double NU;							// The viscosity
	int HYPER_VISC_FLAG;				// Flag to indicate if hyperviscosity is to be used
	double HYPER_VISC_POW;				// The power of the hyper viscosity to use
	double ETA;							// The magnetic diffusivity
	int HYPO_MAG_DIFF_FLAG;				// Flag for indicating if ekman drag is to be used
	double HYPO_MAG_DIFF_POW;			// The power of the hyper drag to be used
	double EPS; 						// Defines the interaction coefficient for the HD shell model equation 
	double EPS_M; 						// Defines the interaction coefficient for the MHD shell model equation
	double k_0;							// Defines the wavenumber constant k_0
	double Lambda;						// Defines the intershell ratio for radius of the shell wavenumber
	int INPUT_FILE_FLAG;				// Flag to indicate if solver is to start by reading in initial condition from file
	double eddy_turnover_time;			// Eddy turnover time is going to be defined as the max eddy turnover time during transient iterations
	double trans_time;					// Integration time it takes to perform the transient iterations
} system_vars_struct;

// Runtime data struct
typedef struct runtime_data_struct {
	double* k;			  				// Array to hold wavenumbers
	double complex* u;	      			// Fourier space velocity
	double complex* b;	      			// Fourier space vorticity
	double complex* rhs; 		  		// Array to hold the RHS of the equation of motion
	double complex* nonlinterm; 		// Array to hold the nonlinear term
	double* a_n;			  			// Fourier vorticity amplitudes
	double* a_n_t_avg;		  			// Fourier vorticity amplitudes
	double* tmp_a_n;		  			// Array to hold the amplitudes of the fourier vorticity before marching forward in time
	double* phi_n;			  			// Fourier velocity phases
	double* b_n;			  			// Fourier velocity amplitudes
	double* b_n_t_avg;		  			// Fourier velocity amplitudes
	double* tmp_b_n;		  			// Array to hold the amplitudes of the fourier vorticity before marching forward in time
	double* psi_n;						// Array to hold the phases of the magnetic modes
	double* time;			  			// Array to hold the simulation times
	double* tot_energy;       			// Array to hold the total energy over the simulation
	double* tot_hel_u;		  			// Array to hold the total helicity in the velocity field
	double* tot_hel_b;		  			// Array to hold the total helicity in the magnetic field
	double* tot_cross_hel;	  			// Array to hold the total cross helicity
	double* tot_diss_u;					// Array to hold the total velocity dissipation 
	double* tot_diss_b;					// Array to hold the total magnetic dissipation 
	double* energy_flux;				// Array to hold the energy flux
	double* energy_diss_u;				// Array to hold the energy dissipation for the velocity field
	double* energy_diss_b;				// Array to hold the energy dissipation for the magnetic field
	double* energy_input_u;				// Array to hold the energy input for the velocity field in the flux balance computation
	double* energy_input_b;				// Array to hold the energy input for the magnetic field in the flux balance computation
	double* tot_energy_flux;			// Array to hold the total energy flux for the current iteration
	double* tot_energy_diss_u;			// Array to hold the total energy dissipation for the velocity field for the current iteration
	double* tot_energy_diss_b;			// Array to hold the total energy dissipation for the magnetic field for the current iteration
	double* tot_energy_input_u;			// Array to hold the total energy input for the velocity field in the flux balance computation for the current iteration
	double* tot_energy_input_b;			// Array to hold the total energy input for the magnetic field in the flux balance computation for the current iteration
	double* energy_spect;				// Array containing the energy spectrum for the current iteration 
	double* diss_spect;					// Array containing the dissipation spectrum for the current iteration 
	double* energy_flux_t_avg;			// Array containing the energy flux time average over simulation
	double* energy_spect_t_avg;			// Array containing the energy spectrum time average over simulation
	double* diss_spect_t_avg;			// Array containing the dissipation spectrum time average over simulation
	double* u_charact;					// Array containing the characteristic velocity
	double* int_scale;					// Array containing the characteristic length scale
	double* taylor_micro_scale;			// Array containing the Taylor microscale
	double* reynolds_no;				// Array containing the Reynolds no.
	double* kolmogorov_scale;			// Array containing the Kolmogorov length scale
	double complex* forcing_u;  		// Array to hold the forcing for the current timestep for the velocity field
	double complex* forcing_b;  		// Array to hold the forcing for the current timestep for the magnetic field
	double* forcing_scaling;  			// Array to hold the initial scaling for the forced modes
	int* forcing_indx;		  			// Array to hold the indices of the forced modes
	int* forcing_k;			  			// Array containg the wavenumbers for the forced modes
	long int num_sys_msr_steps;			// Counts the number of times system measures have been computed, for time averaging.
} runtime_data_struct;

// Runge-Kutta Integration struct
typedef struct RK_data_struct {
	// #if defined(PHASE_ONLY)
	// double* RK1_u;		  		 	// Array to hold the result of the first stage for the velocity field
	// double* RK2_u;		  		 	// Array to hold the result of the second stage for the velocity field
	// double* RK3_u;		  		 	// Array to hold the result of the third stage for the velocity field
	// double* RK4_u;		  		 	// Array to hold the result of the fourth stage for the velocity field
	// double* RK1_b;		  		 	// Array to hold the result of the first stage for the magnetic field
	// double* RK2_b;		  		 	// Array to hold the result of the second stage for the magnetic field
	// double* RK3_b;		  		 	// Array to hold the result of the third stage for the magnetic field
	// double* RK4_b;		  		 	// Array to hold the result of the fourth stage for the magnetic field
	// double* RK_u_tmp;		     	// Array to hold the tempory updates to u - input to Nonlinear term function
	// double* RK_b_tmp;		     	// Array to hold the tempory updates to b - input to Nonlinear term function
	// double* AB_tmp_u;	  			// Array to hold the result of the RHS/Nonlinear term for the Adams Bashforth scheme
	// double* AB_tmp_nonlin_u[3];	 	// Array to hold the previous 3 RHS/Nonlinear terms to update the Adams Bashforth scheme
	// double* AB_tmp_b;	  			// Array to hold the result of the RHS/Nonlinear term for the Adams Bashforth scheme
	// double* AB_tmp_nonlin_b[3];	 	// Array to hold the previous 3 RHS/Nonlinear terms to update the Adams Bashforth scheme
	// #else
	double complex* RK1_u;		  		// Array to hold the result of the first stage for the velocity field
	double complex* RK2_u;		  		// Array to hold the result of the second stage for the velocity field
	double complex* RK3_u;		  		// Array to hold the result of the third stage for the velocity field
	double complex* RK4_u;		  		// Array to hold the result of the fourth stage for the velocity field
	double complex* RK1_b;		  		// Array to hold the result of the first stage for the magnetic field
	double complex* RK2_b;		  		// Array to hold the result of the second stage for the magnetic field
	double complex* RK3_b;		  		// Array to hold the result of the third stage for the magnetic field
	double complex* RK4_b;		  		// Array to hold the result of the fourth stage for the magnetic field
	double complex* RK_u_tmp;		  		// Array to hold the tempory updates to u - input to Nonlinear term function
	double complex* RK_b_tmp;		  		// Array to hold the tempory updates to b - input to Nonlinear term function
	double complex* AB_tmp_u;	  			// Array to hold the result of the RHS/Nonlinear term for the Adams Bashforth scheme
	double complex* AB_tmp_nonlin_u[3];	// Array to hold the previous 3 RHS/Nonlinear terms to update the Adams Bashforth scheme
	double complex* AB_tmp_b;	  			// Array to hold the result of the RHS/Nonlinear term for the Adams Bashforth scheme
	double complex* AB_tmp_nonlin_b[3];	// Array to hold the previous 3 RHS/Nonlinear terms to update the Adams Bashforth scheme
	// #endif
	int AB_pre_steps;				// The number of derivatives to perform for the AB4 scheme
} RK_data_struct;

// HDF5 file info struct
typedef struct HDF_file_info_struct {
	char input_file_name[512];		 // Array holding input file name
	char output_file_name[512];      // Output file name array
	char stats_file_name[512];       // Stats Output file name array
	char phase_sync_file_name[512];  // Phase Sync Output file name array
	char system_msr_file_name[512];  // Output file name array
	char input_dir[512];			 // Inputs directory
	char output_dir[512];			 // Output directory
	char output_tag[64]; 			 // Tag to be added to the output directory
	hid_t input_file_handle;		 // File handle for the input file
	hid_t output_file_handle;		 // Main file handle for the output file 
	hid_t stats_file_handle;		 // Stats file handle for the output file 
	hid_t phase_sync_file_handle;	 // Phase Sync file handle for the output file 
	hid_t sys_msr_file_handle;		 // System Measures file handle for the output file 
	hid_t COMPLEX_DTYPE;			 // Complex datatype handle
	int file_only;					 // Indicates if output should be file only with no output folder created
	hid_t test_file_handle;          // File handle for testing
	char test_file_name[512];        // File name for testing
	hid_t file_space[NUM_DSETS];	 // Identifier for File memory space for the slabbed datasets
	hid_t data_set[NUM_DSETS];  	 // Identifier for Size of the data memory space for the slabbed datasets
	hid_t mem_space[NUM_DSETS]; 	 // Identifier for Size of memory for the slabbed datasets
} HDF_file_info_struct;

// Phase Sync Data struct
typedef struct phase_sync_data_struct {
	double* triads_u;						// Array to hold the time series of each triad for the velocity field
	double* triads_b;						// Array to hold the time series of each triad for the magnetic field
	double* phase_diff_u;					// Array to hold the time series of each phase difference for the velocity field
	double* phase_diff_b;					// Array to hold the time series of each phase difference for the magnetic field
	double complex* phase_diff_u_order;		// Array to hold the time series of the Kuramoto order parameter of the phase difference for the velocity field
	double complex* triad_u_order;			// Array to hold the time series of the Kuramoto order parameter of the triads for the velocity field
	double complex* phase_diff_b_order;     // Array to hold the time series of the Kuramoto order parameter of the phase difference for the magnetic field
	double complex* triad_b_order;			// Array to hold the time series of the Kuramoto order parameter of the triads for the magnetic field
	long int num_phase_sync_steps;	 		// Counter to record the number of time phase sync has been measures, for time averaging
	int num_triads; 						// Records the number of triads in the system
	int num_phase_diff; 					// Records the number of phase difference phases in the system
	gsl_histogram** triad_u_hist;			// Struct to hold the histogram info for the velocity field triads over time
	gsl_histogram** triad_b_hist;			// Struct to hold the histogram info for the magnetic field triads over time
	gsl_histogram** phase_diff_u_hist;		// Struct to hold the histogram info for the velocity field phase differences over time
	gsl_histogram** phase_diff_b_hist;		// Struct to hold the histogram info for the magnetic field phase differences over time
	gsl_histogram** phase_order_Phi_u_hist;	// Struct to hold the histogram info for the velocity field phase orders Phi over time
	gsl_histogram** phase_order_Phi_b_hist;	// Struct to hold the histogram info for the magnetic field phase orders Phi over time
} phase_sync_data_struct;

// Stats data struct
typedef struct stats_data_struct {
	double* vel_str_func[NUM_POW];					// Array to hold the structure functions of the Velocity modes
	double* mag_str_func[NUM_POW];			  		// Array to hold the structure functions of the magnetic modes
	double* vel_flux_str_func[2][NUM_POW];			// Array to hold the structure functions of the value of flux Velocity modes
	double* mag_flux_str_func[2][NUM_POW];			// Array to hold the structure functions of the value of flux magnetic modes
	double* vel_flux_str_func_abs[2][NUM_POW];		// Array to hold the structure functions of the absolute value of flux Velocity modes
	double* mag_flux_str_func_abs[2][NUM_POW];		// Array to hold the structure functions of the absolute value of flux magnetic modes
	gsl_rstat_workspace** real_vel_moments;			// Struct to hold the running stats for the velocity field
	gsl_rstat_workspace** real_mag_moments;			// Struct to hold the running stats for the magnetic field
	gsl_histogram** real_vel_hist;					// Struct to hold the histogram info for the velocity field
	gsl_histogram** real_mag_hist;					// Struct to hold the histogram info for the magnetic field
	long int num_stats_steps;						// Counter for the number of steps statistics have been computed
	int set_stats_flag;								// Flag to indicate if the stats objects such as running stats and histogram limits need to be set
} stats_data_struct;


// Complex datatype struct for HDF5
typedef struct complex_type_tmp {
	double re;   			 // real part 
	double im;   			 // imaginary part 
} complex_type_tmp;


// Declare the global variable pointers across all files
extern system_vars_struct *sys_vars; 		    // Global pointer to system parameters struct
extern runtime_data_struct *run_data; 			// Global pointer to system runtime variables struct 
extern HDF_file_info_struct *file_info; 		// Global pointer to system forcing variables struct 
extern stats_data_struct *stats_data;           // Globale pointer to the statistics struct
extern phase_sync_data_struct *phase_sync;      // Globale pointer to the phase sync data struct

#define __DATA_TYPES
#endif
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------
