/**
* @file solver.c 
* @author Enda Carroll
* @date Jun 2021
* @brief file containing the main functions used in the pseudopectral method
*/
// ---------------------------------------------------------------------
//  Standard Libraries and Headers
// ---------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <complex.h>

// ---------------------------------------------------------------------
//  User Libraries and Headers
// ---------------------------------------------------------------------
#include "data_types.h"
// #include "hdf5_funcs.h"
#include "utils.h"
#include "solver.h"
// #include "sys_msr.h"
// #include "force.h"
// ---------------------------------------------------------------------
//  Global Variables
// ---------------------------------------------------------------------
// Define RK4 variables - Butcher Tableau
#if defined(__RK4)
static const double RK4_C2 = 0.5, 	  RK4_A21 = 0.5, \
				  	RK4_C3 = 0.5,	           					RK4_A32 = 0.5, \
				  	RK4_C4 = 1.0,                      									   RK4_A43 = 1.0, \
				              	 	  RK4_B1 = 1.0/6.0, 		RK4_B2  = 1.0/3.0, 		   RK4_B3  = 1.0/3.0, 		RK4_B4 = 1.0/6.0;
// Define RK5 Dormand Prince variables - Butcher Tableau
#elif defined(__RK5) || defined(__DPRK5)
static const double RK5_C2 = 0.2, 	  RK5_A21 = 0.2, \
				  	RK5_C3 = 0.3,     RK5_A31 = 3.0/40.0,       RK5_A32 = 0.5, \
				  	RK5_C4 = 0.8,     RK5_A41 = 44.0/45.0,      RK5_A42 = -56.0/15.0,	   RK5_A43 = 32.0/9.0, \
				  	RK5_C5 = 8.0/9.0, RK5_A51 = 19372.0/6561.0, RK5_A52 = -25360.0/2187.0, RK5_A53 = 64448.0/6561.0, RK5_A54 = -212.0/729.0, \
				  	RK5_C6 = 1.0,     RK5_A61 = 9017.0/3168.0,  RK5_A62 = -355.0/33.0,     RK5_A63 = 46732.0/5247.0, RK5_A64 = 49.0/176.0,    RK5_A65 = -5103.0/18656.0, \
				  	RK5_C7 = 1.0,     RK5_A71 = 35.0/384.0,								   RK5_A73 = 500.0/1113.0,   RK5_A74 = 125.0/192.0,   RK5_A75 = -2187.0/6784.0,    RK5_A76 = 11.0/84.0, \
				              		  RK5_B1  = 35.0/384.0, 							   RK5_B3  = 500.0/1113.0,   RK5_B4  = 125.0/192.0,   RK5_B5  = -2187.0/6784.0,    RK5_B6  = 11.0/84.0, \
				              		  RK5_Bs1 = 5179.0/57600.0, 						   RK5_Bs3 = 7571.0/16695.0, RK5_Bs4 = 393.0/640.0,   RK5_Bs5 = -92097.0/339200.0, RK5_Bs6 = 187.0/2100.0, RK5_Bs7 = 1.0/40.0;
#endif
// ---------------------------------------------------------------------
//  Function Definitions
// ---------------------------------------------------------------------
/**
 * Main function that performs the solver
 */
void Solve(void) {

	// Initialize variables
	const long int N = 19;

	// Initialize the Runge-Kutta struct
	struct RK_data_struct* RK_data;	   // Initialize pointer to a RK_data_struct
	struct RK_data_struct RK_data_tmp; // Initialize a RK_data_struct
	RK_data = &RK_data_tmp;		       // Point the ptr to this new RK_data_struct

	
	// -------------------------------
	// Allocate memory
	// -------------------------------
	AllocateMemory(N, RK_data);

	// // -------------------------------
	// // FFTW Plans Setup
	// // -------------------------------
	// InitializeFFTWPlans(N);

	// -------------------------------
	// Initialize the System
	// -------------------------------
	// // If in testing / debug mode - create/ open test file
	// #if defined(DEBUG)
	// OpenTestingFile();
	// #endif

	// Initialize the collocation points and wavenumber space 
	InitializeShellWavenumbers(run_data->k, N);

	// Get initial conditions - seed for random number generator is set here
	InitialConditions(run_data->w_hat, run_data->u, run_data->u_hat, N);

	// -------------------------------
	// Integration Variables
	// -------------------------------
	// Initialize integration variables
	double t0;
	double t;
	double dt;
	double T;
	long int trans_steps;
	#if defined(__DPRK5)
	int try = 1;
	double dt_new;
	#endif

	// // Get timestep and other integration variables
	// InitializeIntegrationVariables(&t0, &t, &dt, &T, &trans_steps);

}
/**
 * Function to compute the initial condition for the integration
 */
void InitialConditions(void) {

	// Initialize variables
	
	// ------------------------------------------------
    // Set Seed for RNG
    // ------------------------------------------------
    srand(123456789);

	if(!(strcmp(sys_vars->u0, "N_SCALING"))) {
		// ------------------------------------------------
		// Scaling in N Initial Condition
		// ------------------------------------------------	
	}
	else if(!(strcmp(sys_vars->u0, "RANDOM"))) {
		// ------------------------------------------------
		// Default - Random Initial Conditions
		// ------------------------------------------------	
	}
	else {
		// ------------------------------------------------
		// Default - Pure Random Initial Conditions
		// ------------------------------------------------	
	}
}
/**
 * Function to initialize the shell wavenumber array
 * 
 * @param k Array to contain the shell wavenumber
 * @param N int containging the maximum shell level
 */
void InitializeShellWavenumbers(long int* k, const long int N) {

	// -------------------------------
	// Define Shell Wavenumbers
	// -------------------------------
	for (int i = 0; i < N; ++i) {
		run_data->k[i] = K_0 * pow(LAMBDA, i);
		printf("k[%d]: %ld\n", i, run_data->k[i]);
	}
}
/**
 * Wrapper function used to allocate memory all the nessecary local and global system and integration arrays
 * 
 * @param NBatch  int holding the number of shell levels
 * @param RK_data Pointer to struct containing the integration arrays
 */
void AllocateMemory(const long int N, RK_data_struct* RK_data) {

	// Initialize variables
	

	// -------------------------------
	// Allocate Shell Wavenumbers
	// -------------------------------
	// Allocate the shell wavenumbers
	run_data->k = (long int* )fftw_malloc(sizeof(long int) * N);  // k
	if (run_data->k == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Shell Wavenumber List");
		exit(1);
	}

	// ----------------------------------
	// Allocate Velocity Field Variables
	// ----------------------------------
	// Full velocity field
	run_data->u = (fftw_complex* ) fftw_malloc(sizeof(fftw_complex) * (N + 4));
	if (run_data->u == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Velocity Field");
		exit(1);
	}
	#if defined(PHASE_ONLY)	
	// The Fourier amplitudes
	run_data->a_n = (double* ) fftw_malloc(sizeof(double) * (N + 4));
	if (run_data->a_n == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Velocity Field Amplitude");
		exit(1);
	}
	// The Fourier phases
	run_data->phi_n = (double* ) fftw_malloc(sizeof(double) * (N + 4));
	if (run_data->phi_n == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Velocity Field Phase");
		exit(1);
	}
	#endif

	// ---------------------------------
	// Allocate Magnetic Field Variables
	// ---------------------------------
	#if defined(__MAGNETO)
	// Full velocity field
	run_data->b = (fftw_complex* ) fftw_malloc(sizeof(fftw_complex) * (N + 4));
	if (run_data->b == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Field");
		exit(1);
	}	
	#if defined(PHASE_ONLY)	
	// The Fourier amplitudes
	run_data->b_n = (double* ) fftw_malloc(sizeof(double) * (N + 4));
	if (run_data->b_n == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Field Amplitude");
		exit(1);
	}
	// The Fourier phases
	run_data->psi_n = (double* ) fftw_malloc(sizeof(double) * (N + 4));
	if (run_data->psi_n == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Field Phase");
		exit(1);
	}
	#endif
	#endif
}
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------