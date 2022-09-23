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
#include "hdf5_funcs.h"
#include "utils.h"
#include "solver.h"
#include "sys_msr.h"
// #include "force.h"
// ---------------------------------------------------------------------
//  Global Variables
// ---------------------------------------------------------------------
// Define RK4 variables - Butcher Tableau
#if defined(__RK4) || defined(__INT_FAC_RK4)
static const double RK4_C2 = 0.5, 	  RK4_A21 = 0.5, \
				  	RK4_C3 = 0.5,	           					RK4_A32 = 0.5, \
				  	RK4_C4 = 1.0,                      									   RK4_A43 = 1.0, \
				              	 	  RK4_B1 = 1.0/6.0, 		RK4_B2  = 1.0/3.0, 		   RK4_B3  = 1.0/3.0, 		RK4_B4 = 1.0/6.0;
#endif
// ---------------------------------------------------------------------
//  Function Definitions
// ---------------------------------------------------------------------
/**
 * Main function that performs the solver
 */
void Solve(void) {

	// Initialize variables
	const long int N = sys_vars->N;

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
	InitialConditions(N);

	// -------------------------------
	// Integration Variables
	// -------------------------------
	// Initialize integration variables
	double t0;
	double t;
	double dt;
	double T;
	long int trans_steps;

	// // Get timestep and other integration variables
	InitializeIntegrationVariables(&t0, &t, &dt, &T, &trans_steps);
	
	// -------------------------------
	// Create & Open Output File
	// -------------------------------
	// Inialize system measurables
	InitializeSystemMeasurables(RK_data);	

	// Create and open the output file - also write initial conditions to file
	CreateOutputFilesWriteICs(N);

	// Print update of the initial conditions to the terminal
	PrintUpdateToTerminal(0, t0, dt, T, 0);

	//////////////////////////////
	// Begin Integration
	//////////////////////////////
	t         += dt;
	int iters = 1;
	int save_data_indx;
	if (sys_vars->TRANS_ITERS_FLAG == TRANSIENT_ITERS) {
		save_data_indx = 0;
	}
	else {
		save_data_indx = 1;
	}
	while (t <= T) {

		// -------------------------------	
		// Integration Step
		// -------------------------------
		#if defined(__INT_FAC_RK4)
		IntFacRK4Step(dt, N, RK_data);
		#endif

		// -------------------------------
		// Write To File
		// -------------------------------
		if (iters % sys_vars->SAVE_EVERY == 0) {

			// Record System Measurables
			ComputeSystemMeasurables(t, save_data_indx, RK_data);

			// If and when transient steps are complete write to file
			if (iters > trans_steps) {
				// Write the appropriate datasets to file 
				WriteDataToFile(t, dt, save_data_indx);
				
				// Update saving data index
				save_data_indx++;
			}
		}
		// -------------------------------
		// Print Update To Screen
		// -------------------------------
		#if defined(__PRINT_SCREEN)
		if (sys_vars->TRANS_ITERS_FLAG == TRANSIENT_ITERS) {
			// Print update that transient iters have been complete
			if (iters == trans_steps) {
				printf("\n\n...Transient Iterations Complete!\n\n");
			}
		}
		if (iters % sys_vars->SAVE_EVERY == 0) {
			// Print update of the system to the terminal 
			if (iters <= sys_vars->trans_iters) {
				// If performing transient iterations the system measures are stored in the 0th index
				PrintUpdateToTerminal(iters, t, dt, T, 0);
			}
			else {
				// Print update of the non transient iterations to the terminal 
				PrintUpdateToTerminal(iters, t, dt, T, save_data_indx - 1);
			}
		}
		#endif


		// -------------------------------
		// Update & System Check
		// -------------------------------
		// Update timestep & iteration counter
		iters++;
		if (sys_vars->ADAPT_STEP_FLAG == ADAPTIVE_STEP) {
			// GetTimestep(&dt);
			t += dt; 
		}
		else {
			#if defined(__DPRK5)
			t += dt;
			#else
			t = iters * dt;
			#endif
		}

		// Check System: Determine if system has blown up or integration limits reached
		SystemCheck(dt, iters);
	}
	//////////////////////////////
	// End Integration
	//////////////////////////////
	
	// ------------------------------- 
	// Final Writes to Output File
	// -------------------------------
	FinalWriteAndCloseOutputFile(N, iters, save_data_indx);

	// -------------------------------
	// Clean Up 
	// -------------------------------
	FreeMemory(RK_data);
}
/**
 * Function to perform one step using the 4th order Runge-Kutta method
 * @param dt       The current timestep of the system
 * @param N        int defining the number of shells
 * @param RK_data  Struct pointing the Integration variables: stages, tmp arrays, rhs and arrays needed for NonlinearRHS function
 */
#if defined(__INT_FAC_RK4)
void IntFacRK4Step(const double dt, const long int N, RK_data_struct* RK_data) {

	// Initialize vairables
	int n;
	fftw_complex int_fac_u;
	fftw_complex int_fac_b;
	#if defined(PHASE_ONLY)
	double po_norm_fac_u;
	#if defined(__MAGNETO)
	double po_norm_fac_b;
	#endif
	#endif

	/////////////////////
	/// RK STAGES
	/////////////////////
	// ----------------------- Stage 1
	NonlinearTerm(run_data->u, run_data->b, RK_data->RK1_u, RK_data->RK1_b, N);
	for (int i = 0; i < N; ++i) {
		// Get tmp index
		n = i + 2;

		// Get the integrating factor
		int_fac_u = exp(-sys_vars->NU * dt * run_data->k[i] * run_data->k[i] / 2.0);
		int_fac_b = exp(-sys_vars->ETA * dt * run_data->k[i] * run_data->k[i] / 2.0);

		// Update temorary velocity term
		RK_data->RK_u_tmp[n] = run_data->u[n] * int_fac_u + dt * RK4_A21 * RK_data->RK1_u[n] * int_fac_u;
		#if defined(__MAGNETO)
		// Update temporary magnetic term
		RK_data->RK_b_tmp[n] = run_data->b[n] * int_fac_b + dt * RK4_A21 * RK_data->RK1_b[n] * int_fac_b;		
		#endif
	}

	// ----------------------- Stage 2
	NonlinearTerm(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK2_u, RK_data->RK2_b, N);
	for (int i = 0; i < N; ++i) {
		// Get tmp index
		n = i + 2;

		// Get the integrating factor
		int_fac_u = cexp(-sys_vars->NU * dt * run_data->k[i] * run_data->k[i] / 2.0);
		int_fac_b = cexp(-sys_vars->ETA * dt * run_data->k[i] * run_data->k[i] / 2.0);

		// Update temorary velocity term
		RK_data->RK_u_tmp[n] = run_data->u[n] * int_fac_u + dt * RK4_A32 * RK_data->RK2_u[n];
		#if defined(__MAGNETO)
		// Update temporary magnetic term
		RK_data->RK_b_tmp[n] = run_data->b[n] * int_fac_b + dt * RK4_A32 * RK_data->RK2_b[n];		
		#endif
	}

	// ----------------------- Stage 3
	NonlinearTerm(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK3_u, RK_data->RK3_b, N);
	for (int i = 0; i < N; ++i) {
		// Get tmp index
		n = i + 2;

		// Get the integrating factor
		int_fac_u = cexp(-sys_vars->NU * dt * run_data->k[i] * run_data->k[i]);
		int_fac_b = cexp(-sys_vars->ETA * dt * run_data->k[i] * run_data->k[i]);

		// Update temorary velocity term
		RK_data->RK_u_tmp[n] = run_data->u[n] * int_fac_u + dt * RK4_A43 * RK_data->RK3_u[n] * sqrt(int_fac_u);
		#if defined(__MAGNETO)
		// Update temporary magnetic term
		RK_data->RK_b_tmp[n] = run_data->b[n] * int_fac_b + dt * RK4_A43 * RK_data->RK3_b[n] * sqrt(int_fac_b);		
		#endif
	}

	// ----------------------- Stage 2
	NonlinearTerm(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK4_u, RK_data->RK4_b, N);


	/////////////////////
	/// UPDATE STEP
	/////////////////////
	for (int i = 0; i < N; ++i) {
		// Get tmp index
		n = i + 2;

		// Get the integrating factors
		int_fac_u = cexp(-sys_vars->NU * dt * run_data->k[i] * run_data->k[i]);
		int_fac_b = cexp(-sys_vars->ETA * dt * run_data->k[i] * run_data->k[i]);

		#if defined(PHASE_ONLY)
		// Pre-record the amplitudes for resetting after update step
		po_norm_fac_u = cabs(run_data->u[n]);
		#if defined(__MAGNETO)
		po_norm_fac_b = cabs(run_data->b[n]);
		#endif
		#endif

		///-------------------- Update Step
		// Update the new velocity field
		run_data->u[n] = int_fac_u * run_data->u[n] + dt * RK4_B1 * int_fac_u * RK_data->RK1_u[n] + dt * RK4_B2 * sqrt(int_fac_u) * RK_data->RK2_u[n] + dt * RK4_B3 * sqrt(int_fac_u) * RK_data->RK3_u[n] + dt * RK4_B4 * RK_data->RK4_u[n];
		// run_data->u[n] = int_fac_u * run_data->u[n] + (dt / 6.0) * (int_fac_u * RK_data->RK1_u[n] + 2.0 * sqrt(int_fac_u) * RK_data->RK2_u[n] + 2.0 * sqrt(int_fac_u) * RK_data->RK3_u[n] + RK_data->RK4_u[n]);
		#if defined(__MAGNETO)
		// Update the new magnetic field
		run_data->b[n] = int_fac_b * run_data->b[n] + dt * RK4_B1 * int_fac_b * RK_data->RK1_b[n] + dt * RK4_B2 * sqrt(int_fac_b) * RK_data->RK2_b[n] + dt * RK4_B3 * sqrt(int_fac_b) * RK_data->RK3_b[n] + dt * RK4_B4 * RK_data->RK4_b[n];
		// run_data->b[n] = int_fac_b * run_data->b[n] + (dt / 6.0) * (int_fac_b * RK_data->RK1_b[n] + 2.0 * sqrt(int_fac_b) * RK_data->RK2_b[n] + 2.0 * sqrt(int_fac_b) * RK_data->RK3_b[n] + RK_data->RK4_b[n]);
		#endif

		#if defined(PHASE_ONLY)
		// Reset the amplitudes 
		run_data->u[n] *= (po_norm_fac_u / cabs(run_data->u[n]));
		#if defined(__MAGNETO)
		run_data->b[n] *= (po_norm_fac_b / cabs(run_data->b[n]));
		#endif

		// Record the phases and amplitudes
		run_data->a_n[n]   = cabs(run_data->u[n]);
		run_data->phi_n[n] = carg(run_data->u[n]);
		run_data->b_n[n]   = cabs(run_data->b[n]);
		run_data->psi_n[n] = carg(run_data->b[n]);
		#endif
	}
}
#endif
/**
 * Function that performs the evluation of the nonlinear term
 * @param u        array containing the input velocity modes
 * @param b        array containing the input magnetic modes
 * @param u_nonlin array to hold the result of computing the nonlinear term for velocity field
 * @param b_nonlin array to hold the result of computing the nonlinear term for the magnetic field
 * @param N        int defining the number of shells
 */
void NonlinearTerm(fftw_complex* u, fftw_complex* b, fftw_complex* u_nonlin, fftw_complex* b_nonlin, const long int N) {

	// Initialize variables
	int n;
	fftw_complex u_tmp_1, u_tmp_2, u_tmp_3;
	#if defined(__MAGNETO)
	fftw_complex b_tmp_1, b_tmp_2, b_tmp_3;
	#endif 

	// -----------------------------------
	// Compute The Nonlinear Terms
	// -----------------------------------
	for (int i = 0; i < N; ++i) {
		// Get tmp array index
		n = i + 2;

		// -----------------------------------
		// Compute Temporary Terms
		// -----------------------------------
		// Compute the velocity temporary terms in the nonlinear term
		u_tmp_1 = u[n + 1] * u[n + 2];
		u_tmp_2 = u[n - 1] * u[n + 1];
		u_tmp_3 = u[n - 2] * u[n - 1];
		#if defined(__MAGNETO)
		// Update velocity temporary terms with magnetic field
		u_tmp_1 -= b[n + 1] * b[n + 2];
		u_tmp_2 -= b[n - 1] * b[n + 1];
		u_tmp_3 -= b[n - 2] * b[n - 1];

		// Compute the magnetic temporary terms in the nonlinear term
		b_tmp_1 = u[n + 1] * b[n + 2] - b[n + 1] * u[n + 2];
		b_tmp_2 = u[n - 1] * b[n + 1] - b[n - 1] * u[n + 1];
		b_tmp_3 = u[n - 2] * b[n - 1] - b[n - 2] * u[n - 1];
		#endif
		
		// -----------------------------------
		// Compute Nonlinear Terms
		// -----------------------------------
		// Compute the nonlinear term for the velocity field
		u_nonlin[n] = I * run_data->k[i] * conj(u_tmp_1 - 0.25 * u_tmp_2 - 0.125 * u_tmp_3); 
		#if defined(__MAGNETO)
		// Compute the nonlinear term for the magnetic field
		b_nonlin[n] = I * run_data->k[i] * (1.0 / 6.0) * conj(b_tmp_1 + b_tmp_2 + b_tmp_3); 
		#endif
	}
}
/**
 * Function to compute the initial condition for the integration
 */
void InitialConditions(const long int N) {

	// Initialize variables
	
	// ------------------------------------------------
    // Set Seed for RNG
    // ------------------------------------------------
    srand(123456789);


	if(!(strcmp(sys_vars->u0, "N_SCALING"))) {
		// ------------------------------------------------
		// Scaling in N Initial Condition
		// ------------------------------------------------
		for (int i = 0; i < N; ++i) {
			// Initialize the velocity field
			run_data->u[i + 2] = 1.0 / pow(run_data->k[i], sys_vars->ALPHA) * cexp(I * pow(i + 1, 2.0));

			#if defined(__MAGNETO)
			// Initialize the magnetic field
			run_data->b[i + 2] = 1.0 / pow(run_data->k[i], sys_vars->BETA) * cexp(I * pow(i + 1, 4.0)) * 1e-2;
			#endif
		}	
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
		k[i] = K_0 * pow(LAMBDA, i);
	}
}
/**
 * Function to initialize all the integration time variables
 * @param t0           The initial time of the simulation
 * @param t            The current time of the simulaiton
 * @param dt           The timestep
 * @param T            The final time of the simulation
 * @param trans_steps  The number of iterations to perform before saving to file begins
 */
void InitializeIntegrationVariables(double* t0, double* t, double* dt, double* T, long int* trans_steps) {

	// -------------------------------
	// Get Time variables
	// -------------------------------
	// Compute integration time variables
	(*t0) = sys_vars->t0;
	(*t ) = sys_vars->t0;
	(*dt) = sys_vars->dt;
	(*T ) = sys_vars->T;
	sys_vars->min_dt = 10;
	sys_vars->max_dt = MIN_STEP_SIZE;

	// -------------------------------
	// Integration Counters
	// -------------------------------
	// Number of time steps and saving steps
	sys_vars->num_t_steps = ((*T) - (*t0)) / (*dt);
	if (sys_vars->TRANS_ITERS_FLAG == TRANSIENT_ITERS) {
		// Get the transient iterations
		(* trans_steps)       = (long int)(sys_vars->TRANS_ITERS_FRAC * sys_vars->num_t_steps);
		sys_vars->trans_iters = (* trans_steps);

		// Get the number of steps to perform before printing to file -> allowing for a transient fraction of these to be ignored
		sys_vars->num_print_steps = (sys_vars->num_t_steps >= sys_vars->SAVE_EVERY ) ? (sys_vars->num_t_steps - sys_vars->trans_iters) / sys_vars->SAVE_EVERY : sys_vars->num_t_steps - sys_vars->trans_iters;	 
		printf("Total Iters: %ld\t Saving Iters: %ld\t Transient Steps: %ld\n", sys_vars->num_t_steps, sys_vars->num_print_steps, sys_vars->trans_iters);
	}
	else {
		// Get the transient iterations 
		(* trans_steps)       = 0;
		sys_vars->trans_iters = (* trans_steps);

		// Get the number of steps to perform before printing to file
		sys_vars->num_print_steps = (sys_vars->num_t_steps >= sys_vars->SAVE_EVERY ) ? sys_vars->num_t_steps / sys_vars->SAVE_EVERY + 1 : sys_vars->num_t_steps + 1; // plus one to include initial condition
		printf("Total Iters: %ld\t Saving Iters: %ld\n", sys_vars->num_t_steps, sys_vars->num_print_steps);
	}

	// Variable to control how ofter to print to screen -> set it to half the saving to file steps
	sys_vars->print_every = (sys_vars->num_t_steps >= 10 ) ? (int)sys_vars->SAVE_EVERY : 1;
}
/**
 * Used to print update to screen
 * @param iters          The current iteration of the integration
 * @param t              The current time in the simulation
 * @param dt             The current timestep in the simulation
 * @param T              The final time of the simulation
 * @param save_data_indx The saving index for output data
 */
void PrintUpdateToTerminal(int iters, double t, double dt, double T, int save_data_indx) {

	// Initialize variables
	
	#if defined(__MAGNETO)
	printf("Iter: %d/%ld\tt: %1.6lf/%1.3lf\tdt: %1.6g\tKE: %6.6g\tHEL: %6.6g\tX-HEL: %6.6g\n", iters, sys_vars->num_t_steps, t, T, dt, run_data->tot_energy[save_data_indx], run_data->tot_hel[save_data_indx], run_data->tot_cross_hel[save_data_indx]);
	#else
	printf("Iter: %d/%ld\tt: %1.6lf/%1.3lf\tdt: %1.6g\tKE: %6.6g\n", iters, sys_vars->num_t_steps, t, T, dt, run_data->tot_energy[save_data_indx]);
	#endif

}
/**
 * Function that checks the system to see if it is ok to continue integrations. Checks for blow up, timestep and iteration limits etc
 * @param dt    The updated timestep for the next iteration
 * @param iters The number of iterations for the next iteration
 */
void SystemCheck(double dt, int iters) {

	// -------------------------------
	// Check Stopping Criteria 
	// -------------------------------
	if (dt <= MIN_STEP_SIZE) {
		fprintf(stderr, "\n["YELLOW"SOVLER FAILURE"RESET"]--- Timestep has become too small to continue at Iter: ["CYAN"%d"RESET"]\n-->> Exiting!!!\n", iters);
		exit(1);		
	}
	else if (iters >= MAX_ITERS) {
		fprintf(stderr, "\n["YELLOW"SOVLER FAILURE"RESET"]--- The maximum number of iterations has been reached at Iter: ["CYAN"%d"RESET"]\n-->> Exiting!!!\n", iters);
		exit(1);		
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

	// -------------------------------
	// Allocate Integration Variables 
	// -------------------------------
	// Runge-Kutta Integration arrays
	RK_data->RK1_u       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	if (RK_data->RK1_u == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK1 Velocity");
		exit(1);
	}
	RK_data->RK2_u       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	if (RK_data->RK2_u == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK2 Velocity");
		exit(1);
	}
	RK_data->RK3_u       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	if (RK_data->RK3_u == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK3 Velocity");
		exit(1);
	}
	RK_data->RK4_u       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	if (RK_data->RK4_u == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK4 Velocity");
		exit(1);
	}
	#if defined(__MAGNETO)
	RK_data->RK1_b       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	if (RK_data->RK1_b == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK1 Magnetic");
		exit(1);
	}
	RK_data->RK2_b       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	if (RK_data->RK2_b == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK2 Magnetic");
		exit(1);
	}
	RK_data->RK3_b       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	if (RK_data->RK3_b == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK3 Magnetic");
		exit(1);
	}
	RK_data->RK4_b       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	if (RK_data->RK4_b == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK4 Magnetic");
		exit(1);
	}
	#endif
	RK_data->RK_u_tmp    = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	if (RK_data->RK_u_tmp == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK_u_tmp");
		exit(1);
	}
	RK_data->RK_b_tmp    = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	if (RK_data->RK_b_tmp == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK_b_tmp");
		exit(1);
	}

	// -------------------------------
	// Initialize All Data 
	// -------------------------------
	for (int i = 0; i < N + 4; ++i) {
		// Initialize the arrays
		if (i < N) {
			run_data->k[i] = 0;
		}
		run_data->u[i]        = 0.0 + 0.0 * I;
		RK_data->RK1_u[i]    = 0.0 + 0.0 * I;
		RK_data->RK2_u[i]    = 0.0 + 0.0 * I;
		RK_data->RK3_u[i]    = 0.0 + 0.0 * I;
		RK_data->RK4_u[i]    = 0.0 + 0.0 * I;
		RK_data->RK_u_tmp[i] = 0.0 + 0.0 * I;
		#if defined(PHASE_ONLY)
		run_data->a_n[i]   = 0.0;
		run_data->phi_n[i] = 0.0;
		#endif
		#if defined(__MAGNETO)
		run_data->b[i]     = 0.0 + 0.0 * I;
		RK_data->RK1_b[i] = 0.0 + 0.0 * I;
		RK_data->RK2_b[i] = 0.0 + 0.0 * I;
		RK_data->RK3_b[i] = 0.0 + 0.0 * I;
		RK_data->RK4_b[i] = 0.0 + 0.0 * I;
		RK_data->RK_b_tmp[i] = 0.0 + 0.0 * I;
		#if defined(PHASE_ONLY)
		run_data->b_n[i]   = 0.0;
		run_data->psi_n[i] = 0.0;
		#endif
		#endif
	}
}
/**
 * Wrapper function that frees any memory dynamcially allocated in the programme
 * @param RK_data Pointer to a struct contaiing the integraiont arrays
 */
void FreeMemory(RK_data_struct* RK_data) {

	// ------------------------
	// Free memory 
	// ------------------------
	// Free shell variables
	fftw_free(run_data->k);

	// Free system variables
	fftw_free(run_data->u);
	#if defined(PHASE_ONLY)
	fftw_free(run_data->a_n);
	fftw_free(run_data->phi_n);
	#endif
	#if defined(__MAGNETO)
	fftw_free(run_data->b);
	#if defined(PHASE_ONLY)
	fftw_free(run_data->b_n);
	fftw_free(run_data->psi_n);
	#endif
	#endif
	#if defined(__SYS_MEASURES)
	fftw_free(run_data->tot_energy);
	#if defined(__MAGNETO)
	fftw_free(run_data->tot_hel);
	fftw_free(run_data->tot_cross_hel);
	#endif
	#endif
	#if defined(__TIME)
	fftw_free(run_data->time);
	#endif
	#if defined(__ENRG_FLUX)
	fftw_free(run_data->energy_flux);
	fftw_free(run_data->energy_diss);
	#endif

	// Free integration variables
	fftw_free(RK_data->RK1_u);
	fftw_free(RK_data->RK2_u);
	fftw_free(RK_data->RK3_u);
	fftw_free(RK_data->RK4_u);
	fftw_free(RK_data->RK_u_tmp);
	#if defined(__MAGNETO)
	fftw_free(RK_data->RK1_b);
	fftw_free(RK_data->RK2_b);
	fftw_free(RK_data->RK3_b);
	fftw_free(RK_data->RK4_b);
	fftw_free(RK_data->RK_b_tmp);
	#endif
}
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------