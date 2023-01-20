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
#include "stats.h"
#include "hdf5_funcs.h"
#include "utils.h"
#include "solver.h"
#include "sys_msr.h"
#include "phase_sync.h"
// ---------------------------------------------------------------------
//  Global Variables
// ---------------------------------------------------------------------
// Define RK4 variables - Butcher Tableau
#if defined(RK4) || defined(INT_FAC_RK4) || defined(AB4CN)
static const double RK4_C2 = 0.5, 	  RK4_A21 = 0.5, \
				  	RK4_C3 = 0.5,	           					RK4_A32 = 0.5, \
				  	RK4_C4 = 1.0,                      									   RK4_A43 = 1.0, \
				              	 	  RK4_B1 = 1.0/6.0, 		RK4_B2  = 1.0/3.0, 		   RK4_B3  = 1.0/3.0, 		RK4_B4 = 1.0/6.0;
#endif
#if defined(AB4CN)
static const double AB4_1 = 55.0/24.0, AB4_2 = -59.0/24.0,		AB4_3 = 37.0/24.0,			AB4_4 = -3.0/8.0;
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

	// -------------------------------
	// Initialize the System
	// -------------------------------
	// Initialize the collocation points and wavenumber space 
	InitializeShellWavenumbers(run_data->k, N);

	// Initialize stats objects if required
	#if defined(STATS)
	InitializeStats();
	#endif
	
	// Get initial conditions - seed for random number generator is set here
	// If input file is selected the function to read input file is called from here
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
	// Forcing
	// -------------------------------
	// Initialize the forcing
	InitializeForicing(N, dt);
	
	// -------------------------------
	// Create & Open Output File
	// -------------------------------
	// Inialize system measurables
	InitializeSystemMeasurables(RK_data);

	// Initialize the phase sync objects
	#if defined(__CONSERVED_PHASES) || defined(__PHASE_SYNC) || defined(__PHASE_SYNC_STATS)
	InitializePhaseSyncObjects();
	#endif

	// Create and open the output file - also write initial conditions to file
	CreateOutputFilesWriteICs(N);

	// Print update of the initial conditions to the terminal
	PrintUpdateToTerminal(0, t0, dt, T, 0);


	//////////////////////////////
	// Begin Integration
	//////////////////////////////
	t         += dt;
	long int iters = 1;
	long int save_data_indx;
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
		#if defined(INT_FAC_RK4)
		IntFacRK4Step(dt, N, RK_data);
		#endif
		#if defined(RK4)
		RK4Step(dt, N, RK_data);
		#endif
		#if defined(AB4CN)
		AB4CNStep(dt, (long int)iters, N, RK_data);
		#endif

		// -------------------------------
		// Compute System Quantities
		// -------------------------------
		// Compute stats
		#if defined(STATS)
		if (iters < trans_steps) {
			ComputeStats(iters, save_data_indx);
		}
		#endif

		// -------------------------------
		// Write To File
		// -------------------------------
		if (iters % sys_vars->SAVE_EVERY == 0) {

			// Record System Measurables
			ComputeSystemMeasurables(t, save_data_indx, RK_data);

			// printf("Eddy Turnover: %lf\n", pow(run_data->k[2] / K_0, 2.0 / 3.0));

			// Compute the conserved phases and phase order parameters
			#if defined(__CONSERVED_PHASES) || defined(__PHASE_SYNC) || defined(__PHASE_SYNC_STATS)
			ComputePhaseSyncData(iters);
			#endif

			// Compute stats
			#if defined(STATS)
			ComputeStats(iters, save_data_indx);
			#endif

			// If and when transient steps are complete write to file
			if (iters >= trans_steps) {
				// Write the appropriate datasets to file 
				WriteDataToFile(t, iters, save_data_indx);
				
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
		SystemCheck(dt, iters, save_data_indx);
	}
	//////////////////////////////
	// End Integration
	//////////////////////////////
	
	// ------------------------------- 
	// Final Writes to Output File
	// -------------------------------
	// Write the final state to file if it hasn't been printed already
	if (save_data_indx < sys_vars->num_print_steps) {
		printf("\n\n...Writing Final State To File!\n\n");
		WriteDataToFile(t, iters, sys_vars->num_print_steps - 1);
	}
	// Compute System Measures on final state
	ComputeSystemMeasurables(t, sys_vars->num_print_steps - 1, RK_data);

	// Write the stats and system measures to file
	FinalWriteAndCloseOutputFile(N, iters, sys_vars->num_print_steps - 1);

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
#if defined(INT_FAC_RK4)
void IntFacRK4Step(const double dt, const long int N, RK_data_struct* RK_data) {

	// Initialize vairables
	int n;
	double complex int_fac_u;
	double complex int_fac_u_1;
	#if defined(__MAGNETO)
	double complex int_fac_b;
	double complex int_fac_b_1;
	#endif
	#if defined(PHASE_ONLY_FXD_AMP)
	double po_norm_fac_u;
	#if defined(__MAGNETO)
	double po_norm_fac_b;
	#endif
	#endif

	///----------------------- Forcing
	// Compute the forcing for the current iteration
	ComputeForcing(dt, N);

	///----------------------- Input to RHS/Nonlinear Term
	// Get the input velocity and magnetic fields for the nonlinear term
	for (int i = 0; i < N; ++i) {
		// Get proper index
		n = i + 2;

		// Get the input fields
		#if defined(PHASE_ONLY)
		// Get the Fourier velocity from the Fourier phases and amplitudes
		RK_data->RK_u_tmp[n] = run_data->a_n[n] * cexp(I * run_data->phi_n[n]);
		#if defined(__MAGNETO)
		RK_data->RK_b_tmp[n] = run_data->b_n[n] * cexp(I * run_data->psi_n[n]);
		#endif
		#else
		// Get the Fourier velocity field
		RK_data->RK_u_tmp[n] = run_data->u[n];
		#if defined(__MAGNETO)
		RK_data->RK_b_tmp[n] = run_data->b[n];
		#endif
		#endif
	}

	/////////////////////
	/// RK STAGES
	/////////////////////
	// ----------------------- Stage 1
	NonlinearTermWithForcing(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK1_u, RK_data->RK1_b);
	for (int i = 0; i < N; ++i) {
		// Get proper indx
		n = i + 2;

		// Get the integrating factor
		int_fac_u = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] / 2.0);

		// Update temorary velocity term
		#if defined(PHASE_ONLY)
		RK_data->RK_u_tmp[n] = run_data->a_n[n] * cexp(I * (run_data->phi_n[n] * int_fac_u + dt * RK4_A21 * RK_data->RK1_u[n] * int_fac_u));
		#else
		RK_data->RK_u_tmp[n] = run_data->u[n] * int_fac_u + dt * RK4_A21 * RK_data->RK1_u[n] * int_fac_u;
		#endif
		#if defined(__MAGNETO)
		// Get the integrating factor
		int_fac_b = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] / 2.0);
		
		// Update temporary magnetic term
		#if defined(PHASE_ONLY)
		RK_data->RK_b_tmp[n] = run_data->b_n[n] * cexp(I * (run_data->psi_n[n] * int_fac_b + dt * RK4_A21 * RK_data->RK1_b[n] * int_fac_b));
		#else
		RK_data->RK_b_tmp[n] = run_data->b[n] * int_fac_b + dt * RK4_A21 * RK_data->RK1_b[n] * int_fac_b;
		#endif
		#endif
	}

	// ----------------------- Stage 2
	NonlinearTermWithForcing(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK2_u, RK_data->RK2_b);
	for (int i = 0; i < N; ++i) {
		// Get proper indx
		n = i + 2;

		// Get the integrating factor
		int_fac_u = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] / 2.0);

		// Update temorary velocity term
		#if defined(PHASE_ONLY)
		RK_data->RK_u_tmp[n] = run_data->a_n[n] * cexp(I * (run_data->phi_n[n] * int_fac_u + dt * RK4_A32 * RK_data->RK2_u[n]));
		#else
		RK_data->RK_u_tmp[n] = run_data->u[n] * int_fac_u + dt * RK4_A32 * RK_data->RK2_u[n];
		#endif
		#if defined(__MAGNETO)
		// Get the integrating factor
		int_fac_b = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] / 2.0);

		// Update temporary magnetic term
		#if defined(PHASE_ONLY)
		RK_data->RK_b_tmp[n] = run_data->b_n[n] * cexp(I * (run_data->b[n] * int_fac_b + dt * RK4_A32 * RK_data->RK2_b[n]));
		#else
		RK_data->RK_b_tmp[n] = run_data->b[n] * int_fac_b + dt * RK4_A32 * RK_data->RK2_b[n];		
		#endif
		#endif
	}

	// ----------------------- Stage 3
	NonlinearTermWithForcing(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK3_u, RK_data->RK3_b);
	for (int i = 0; i < N; ++i) {
		// Get proper indx
		n = i + 2;

		// Get the integrating factor
		int_fac_u   = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n]);
		int_fac_u_1 = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] / 2.0);

		// Update temorary velocity term
		#if defined(PHASE_ONLY)
		RK_data->RK_u_tmp[n] = run_data->a_n[n] * cexp(I * (run_data->phi_n[n] * int_fac_u + dt * RK4_A43 * RK_data->RK3_u[n] * int_fac_u_1));
		#else
		RK_data->RK_u_tmp[n] = run_data->u[n] * int_fac_u + dt * RK4_A43 * RK_data->RK3_u[n] * int_fac_u_1;
		#endif
		#if defined(__MAGNETO)
		// Get the integrating factor
		int_fac_b   = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n]);
		int_fac_b_1 = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] / 2.0);
		
		// Update temporary magnetic term
		#if defined(PHASE_ONLY)
		RK_data->RK_b_tmp[n] = run_data->b_n[n] * cexp(I * (run_data->b[n] * int_fac_b + dt * RK4_A43 * RK_data->RK3_b[n] * int_fac_b_1));
		#else
		RK_data->RK_b_tmp[n] = run_data->b[n] * int_fac_b + dt * RK4_A43 * RK_data->RK3_b[n] * int_fac_b_1;		
		#endif
		#endif
	}

	// ----------------------- Stage 4
	NonlinearTermWithForcing(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK4_u, RK_data->RK4_b);


	/////////////////////
	/// UPDATE STEP
	/////////////////////
	for (int i = 0; i < N; ++i) {
		// Get temporary index
		n = i + 2;

		// Pre-record the amplitudes for resetting after update step
		#if defined(PHASE_ONLY_FXD_AMP)
		po_norm_fac_u = cabs(run_data->u[n]);
		#if defined(__MAGNETO)
		po_norm_fac_b = cabs(run_data->b[n]);
		#endif
		#endif

		///-------------------- Update Step
		// Get the integrating factors
		int_fac_u   = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n]);
		int_fac_u_1 = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] / 2.0);

		// Update the new velocity field
		#if defined(PHASE_ONLY)
		run_data->phi_n[n] = int_fac_u * run_data->phi_n[n] + dt * RK4_B1 * int_fac_u * RK_data->RK1_u[n] + dt * RK4_B2 * int_fac_u_1 * RK_data->RK2_u[n] + dt * RK4_B3 * int_fac_u_1 * RK_data->RK3_u[n] + dt * RK4_B4 * RK_data->RK4_u[n];
		#else
		run_data->u[n] = int_fac_u * run_data->u[n] + dt * RK4_B1 * int_fac_u * RK_data->RK1_u[n] + dt * RK4_B2 * int_fac_u_1 * RK_data->RK2_u[n] + dt * RK4_B3 * int_fac_u_1 * RK_data->RK3_u[n] + dt * RK4_B4 * RK_data->RK4_u[n];
		#endif
		#if defined(__MAGNETO)
		// Get the integrating factors
		int_fac_b   = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n]);
		int_fac_b_1 = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] / 2.0);

		// Update the new magnetic field
		#if defined(PHASE_ONLY)
		run_data->psi_n[n] = int_fac_b * run_data->psi_n[n] + dt * RK4_B1 * int_fac_b * RK_data->RK1_b[n] + dt * RK4_B2 * int_fac_b_1 * RK_data->RK2_b[n] + dt * RK4_B3 * int_fac_b_1 * RK_data->RK3_b[n] + dt * RK4_B4 * RK_data->RK4_b[n];
		#else
		run_data->b[n] = int_fac_b * run_data->b[n] + dt * RK4_B1 * int_fac_b * RK_data->RK1_b[n] + dt * RK4_B2 * int_fac_b_1 * RK_data->RK2_b[n] + dt * RK4_B3 * int_fac_b_1 * RK_data->RK3_b[n] + dt * RK4_B4 * RK_data->RK4_b[n];
		#endif
		#endif


		///-------------------- Phase Only resetting
		#if defined(PHASE_ONLY_FXD_AMP)
		// Reset the amplitudes 
		run_data->u[n] *= (po_norm_fac_u / cabs(run_data->u[n]));
		
		// Record the phases and amplitudes
		run_data->a_n[n]   = cabs(run_data->u[n]);
		run_data->phi_n[n] = carg(run_data->u[n]);

		#if defined(__MAGNETO)
		run_data->b[n] *= (po_norm_fac_b / cabs(run_data->b[n]));
		
		// Record the phases and amplitudes
		run_data->b_n[n]   = cabs(run_data->b[n]);
		run_data->psi_n[n] = carg(run_data->b[n]);
		#endif
		#endif

		// #if defined(__MAGNETO)
		// printf("u[%d]:\t%1.16lf\t%1.16lf i\tb[%d]:\t%1.16lf\t%1.16lf i\n", i - 1, creal(run_data->u[i]), cimag(run_data->u[i]),  i - 1, creal(run_data->b[i]), cimag(run_data->b[i]));
		// #else
		// printf("u[%d]:\t%1.16lf\t%1.16lf i\n", i - 1, creal(run_data->u[i]), cimag(run_data->u[i]));		
		// #endif
	}
}
#endif
/**
 * Function to perform one step using the 4th order Runge-Kutta method
 * @param dt       The current timestep of the system
 * @param N        int defining the number of shells
 * @param RK_data  Struct pointing the Integration variables: stages, tmp arrays, rhs and arrays needed for NonlinearRHS function
 */
#if defined(RK4) || defined(AB4CN)
void RK4Step(const double dt, const long int N, RK_data_struct* RK_data) {

	// Initialize vairables
	int n;
	#if defined(PHASE_ONLY_FXD_AMP)
	double po_norm_fac_u;
	#if defined(__MAGNETO)
	double po_norm_fac_b;
	#endif
	#endif

	///----------------------- Forcing
	// Compute the forcing for the current iteration
	ComputeForcing(dt, N);

	///----------------------- Input to RHS/Nonlinear Term
	// Get the input velocity and magnetic fields for the nonlinear term
	for (int i = 0; i < N; ++i) {
		// Get proper index
		n = i + 2;

		// Get the input fields
		#if defined(PHASE_ONLY)
		// Get the Fourier velocity from the Fourier phases and amplitudes
		RK_data->RK_u_tmp[n] = run_data->a_n[n] * cexp(I * run_data->phi_n[n]);
		#if defined(__MAGNETO)
		RK_data->RK_b_tmp[n] = run_data->b_n[n] * cexp(I * run_data->psi_n[n]);
		#endif
		#else
		// Get the Fourier velocity field
		RK_data->RK_u_tmp[n] = run_data->u[n];
		#if defined(__MAGNETO)
		RK_data->RK_b_tmp[n] = run_data->b[n];
		#endif
		#endif
	}

	/////////////////////
	/// RK STAGES
	/////////////////////
	// ----------------------- Stage 1	
	NonlinearTermWithForcing(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK1_u, RK_data->RK1_b);
	for (int i = 0; i < N; ++i) {
		// Get proper index
		n = i + 2;

		// Update temporary input for nonlinear term
		#if defined(PHASE_ONLY)
		RK_data->RK_u_tmp[n] = run_data->a_n[n] * cexp(I * (run_data->phi_n[n] + dt * RK4_A21 * RK_data->RK1_u[n]));
		#if defined(__MAGNETO)
		RK_data->RK_b_tmp[n] = run_data->a_n[n] * cexp(I * (run_data->psi_n[n] + dt * RK4_A21 * RK_data->RK1_b[n]));
		#endif
		#else
		// Add dissipative and forcing terms & Update temorary velocity term
		RK_data->RK1_u[n] -= sys_vars->NU * run_data->k[n] * run_data->k[n] * run_data->u[n];
		RK_data->RK_u_tmp[n] = run_data->u[n] + dt * RK4_A21 * RK_data->RK1_u[n];
		#if defined(__MAGNETO)
		// Add dissipative and forcing terms & Update temporary magnetic term
		RK_data->RK1_b[n] -= sys_vars->ETA * run_data->k[n] * run_data->k[n] * run_data->b[n];
		RK_data->RK_b_tmp[n] = run_data->b[n] + dt * RK4_A21 * RK_data->RK1_b[n];
		#endif
		#endif
	}

	// ----------------------- Stage 2
	NonlinearTermWithForcing(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK2_u, RK_data->RK2_b);
	for (int i = 0; i < N; ++i) {
		// Get proper index
		n = i + 2;

		// Update temporary input for nonlinear term
		#if defined(PHASE_ONLY)
		RK_data->RK_u_tmp[n] = run_data->a_n[n] * cexp(I * (run_data->phi_n[n] + dt * RK4_A32 * RK_data->RK2_u[n]));
		#if defined(__MAGNETO)
		RK_data->RK_b_tmp[n] = run_data->b_n[n] * cexp(I * (run_data->psi_n[n] + dt * RK4_A32 * RK_data->RK2_b[n]));
		#endif
		#else
		// Add dissipative & Update temorary velocity term
		RK_data->RK2_u[n] -= sys_vars->NU * run_data->k[n] * run_data->k[n] * run_data->u[n];
		RK_data->RK_u_tmp[n] = run_data->u[n] + dt * RK4_A32 * RK_data->RK2_u[n];
		#if defined(__MAGNETO)
		// Add dissipative & Update temporary magnetic term
		RK_data->RK2_b[n] -= sys_vars->ETA * run_data->k[n] * run_data->k[n] * run_data->b[n];
		RK_data->RK_b_tmp[n] = run_data->b[n] + dt * RK4_A32 * RK_data->RK2_b[n];
		#endif
		#endif
	}

	// ----------------------- Stage 3
	NonlinearTermWithForcing(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK3_u, RK_data->RK3_b);
	for (int i = 0; i < N; ++i) {
		// Get proper index
		n = i + 2;

		#if defined(PHASE_ONLY)
		RK_data->RK_u_tmp[n] = run_data->a_n[n] * cexp(I * (run_data->phi_n[n] + dt * RK4_A43 * RK_data->RK3_u[n]));
		#if defined(__MAGNETO)
		RK_data->RK_b_tmp[n] = run_data->b_n[n] * cexp(I * (run_data->psi_n[n] + dt * RK4_A43 * RK_data->RK3_b[n]));
		#endif
		#else
		// Add dissipative & Update temorary velocity term
		RK_data->RK3_u[n] -= sys_vars->NU * run_data->k[n] * run_data->k[n] * run_data->u[n];
		RK_data->RK_u_tmp[n] = run_data->u[n] + dt * RK4_A43 * RK_data->RK3_u[n];
		#if defined(__MAGNETO)
		// Add dissipative & Update temporary magnetic term
		RK_data->RK3_b[n] -= sys_vars->ETA * run_data->k[n] * run_data->k[n] * run_data->b[n];
		RK_data->RK_b_tmp[n] = run_data->b[n] + dt * RK4_A43 * RK_data->RK3_b[n];
		#endif
		#endif
	}

	// ----------------------- Stage 4
	NonlinearTermWithForcing(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK4_u, RK_data->RK4_b);
	// Add the linear term
	for (int i = 0; i < N; ++i) {
		// Get proper index
		n = i + 2;

		#if !defined(PHASE_ONLY)
		RK_data->RK4_u[n] -= sys_vars->NU * run_data->k[n] * run_data->k[n] * run_data->u[n];
		#if defined(__MAGNETO)
		RK_data->RK4_b[n] -= sys_vars->ETA * run_data->k[n] * run_data->k[n] * run_data->b[n];
		#endif
		#endif
	}

	/////////////////////
	/// UPDATE STEP
	/////////////////////
	for (int i = 0; i < N; ++i) {
		// Get tmp index
		n = i + 2;

		#if defined(PHASE_ONLY_FXD_AMP)
		// Pre-record the amplitudes for resetting after update step
		po_norm_fac_u = cabs(run_data->u[n]);
		#if defined(__MAGNETO)
		po_norm_fac_b = cabs(run_data->b[n]);
		#endif
		#endif

		///-------------------- Update Step
		#if defined(PHASE_ONLY)
		// Update the new velocity field
		run_data->phi_n[n] = run_data->phi_n[n] + dt * RK4_B1 * RK_data->RK1_u[n] + dt * RK4_B2 * RK_data->RK2_u[n] + dt * RK4_B3 * RK_data->RK3_u[n] + dt * RK4_B4 * RK_data->RK4_u[n];
		#if defined(__MAGNETO)
		// Update the new magnetic field
		run_data->psi_n[n] = run_data->psi_n[n] + dt * RK4_B1 * RK_data->RK1_b[n] + dt * RK4_B2 * RK_data->RK2_b[n] + dt * RK4_B3 * RK_data->RK3_b[n] + dt * RK4_B4 * RK_data->RK4_b[n];
		#endif
		#else
		// Update the new velocity field
		run_data->u[n] = run_data->u[n] + dt * RK4_B1 * RK_data->RK1_u[n] + dt * RK4_B2 * RK_data->RK2_u[n] + dt * RK4_B3 * RK_data->RK3_u[n] + dt * RK4_B4 * RK_data->RK4_u[n];
		#if defined(__MAGNETO)
		// Update the new magnetic field
		run_data->b[n] = run_data->b[n] + dt * RK4_B1 * RK_data->RK1_b[n] + dt * RK4_B2 * RK_data->RK2_b[n] + dt * RK4_B3 * RK_data->RK3_b[n] + dt * RK4_B4 * RK_data->RK4_b[n];
		#endif
		#endif

		#if defined(PHASE_ONLY_FXD_AMP)
		// Reset the amplitudes 
		run_data->u[n] *= (po_norm_fac_u / cabs(run_data->u[n]));

		// Record the phases and amplitudes
		run_data->a_n[n]   = cabs(run_data->u[n]);
		run_data->phi_n[n] = carg(run_data->u[n]);
		#if defined(__MAGNETO)
		run_data->b[n] *= (po_norm_fac_b / cabs(run_data->b[n]));

		// Record the phases and amplitudes
		run_data->b_n[n]   = cabs(run_data->b[n]);
		run_data->psi_n[n] = carg(run_data->b[n]);
		#endif
		#endif
	}
	// for (int i = 0; i < N + 4; ++i) {
	// 	#if defined(__MAGNETO)
	// 	printf("u[%d]:\t%1.16lf\t%1.16lf i\tb[%d]:\t%1.16lf\t%1.16lf i\n", i - 1, creal(run_data->u[i]), cimag(run_data->u[i]),  i - 1, creal(run_data->b[i]), cimag(run_data->b[i]));
	// 	#elif defined(PHASE_ONLY)
	// 	printf("a[%d]:\t%1.16lf\t--\tp[%d]:\t%1.16lf\t\n", i - 1, run_data->a_n[i], i - 1, RK_data->RK1_u[i]);
	// 	#else
	// 	printf("u[%d]:\t%1.16lf\t%1.16lf i\n", i - 1, creal(run_data->u[i]), cimag(run_data->u[i]));		
	// 	#endif
	// }	
}
#endif
#if defined(AB4CN)
/**
 * Function to perform an integration step using the 4th order Adams Bashforth shceme with Crank Nicholson differencing on the viscosity term.
 * @param dt      The current timestep in the simulation
 * @param iters   The current iteration of the simulation
 * @param N       The number of shell modes	
 * @param RK_data Struct containing the integration arrays
 */
void AB4CNStep(const double dt, const long iters, const long int N, RK_data_struct* RK_data) {

	// Initialize variables
	int n;
	double D_fac_u;
	#if defined(__MAGNETO)	
	double D_fac_b;
	#endif
	#if defined(PHASE_ONLY_FXD_AMP)
	double po_norm_fac_u;
	#if defined(__MAGNETO)
	double po_norm_fac_b;
	#endif
	#endif

	/////////////////////
	/// AB Pre Steps
	/////////////////////
	if (iters <= RK_data->AB_pre_steps) {
		// -----------------------------------
		// Perform RK4 Step
		// -----------------------------------
		// March the field forward in time using RK4 step
		RK4Step(dt, N, RK_data);

		// Save the nonlinear term for each pre step for use in the update step of the AB4CN scheme
		memcpy(&(RK_data->AB_tmp_nonlin_u[iters - 1][0]), RK_data->RK1_u, sizeof(double complex) * (N + 4));
		#if defined(__MAGNETO)
		memcpy(&(RK_data->AB_tmp_nonlin_b[iters - 1][0]), RK_data->RK1_b, sizeof(double complex) * (N + 4));		
		#endif
	}
	else {
		// -----------------------------------
		// Compute Forcing
		// -----------------------------------
		// Compute the forcing for the current iteration
		ComputeForcing(dt, N);

		// -----------------------------------
		// Compute RHS
		// -----------------------------------
		// Get the input velocity and magnetic fields for the nonlinear term
		for (int i = 0; i < N; ++i) {
			// Get proper index
			n = i + 2;

			// Get the input fields
			#if defined(PHASE_ONLY)
			// Get the Fourier velocity from the Fourier phases and amplitudes
			RK_data->RK_u_tmp[n] = run_data->a_n[n] * cexp(I * run_data->phi_n[n]);
			#if defined(__MAGNETO)
			RK_data->RK_b_tmp[n] = run_data->b_n[n] * cexp(I * run_data->psi_n[n]);
			#endif
			#else
			// Get the Fourier velocity field
			RK_data->RK_u_tmp[n] = run_data->u[n];
			#if defined(__MAGNETO)
			RK_data->RK_b_tmp[n] = run_data->b[n];
			#endif
			#endif
		}

		// Get the nonlinear term for the current step
		NonlinearTermWithForcing(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->AB_tmp_u, RK_data->AB_tmp_b);
		

		/////////////////////
		/// AB Update Step
		/////////////////////
		for (int i = 0; i < N; ++i) {
			// Get tmp index
			n = i + 2;

			// Pre-record the amplitudes for resetting after update step
			#if defined(PHASE_ONLY_FXD_AMP)
			po_norm_fac_u = cabs(run_data->u[n]);
			#if defined(__MAGNETO)
			po_norm_fac_b = cabs(run_data->b[n]);
			#endif
			#endif

			#if defined(PHASE_ONLY)
			run_data->phi_n[n] = run_data->phi_n[n] + dt * (AB4_1 * RK_data->AB_tmp_u[n] + AB4_2 * RK_data->AB_tmp_nonlin_u[2][n] + AB4_3 * RK_data->AB_tmp_nonlin_u[1][n] + AB4_4 * RK_data->AB_tmp_nonlin_u[0][n]);
			#if defined(__MAGNETO)
			run_data->psi_n[n] = run_data->psi_n[n] + dt * (AB4_1 * RK_data->AB_tmp_b[n] + AB4_2 * RK_data->AB_tmp_nonlin_b[2][n] + AB4_3 * RK_data->AB_tmp_nonlin_b[1][n] + AB4_4 * RK_data->AB_tmp_nonlin_b[0][n]);
			#endif
			#else
			// Compute the D factor for the velocity field
			D_fac_u = dt * (sys_vars->NU * run_data->k[i] * run_data->k[i]);

			// Update the new velocity field
			run_data->u[n] = run_data->u[n] * ((2.0 - D_fac_u) / (2.0 + D_fac_u)) + (2.0 * dt / (2.0 + D_fac_u)) * (AB4_1 * RK_data->AB_tmp_u[n] + AB4_2 * RK_data->AB_tmp_nonlin_u[2][n] + AB4_3 * RK_data->AB_tmp_nonlin_u[1][n] + AB4_4 * RK_data->AB_tmp_nonlin_u[0][n]);
			#if defined(__MAGNETO)
			// Compute the D factor for the velocity field
			D_fac_b = dt * (sys_vars->ETA * run_data->k[i] * run_data->k[i]);

			// Update the new magnetic field
			run_data->b[n] = run_data->b[n] * ((2.0 - D_fac_b) / (2.0 + D_fac_b)) + (2.0 * dt / (2.0 + D_fac_b)) * (AB4_1 * RK_data->AB_tmp_b[n] + AB4_2 * RK_data->AB_tmp_nonlin_b[2][n] + AB4_3 * RK_data->AB_tmp_nonlin_b[1][n] + AB4_4 * RK_data->AB_tmp_nonlin_b[0][n]);
			#endif
			#endif

			// Reset the amplitudes 
			#if defined(PHASE_ONLY_FXD_AMP)
			run_data->u[n] *= (po_norm_fac_u / cabs(run_data->u[n]));

			// Record the phases and amplitudes
			run_data->a_n[n]   = cabs(run_data->u[n]);
			run_data->phi_n[n] = carg(run_data->u[n]);
			#if defined(__MAGNETO)
			run_data->b[n] *= (po_norm_fac_b / cabs(run_data->b[n]));

			// Record the phases and amplitudes
			run_data->b_n[n]   = cabs(run_data->b[n]);
			run_data->psi_n[n] = carg(run_data->b[n]);
			#endif
			#endif
		}

		// -----------------------------------
		// Update Previous Nonlinear Terms
		// -----------------------------------
		// Update the previous Nonlinear term arrays for next iteration
		memcpy(&(RK_data->AB_tmp_nonlin_u[0][0]), &(RK_data->AB_tmp_nonlin_u[1][0]), sizeof(double complex) * (N + 4));
		memcpy(&(RK_data->AB_tmp_nonlin_u[1][0]), &(RK_data->AB_tmp_nonlin_u[2][0]), sizeof(double complex) * (N + 4));
		memcpy(&(RK_data->AB_tmp_nonlin_u[2][0]), RK_data->AB_tmp_u, sizeof(double complex) * (N + 4));
		#if defined(__MAGNETO)
		memcpy(&(RK_data->AB_tmp_nonlin_b[0][0]), &(RK_data->AB_tmp_nonlin_b[1][0]), sizeof(double complex) * (N + 4));
		memcpy(&(RK_data->AB_tmp_nonlin_b[1][0]), &(RK_data->AB_tmp_nonlin_b[2][0]), sizeof(double complex) * (N + 4));
		memcpy(&(RK_data->AB_tmp_nonlin_b[2][0]), RK_data->AB_tmp_b, sizeof(double complex) * (N + 4));
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
void NonlinearTerm(double complex* u, double complex* b, double complex* u_nonlin, double complex* b_nonlin, const long int N) {

	// Initialize variables
	int n;
	const double lambda_pow         = sys_vars->Lambda * sys_vars->Lambda;
	const double interact_coeff_u_1 = sys_vars->EPS / sys_vars->Lambda;
	const double interact_coeff_u_2 = (1.0 - sys_vars->EPS) / lambda_pow;
	#if defined(__MAGNETO)
	const double interact_coeff_b_1 = 1.0 - sys_vars->EPS - sys_vars->EPS_M;
	const double interact_coeff_b_2 = sys_vars->EPS_M / sys_vars->Lambda;
	const double interact_coeff_b_3 = (1.0 - sys_vars->EPS_M) / lambda_pow;
	#endif
	double complex u_tmp_1, u_tmp_2, u_tmp_3;
	#if defined(__MAGNETO)
	double complex b_tmp_1, b_tmp_2, b_tmp_3;
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
		u_nonlin[n] = I * run_data->k[n] * conj(u_tmp_1 - interact_coeff_u_1 * u_tmp_2 - interact_coeff_u_2 *  u_tmp_3);
		#if defined(__MAGNETO)
		// Compute the nonlinear term for the magnetic field
		b_nonlin[n] = I * run_data->k[n] * conj(interact_coeff_b_1 * b_tmp_1 + interact_coeff_b_2 * b_tmp_2 + interact_coeff_b_3 * b_tmp_3); 
		#endif
		// printf("u[%d]:\t%lf\t%lfi\t-\tf[%d]:\t%lf\t%lfi\t----\t%lf\t%lf\n", i + 1, creal(u_nonlin[n]), cimag(u_nonlin[n]), i + i, creal(run_data->forcing_u[n]), cimag(run_data->forcing_u[n]), creal(u_nonlin[n]) / creal(run_data->forcing_u[n]), cimag(u_nonlin[n]) / cimag(run_data->forcing_u[n]));
	}
}
/**
 * Wrapper function to get the approriate nonlinear term for the system being solver -> either phase only or full model
 * @param input_u  Input array for the velocity field, either velocity or the phases
 * @param input_b  Input array for the magnetic field, the magnetic field or the phases
 * @param output_u Output array for the velocity field
 * @param output_b Output array for the magnetic field
 */
void NonlinearTermWithForcing(double complex* input_u, double complex* input_b, double complex* output_u, double complex* output_b) {

	// Initialize variables
	int n;
	double complex tmp_output_u;
	#if defined(__MAGNETO)
	double complex tmp_output_b;
	#endif

	// ----------------------------------------
	// Compute the Nonlinear Term with Forcing
	// ----------------------------------------
	NonlinearTerm(input_u, input_b, output_u, output_b, sys_vars->N);

	// -----------------------------------
	// Add Forcing
	// -----------------------------------
	// Add forcing here
	AddForcing(output_u, output_b);

	// ----------------------------------------
    // Get the Output Array
    // ----------------------------------------
    for (int i = 0; i < sys_vars->N; ++i) {
		// Get tmp array index
		n = i + 2;

		// Store tmp output value
		tmp_output_u = output_u[n];
		#if defined(__MAGNETO)
		tmp_output_b = output_b[n];
		#endif

		// Get the appropriate output value depending on the system being solved
		#if defined(PHASE_ONLY)
		// Get the Phase only RHS
		if (run_data->a_n[n] != 0.0) {
			output_u[n] = cimag(tmp_output_u * cexp(-I * run_data->phi_n[n])) / run_data->a_n[n];
			#if defined(__MAGNETO)
			output_b[n] = cimag(tmp_output_b * cexp(-I * run_data->psi_n[n])) / run_data->b_n[n];
			#endif
		}
		else {
			output_u[n] = 0.0;
			#if defined(__MAGNETO)
			output_b[n] = 0.0;
			#endif
		}
		#else
		// Get the full systems Nonlinear term and forcing
		output_u[n] = tmp_output_u;
		#if defined(__MAGNETO)
		output_b[n] = tmp_output_b;
		#endif
		#endif
	}
}
// #endif
/**
 * Function to compute the initial condition for the integration
 */
void InitialConditions(const long int N) {

	// Initialize variables
	double r1, r3;
	#if defined(__MAGNETO)
	double r2, r4;
	#endif

	// ------------------------------------------------
    // Set Seed for RNG
    // ------------------------------------------------
    srand(123456789);


	// ------------------------------------------------
    // Check if Reading From Input File
    // ------------------------------------------------
    if (sys_vars->INPUT_FILE_FLAG == INPUT_FILE) {
		// ------------------------------------------------
	    // Read in Initial Condition From File
	    // ------------------------------------------------
	    ReadInputFile(N);
    }	
    else {
		for (int i = 0; i < N + 4; ++i) {

			// Initialize the edges shells
			if (i < 2 || i > N + 1) {
				// Initialize the velocity field
				run_data->u[i] = 0.0 + 0.0 * I;
				// Initialize the phases and amplitudes
				#if defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY) 
				run_data->a_n[i]   = 0.0;
				run_data->phi_n[i] = 0.0;
				#endif

				// Initialize the magnetic field
				#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
				run_data->b[i] = 0.0 + 0.0 * I;
				// Initialize the phases and amplitudes
				#if defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY) 
				run_data->b_n[i]   = 0.0;
				run_data->psi_n[i] = 0.0;
				#endif
				#endif
			}
			// Initialize the interior shells
			else {
				if(!(strcmp(sys_vars->u0, "N_SCALING"))) {
					// ------------------------------------------------
					// Scaling in N Initial Condition
					// ------------------------------------------------
					// Initialize the velocity field
					run_data->u[i] = 1.0 / pow(run_data->k[i], sys_vars->ALPHA) * cexp(I * pow(i - 1, 2.0)) / sqrt(75);
					// Record the phases and amplitudes
					#if defined(PHASE_ONLY) || defined(PHASE_ONLY_FXD_AMP)
					run_data->a_n[i]   = (1.0 / pow(run_data->k[i], sys_vars->ALPHA));
					run_data->phi_n[i] = pow(i - 1, 2.0);
					#endif

					// Initialize the magnetic field
					#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
					run_data->b[i] = 1.0 / pow(run_data->k[i], sys_vars->BETA) * cexp(I * pow(i - 1, 4.0)) * 1e-2 / sqrt(75);
					// Record the phases and amplitudes
					#if defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY)
					run_data->b_n[i]   = 1.0 / pow(run_data->k[i], sys_vars->BETA) * 1e-2;
					run_data->psi_n[i] = cexp(I * pow(i - 1, 4.0));
					#endif
					#endif
				}
				else if(!(strcmp(sys_vars->u0, "RANDOM"))) {
					// ------------------------------------------------
					// Default - Random Initial Conditions
					// ------------------------------------------------	
					// Get random uniform number
					r1 = (double)rand() / (double)RAND_MAX;
					#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
					r2 = (double)rand() / (double)RAND_MAX;
					#endif

					// Initialize the velocity field
					run_data->u[i] = 1.0 / pow(run_data->k[i], sys_vars->ALPHA) * cexp(I * r1 * 2.0 * M_PI); // *1e-3
					// Record the phases and amplitudes
					#if defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY)
					run_data->a_n[i]   = 1.0 / pow(run_data->k[i], sys_vars->ALPHA); // *1e-3
					run_data->phi_n[i] = r1 * 2.0 * M_PI;
					#endif

					// Initialize the magnetic field
					#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
					run_data->b[i] = 1.0 / pow(run_data->k[i], sys_vars->BETA) * cexp(I * r2 * 2.0 * M_PI); // *1e-3

					// Record the phases and amplitudes
					#if defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY)
					run_data->b_n[i]   = 1.0 / pow(run_data->k[i], sys_vars->BETA); // *1e-3
					run_data->psi_n[i] = r2 * 2.0 * M_PI;
					#endif
					#endif
				}
				else if(!(strcmp(sys_vars->u0, "ZERO_PHASE"))) {
					// ------------------------------------------------
					// Zero Initial Conditions
					// ------------------------------------------------	
					// Initialize the velocity field
					run_data->u[i] = 1.0 / pow(run_data->k[i], sys_vars->ALPHA);

					// Record the phases and amplitudes
					#if defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY)
					run_data->a_n[i]   = 1.0 / pow(run_data->k[i], sys_vars->ALPHA);
					run_data->phi_n[i] = 0.0;
					#endif

					// Initialize the magnetic field
					#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
					run_data->b[i] = 1.0 / pow(run_data->k[i], sys_vars->BETA) * 1e-4;

					// Record the phases and amplitudes
					#if defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY)
					run_data->b_n[i]   = 1.0 / pow(run_data->k[i], sys_vars->BETA) * 1e-4;
					run_data->psi_n[i] = 0.0;
					#endif
					#endif
				}
				else if(!(strcmp(sys_vars->u0, "ZERO"))) {
					// ------------------------------------------------
					// Zero Initial Conditions
					// ------------------------------------------------	
					// Initialize the velocity field
					run_data->u[i] = 0.0 + 0.0 * I;

					// Record the phases and amplitudes
					#if defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY)
					run_data->a_n[i]   = 0.0;
					run_data->phi_n[i] = 0.0;
					#endif

					// Initialize the magnetic field
					#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
					run_data->b[i] = 0.0 + 0.0 * I;

					// Record the phases and amplitudes
					#if defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY)
					run_data->b_n[i]   = 0.0;
					run_data->psi_n[i] = 0.0;
					#endif
					#endif
				}
				else {
					// ------------------------------------------------
					// Default - Pure Random Initial Conditions
					// ------------------------------------------------
					// Get random uniform number
					r1 = (double)rand() / (double)RAND_MAX;
					r3 = (double)rand() / (double)RAND_MAX;
					#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
					r2 = (double)rand() / (double)RAND_MAX;
					r4 = (double)rand() / (double)RAND_MAX;
					#endif

					// Initialize the velocity field
					run_data->u[i] = r1 + r3 * I;

					// Record the phases and amplitudes
					#if defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY)
					run_data->a_n[i]   = r1;
					run_data->phi_n[i] = r3;
					#endif

					// Initialize the magnetic field
					#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
					run_data->b[i] = (r2 + r4) * 1e-2;
					// Record the phases and amplitudes
					#if defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY)
					run_data->b_n[i]   = r2;
					run_data->psi_n[i] = r4;
					#endif
					#endif
				}
			}
			// ------------------------------------------------
			// Get the Elsassar Variables 
			// ------------------------------------------------
			#if defined(__ELSASSAR_MHD) 
			run_data->z_plus[i]  = run_data->u[i] + run_data->b[i];
			run_data->z_minus[i] = run_data->u[i] - run_data->b[i];
			#endif
			// printf("u[%d]:\t%1.16lf\t%1.16lf i\tb[%d]:\t%1.16lf\t%1.16lf i\n", i - 1, creal(run_data->u[i]), cimag(run_data->u[i]),  i - 1, creal(run_data->b[i]), cimag(run_data->b[i]));		
			// printf("a_n[%d]:\t%1.16lf\tphi[%d]:\t%1.16lf\tb_n[%d]:\t%1.16lf\tpsi[%d]:\t%1.16lf\n", i, run_data->a_n[i], i, run_data->phi_n[i], i, run_data->b_n[i], i, run_data->psi_n[i]);
			// printf("a_n[%d]:\t%1.16lf\tphi[%d]:\t%1.16lf\n", i - 1, run_data->a_n[i], i - 1, run_data->phi_n[i]);
		}
		// printf("\n");
    }
}
/**
 * Function to initialize the shell wavenumber array
 * 
 * @param k Array to contain the shell wavenumber
 * @param N int containging the maximum shell level
 */
void InitializeShellWavenumbers(double* k, const long int N) {

	// -------------------------------
	// Define Shell Wavenumbers
	// -------------------------------
	for (int i = 0; i < N + 4; ++i) {
		if (i >= 2 && i < N + 2) {
			k[i] = sys_vars->k_0 * pow(sys_vars->Lambda, (i - 2));
		}
		else {
			k[i] = 0.0;
		}
	}
}
/**
 * Function to initialize the forcing 
 * @param N Number of shells
 */
void InitializeForicing(const long int N, double dt) {

	// Initialize variables
	int n; 
	double tau_0;
	
	// ------------------------------------------------
	// Allocate Memory for Forcing Data
	// ------------------------------------------------
	run_data->forcing_u = (double complex* )malloc(sizeof(double complex) * (N + 4));
	run_data->forcing_b = (double complex* )malloc(sizeof(double complex) * (N + 4));

	// ------------------------------------------------
	// Initialize Forcing Data
	// ------------------------------------------------
	for (int i = 0; i < N + 4; ++i) {
		// Get the proper indx
		n = i + 2;

		// Compute the forcing and intialize
		if(!(strcmp(sys_vars->forcing, "DELTA"))) {
			// ------------------------------------------------
			// Delta function on the zero mode
			// ------------------------------------------------
			// Record the forcing 
			if (i - 1 == sys_vars->force_k) {
				run_data->forcing_u[i] = sys_vars->force_scale_var * (1.0 + 1.0 * I) * my_delta(0.0, 0.0);			
			}
			else {
				run_data->forcing_u[i] = 0.0 + 0.0 * I;			
			}
			#if defined(__MAGNETO)
			run_data->forcing_b[i] = 0.0 + 0.0 * I;			
			#endif
		}
		else if(!(strcmp(sys_vars->forcing, "FXD_AMP"))) {
			// -----------------------------------------------
			// Fixed Amplitude Forcing
			// -----------------------------------------------
			// Initialize variables
			double amp_u_1;
			#if defined(__MAGNETO)
			double amp_b_1;			
			#endif

			// Record the amplitude of the first mode
			if (i == 2) {
				#if defined(PHASE_ONLY)
				amp_u_1                = run_data->a_n[i];
				run_data->forcing_u[i] = run_data->a_n[i] * cexp(I * run_data->phi_n[i]);
				#if defined(__MAGNETO)
				amp_b_1                = run_data->b_n[i];
				run_data->forcing_b[i] = run_data->b_n[i] * cexp(I * run_data->psi_n[i]);
				#endif				
				#else
				amp_u_1                = cabs(run_data->u[i]);
				run_data->forcing_u[i] = run_data->u[i];
				#if defined(__MAGNETO)
				amp_b_1                = cabs(run_data->b[i]);
				run_data->forcing_b[i] = run_data->b[i];
				#endif				
				#endif
			}
			else if (i > 2 && i <= sys_vars->force_k + 1){
				// Set the amplitudes of the fixed modes
				#if defined(PHASE_ONLY)
				run_data->a_n[i]       = amp_u_1;
				run_data->forcing_u[i] = run_data->a_n[i] * cexp(I * run_data->phi_n[i]);
				#if defined(__MAGNETO)
				run_data->b_n[i]       = amp_b_1;
				run_data->forcing_b[i] = run_data->b_n[i] * cexp(I * run_data->psi_n[i]);
				#endif				
				#else
				run_data->u[i]         = amp_u_1;
				run_data->forcing_u[i] = run_data->u[i];
				#if defined(__MAGNETO)
				run_data->b[i]         = amp_b_1;
				run_data->forcing_b[i] = run_data->b[i];
				#endif		
				#endif	
			}
		}
		else if(!(strcmp(sys_vars->forcing, "EXP_STOC"))) {
			// ------------------------------------------------
			// Exponentially Correlated Stochastic Forcing - Satisfying a Langevin Eqn -> \dot{f}_n = -1/\tau_1 f_n + \sigma\zeta; \zeta ~ N(0, 1); \sigma is noise amplitude and \tau_1 is the largest scale eddy turnover time
			// ------------------------------------------------
			// Loop over the forced modes and compute the intial forcing forcing
			if (i >= 2 && i <= sys_vars->force_k + 1){
				// Generate uniform random numbers
				double rand1 = ((double)rand() / (double) RAND_MAX);
				double rand2 = ((double)rand() / (double) RAND_MAX);

				// Compute the forcing timescale
				tau_0 = 1.0 / pow(run_data->k[i] / K_0, 2.0/3.0);

				// Compute the forcing
				run_data->forcing_u[i] = sqrt(- 2.0 * (dt / tau_0) * log10(rand1)) * cexp(I * 2.0 * M_PI * rand2);
				#if defined(__MAGNETO)
				run_data->forcing_b[i] = 0.0 + 0.0 * I;
				#endif
			}
			else {
				run_data->forcing_u[i] = 0.0 + 0.0 * I;	
				#if defined(__MAGNETO)
				run_data->forcing_b[i] = 0.0 + 0.0 * I;
				#endif
			}
		}
		else if(!(strcmp(sys_vars->forcing, "DELTA_STOC"))) {
			// ------------------------------------------------
			// Delta Correlated Stochastic Forcing - Satisfying a Langevin Eqn -> \dot{f}_n = -1/\tau_1 f_n + \sigma\zeta; \zeta ~ N(0, 1); \sigma is noise amplitude and \tau_1 is the largest scale eddy turnover time
			// ------------------------------------------------
			// Loop over the forced modes and compute the intial forcing forcing
			if (i >= 2 && i <= sys_vars->force_k + 1){
				// Generate uniform random numbers
				double rand1 = ((double)rand() / (double) RAND_MAX);
				double rand2 = ((double)rand() / (double) RAND_MAX);

				// Compute the forcing timescale
				tau_0 = 1.0 / pow(run_data->k[i] / K_0, 2.0/3.0);

				// Compute the forcing
				run_data->forcing_u[i] = FORC_STOC_SIGMA * sqrt(-2.0 * (dt / tau_0) * log10(rand1)) * cexp(I * 2.0 * M_PI * rand2);
				#if defined(__MAGNETO)
				run_data->forcing_b[i] = 0.0 + 0.0 * I;
				#endif
			}
			else {
				run_data->forcing_u[i] = 0.0 + 0.0 * I;	
				#if defined(__MAGNETO)
				run_data->forcing_b[i] = 0.0 + 0.0 * I;
				#endif
			}
		}
		else if(!(strcmp(sys_vars->forcing, "NONE"))) {
			// ------------------------------------------------
			// No Forcing
			// ------------------------------------------------
			run_data->forcing_u[i] = 0.0 + 0.0 * I;
			#if defined(__MAGNETO)
			run_data->forcing_b[i] = 0.0 + 0.0 * I;
			#endif
		}
		else {
			if (i == 0) {
				// Print warning to screen that no valid forcing (incl. NONE) was selected
				printf("\n["MAGENTA"WARNING"RESET"] --- No valid forcing was selected!!!\n");
			}
		}
		// printf("f_u[%d]:\t%1.16lf\t%1.16lf i\t\tf_b[%d]:\t%1.16lf\t%1.16lf i\n", i - 1, creal(run_data->forcing_u[i]), creal(run_data->forcing_u[i]), i - 1, creal(run_data->forcing_b[i]), creal(run_data->forcing_b[i]));
		// printf("u[%d]:\t%1.16lf\t%1.16lf i\t\tf_u[%d]:\t%1.16lf\t%1.16lf i\n", i - 1, creal(run_data->forcing_u[i]), creal(run_data->forcing_u[i]), i - 1, creal(run_data->forcing_b[i]), creal(run_data->forcing_b[i]));
	}
	// printf("\n");
}
/**
 * Function compute the forcing for the current timestep
 * @param N Number of shells
 */
void ComputeForcing(double dt, const long int N) {

	// Initialize varaibles
	double rand1, rand2;
	double tau_0, exp_fac;

	// ------------------------------------------------
	// Exponentially Correlated Stochastic Forcing
	// ------------------------------------------------
 	if(!(strcmp(sys_vars->forcing, "EXP_STOC"))) {
		// Loop over the forced modes and compute the forcing -> note shell n = 1 is at index 2
		for (int i = 2; i <= sys_vars->force_k + 1; ++i) {
 			// Generate two uniform random numbers
			rand1 = ((double)rand() / (double) RAND_MAX);
			rand2 = ((double)rand() / (double) RAND_MAX);

			// Compute the forcing timescale
			tau_0 = 1.0 / pow(run_data->k[i] / K_0, 2.0/3.0);

			// Compute the exponential prefactor
			exp_fac = cexp(-dt / tau_0);

			// Compute forcing
			run_data->forcing_u[i] = sys_vars->force_scale_var * (exp_fac * run_data->forcing_u[i] + sqrt(-2.0 * (1.0 - pow(exp_fac, 2.0)) * log10(rand1)) * cexp(I * 2.0 * M_PI * rand2));
		}
	}
	else if(!(strcmp(sys_vars->forcing, "DELTA_STOC"))) {
		// Loop over the forced modes and compute the forcing -> note shell n = 1 is at index 2
		for (int i = 2; i <= sys_vars->force_k + 1; ++i) {
 			// Generate two uniform random numbers
			rand1 = ((double)rand() / (double) RAND_MAX);
			rand2 = ((double)rand() / (double) RAND_MAX);

			// Compute the forcing timescale
			tau_0 = 1.0 / pow(run_data->k[i] / K_0, 2.0/3.0);

			// printf("|u|: %lf \t|k|: %lf \t tau: %lf\n", cabs(run_data->u[i]), cabs(run_data->k[i]), tau_0);

			// Compute forcing 
			run_data->forcing_u[i] = run_data->forcing_u[i] - dt * run_data->forcing_u[i] / tau_0 + sys_vars->force_scale_var * sqrt(-2.0 * (dt / tau_0) * log10(rand1)) * cexp(I * 2.0 * M_PI * rand2);
			// run_data->forcing_u[i] = run_data->forcing_u[i] - dt * run_data->forcing_u[i] / tau_0 + FORC_STOC_SIGMA * sqrt(-2.0 * (dt / tau_0) * log10(rand1)) * cexp(I * 2.0 * M_PI * rand2);
			// run_data->forcing_u[i] = run_data->forcing_u[i] + sys_vars->force_scale_var * (- run_data->forcing_u[i] / tau_0 + FORC_STOC_SIGMA * sqrt(-2.0 * (dt / tau_0) * log10(rand1)) * cexp(I * 2.0 * M_PI * rand2));
		}
	}
}
/**
 * Wrapper function to add the forcing term to the nonlinear term
 */
void AddForcing(double complex* u_nonlin, double complex* b_nonlin) {

	// Initialize variables
	const int N = sys_vars->N;
	int n;

	for (int i = 0; i < N; ++i) {
		// Get the propper index
		n = i + 2;

		// ------------------------------------------------
		// Add Forcing
		// ------------------------------------------------
		//------------------------------------------------------- Fixed First k Amplitude Forcing
		if(!(strcmp(sys_vars->forcing, "FXD_AMP"))) {
			// If fixed amplitude forcing reset the amplitudes to the forced modes
			if (i >= 2 && i <= sys_vars->force_k + 1) {
				u_nonlin[i] *= cabs(run_data->forcing_u[i]) / cabs(u_nonlin[i]);
				#if defined(__MAGNETO)
				// Add forcing here for the magnetic field
				b_nonlin[i] *= cabs(run_data->forcing_b[i]) / cabs(b_nonlin[i]); 
				#endif
			}
		}
		//------------------------------------------------------- Add Stochastic forcing and also update the forcing for the next iteration
		else if(!(strcmp(sys_vars->forcing, "EXP_STOC")) || !(strcmp(sys_vars->forcing, "DELTA_STOC"))) {
			if (n >= 2 && n <= sys_vars->force_k + 1){
				// Add the forcing 
				u_nonlin[n] += run_data->forcing_u[n];
			}
		}
		//------------------------------------------------------- Else if normal forcing is selected just add the forcing to the forced modes
		else {
			u_nonlin[n] += run_data->forcing_u[n];
			#if defined(__MAGNETO)
			// Add forcing here for the magnetic field
			b_nonlin[n] += run_data->forcing_b[n]; 
			#endif
		}
	// if (n >= 2 && n <= N - 1) {
	// 	printf("u[%d]: %lf %lf\ta_n[%d]: %lf\tf[%d]: %lf %lf\tf_n[%d]: %lf\n", n, creal(run_data->u[n]), cimag(run_data->u[n]), n, cabs(run_data->u[n]), n, creal(run_data->forcing_u[n]), cimag(run_data->forcing_u[n]), n, cabs(run_data->forcing_u[n]));
	// }
	}
	// for (n = 0; n < N + 4; ++n) {
	// 	printf("u[%d]: %lf %lf\ta_n[%d]: %lf\tf[%d]: %lf %lf\tf_n[%d]: %lf\n", n, creal(run_data->u[n]), cimag(run_data->u[n]), n, cabs(run_data->u[n]), n, creal(run_data->forcing_u[n]), cimag(run_data->forcing_u[n]), n, cabs(run_data->forcing_u[n]));
	// }
	// printf("\n");
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
	sys_vars->min_dt = MIN_STEP_SIZE;
	sys_vars->max_dt = 10;
	double tmp_num_t_steps;

	// -------------------------------
	// Integration Counters
	// -------------------------------
	// Number of time steps and saving steps
	tmp_num_t_steps = ((*T) - (*t0)) / (*dt);
	sys_vars->num_t_steps = (long int)round(tmp_num_t_steps);
	if (sys_vars->TRANS_ITERS_FLAG == TRANSIENT_ITERS) {
		// Get the transient iterations
		(* trans_steps)       = (long int)round(sys_vars->TRANS_ITERS_FRAC * tmp_num_t_steps);
		sys_vars->trans_iters = (* trans_steps);
		sys_vars->trans_time  = (* trans_steps) * (* dt);

		// Get the number of steps to perform before printing to file -> allowing for a transient fraction of these to be ignored
		sys_vars->num_print_steps = (sys_vars->num_t_steps >= sys_vars->SAVE_EVERY ) ? (sys_vars->num_t_steps - sys_vars->trans_iters) / sys_vars->SAVE_EVERY  + 1: sys_vars->num_t_steps - sys_vars->trans_iters + 1;	 
		printf("Total Iters: %ld\t Saving Iters: %ld\t Transient Steps: %ld\n", sys_vars->num_t_steps, sys_vars->num_print_steps, sys_vars->trans_iters);
	}
	else {
		// Get the transient iterations 
		(* trans_steps)       = 0;
		sys_vars->trans_iters = (* trans_steps);
		sys_vars->trans_time  = (* trans_steps) * (* dt);


		// Get the number of steps to perform before printing to file
		sys_vars->num_print_steps = (sys_vars->num_t_steps >= sys_vars->SAVE_EVERY ) ? sys_vars->num_t_steps / sys_vars->SAVE_EVERY + 1: sys_vars->num_t_steps + 1; // plus one to include initial condition
		printf("Total Iters: %ld\t Saving Iters: %ld\n", sys_vars->num_t_steps, sys_vars->num_print_steps);
		// printf("\n\nT:%1.16lf,t0:%1.16lf,dt:%1.16lf--numt:%1.16lf\tround(numt):%1.16lf\ttrans_frac:%1.16lf\ttans:%1.16lf\ttrans:%ld NUm_print: %ld\n\n", (*T), (*t0), (*dt), ((*T) - (*t0)) / (*dt), round(((*T) - (*t0)) / (*dt)), sys_vars->TRANS_ITERS_FRAC, (sys_vars->TRANS_ITERS_FRAC * (double)tmp_num_t_steps), (*trans_steps), sys_vars->num_print_steps);
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
void PrintUpdateToTerminal(long int iters, double t, double dt, double T, long int save_data_indx) {

	// Initialize variables
	
	if (iters < sys_vars->trans_iters) {
		#if defined(__MAGNETO)
		printf("Iter: %1.3g / %1.1g\tt: %1.6lf / %1.3lf\ttau_0: %1.2lf\tdt: %1.6g\tKE: %6.6g\tHEL_U: %6.6g\tHEL_B: %6.6g\tX-HEL: %6.6g\n", (double)iters, (double)sys_vars->num_t_steps, t, T, sys_vars->eddy_turnover_time, dt, run_data->tot_energy[save_data_indx], run_data->tot_hel_u[save_data_indx], run_data->tot_hel_b[save_data_indx], run_data->tot_cross_hel[save_data_indx]);
		#else
		printf("Iter: %1.3g / %1.1g\tt: %1.6lf / %1.3lf\ttau_0: %1.2lf\tdt: %1.6g\tKE: %6.6g\tHEL: %6.6g\t\n", (double)iters, (double)sys_vars->num_t_steps, t, T, sys_vars->eddy_turnover_time, dt, run_data->tot_energy[save_data_indx], run_data->tot_hel_u[save_data_indx]);
		#endif
	}
	else {
		#if defined(__MAGNETO)
		printf("Iter: %1.3g / %1.1g\tt: %1.6lf / %1.3lf\ttau: %1.2lftau_0 / %1.2lftau_0\ttau_0: %1.2lf\tdt: %1.6g\tKE: %6.6g\tHEL_U: %6.6g\tHEL_B: %6.6g\tX-HEL: %6.6g\n", (double)iters, (double)sys_vars->num_t_steps, t, T, (t - sys_vars->trans_time)/sys_vars->eddy_turnover_time, (T - sys_vars->trans_time)/sys_vars->eddy_turnover_time, sys_vars->eddy_turnover_time, dt, run_data->tot_energy[save_data_indx], run_data->tot_hel_u[save_data_indx], run_data->tot_hel_b[save_data_indx], run_data->tot_cross_hel[save_data_indx]);
		#else
		printf("Iter: %1.3g / %1.1g\tt: %1.6lf / %1.3lf\ttau: %1.2lftau_0 / %1.2lftau_0\ttau_0: %1.2lf\tdt: %1.6g\tKE: %6.6g\tHEL: %6.6g\t\n", (double)iters, (double)sys_vars->num_t_steps, t, T, (t - sys_vars->trans_time)/sys_vars->eddy_turnover_time, (T - sys_vars->trans_time)/sys_vars->eddy_turnover_time, sys_vars->eddy_turnover_time, dt, run_data->tot_energy[save_data_indx], run_data->tot_hel_u[save_data_indx]);
		#endif
	}
}
/**
 * Function that checks the system to see if it is ok to continue integrations. Checks for blow up, timestep and iteration limits etc
 * @param dt    		 The updated timestep for the next iteration
 * @param iters 		 The number of iterations for the next iteration
 * @param save_data_indx The current index for saving data to
 */
void SystemCheck(double dt, long int iters, long int save_data_indx) {

	// Initialize variables
	int n;
	double max_vel = 0.0;
	double max_mag = 0.0;
	
	// -------------------------------
	// Check Stopping Criteria 
	// -------------------------------	
	// Get the max field value
	for (int i = 0; i < sys_vars	->N; ++i) {
		// Get temp index
		n = i + 2;

		// Check for max
		#if defined(PHASE_ONLY)
		max_vel = fmax(max_vel, cabs(run_data->a_n[n]));
		#if defined(__MAGNETO)
		max_mag = fmax(max_vel, cabs(run_data->b_n[n]));
		#endif
		#else
		max_vel = fmax(max_vel, cabs(run_data->u[n]));
		#if defined(__MAGNETO)
		max_mag = fmax(max_vel, cabs(run_data->b[n]));
		#endif
		#endif
	}

	// -------------------------------
	// Check Stopping Criteria 
	// -------------------------------
	if (dt <= MIN_STEP_SIZE) {
		// Print error message to error stream
		fprintf(stderr, "\n["YELLOW"SOVLER FAILURE"RESET"]--- Timestep has become too small to continue at Iter: ["CYAN"%ld"RESET"]\n-->> Exiting!!!\n", iters);

		// Write current state of the system to file
		FinalWriteAndCloseOutputFile(sys_vars->N, iters, save_data_indx);

		// Exit program
		exit(1);		
	}
	else if (iters >= MAX_ITERS) {
		// Print error message to error stream
		fprintf(stderr, "\n["YELLOW"SOVLER FAILURE"RESET"]--- The maximum number of iterations has been reached at Iter: ["CYAN"%ld"RESET"]\n-->> Exiting!!!\n", iters);

		// Write current state of the system to file
		FinalWriteAndCloseOutputFile(sys_vars->N, iters, save_data_indx);

		// Exit program
		exit(1);		
	}
	else if ((max_vel >= MAX_FIELD_LIM) || (max_mag >= MAX_FIELD_LIM)) {
		// Print error message to error stream
		fprintf(stderr, "\n["YELLOW"SOVLER FAILURE"RESET"]--- The maximum field value has been reached at Iter: ["CYAN"%ld"RESET"]\n-->> Exiting!!!\n", iters);

		// Write current state of the system to file
		FinalWriteAndCloseOutputFile(sys_vars->N, iters, save_data_indx);

		// Exit program
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
	run_data->k = (double* )malloc(sizeof(double) * (N + 4));  // k
	if (run_data->k == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Shell Wavenumber List");
		exit(1);
	}

	// ----------------------------------
	// Allocate Velocity Field Variables
	// ----------------------------------
	// Full velocity field
	run_data->u = (double complex* ) malloc(sizeof(double complex) * (N + 4));
	if (run_data->u == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Velocity Field");
		exit(1);
	}
	#if defined(PHASE_ONLY_FXD_AMP)	|| defined(PHASE_ONLY)
	// The Fourier amplitudes
	run_data->a_n = (double* ) malloc(sizeof(double) * (N + 4));
	if (run_data->a_n == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Velocity Field Amplitude");
		exit(1);
	}
	// The Fourier phases
	run_data->phi_n = (double* ) malloc(sizeof(double) * (N + 4));
	if (run_data->phi_n == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Velocity Field Phase");
		exit(1);
	}
	#endif

	// ---------------------------------
	// Allocate Magnetic Field Variables
	// ---------------------------------
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	// Full magnetic field
	run_data->b = (double complex* ) malloc(sizeof(double complex) * (N + 4));
	if (run_data->b == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Field");
		exit(1);
	}	
	#if defined(__ELSASSAR_MHD)
	// Elsassar Variables
	run_data->z_plus = (double complex* ) malloc(sizeof(double complex) * (N + 4));
	if (run_data->z_plus == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Elsassar Z+");
		exit(1);
	}	
	run_data->z_minus = (double complex* ) malloc(sizeof(double complex) * (N + 4));
	if (run_data->z_minus == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Elsassar Z-");
		exit(1);
	}	
	#endif
	#if defined(PHASE_ONLY_FXD_AMP)	|| defined(PHASE_ONLY)
	// The Fourier amplitudes
	run_data->b_n = (double* ) malloc(sizeof(double) * (N + 4));
	if (run_data->b_n == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Field Amplitude");
		exit(1);
	}
	// The Fourier phases
	run_data->psi_n = (double* ) malloc(sizeof(double) * (N + 4));
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
	RK_data->RK1_u       = (double complex* )malloc(sizeof(double complex) * (N + 4));
	RK_data->RK2_u       = (double complex* )malloc(sizeof(double complex) * (N + 4));
	RK_data->RK3_u       = (double complex* )malloc(sizeof(double complex) * (N + 4));
	RK_data->RK4_u       = (double complex* )malloc(sizeof(double complex) * (N + 4));
	RK_data->RK_u_tmp    = (double complex* )malloc(sizeof(double complex) * (N + 4));
	#if defined(AB4CN)
	RK_data->AB_tmp_u	 = (double complex* )malloc(sizeof(double complex) * (N + 4));
	for (int i = 0; i < 3; ++i) {
		RK_data->AB_tmp_nonlin_u[i] = (double complex* )malloc(sizeof(double complex) * (N + 4));
	}
	#endif
	#if defined(__MAGNETO)
	RK_data->RK1_b       = (double complex* )malloc(sizeof(double complex) * (N + 4));
	RK_data->RK2_b       = (double complex* )malloc(sizeof(double complex) * (N + 4));
	RK_data->RK3_b       = (double complex* )malloc(sizeof(double complex) * (N + 4));
	RK_data->RK4_b       = (double complex* )malloc(sizeof(double complex) * (N + 4));
	RK_data->RK_b_tmp    = (double complex* )malloc(sizeof(double complex) * (N + 4));
	#if defined(AB4CN)
	RK_data->AB_tmp_b	 = (double complex* )malloc(sizeof(double complex) * (N + 4));
	for (int i = 0; i < 3; ++i) {
		RK_data->AB_tmp_nonlin_b[i] = (double complex* )malloc(sizeof(double complex) * (N + 4));
	}
	#endif
	#endif
	// #endif
	if (RK_data->RK1_u == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK1 Velocity");
		exit(1);
	}
	if (RK_data->RK2_u == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK2 Velocity");
		exit(1);
	}
	if (RK_data->RK3_u == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK3 Velocity");
		exit(1);
	}
	if (RK_data->RK4_u == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK4 Velocity");
		exit(1);
	}
	if (RK_data->RK_u_tmp == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK_u_tmp");
		exit(1);
	}

	#if defined(AB4CN)
	if (RK_data->AB_tmp_u == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "AB_tmp_u");
		exit(1);
	}
	for (int i = 0; i < 3; ++i) {
		if (RK_data->AB_tmp_nonlin_u[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "AB_tmp_nonlin_u");
			exit(1);
		}
	}
	#endif
	#if defined(__MAGNETO)
	if (RK_data->RK1_b == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK1 Magnetic");
		exit(1);
	}
	if (RK_data->RK2_b == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK2 Magnetic");
		exit(1);
	}
	if (RK_data->RK3_b == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK3 Magnetic");
		exit(1);
	}
	if (RK_data->RK4_b == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK4 Magnetic");
		exit(1);
	}
	if (RK_data->RK_b_tmp == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "RK_b_tmp");
		exit(1);
	}
	#if defined(AB4CN)
	if (RK_data->AB_tmp_b == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "AB_tmp_b");
		exit(1);
	}
	for (int i = 0; i < 3; ++i) {
		if (RK_data->AB_tmp_nonlin_b[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Integration Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "AB_tmp_nonlin_b");
			exit(1);
		}
	}
	#endif
	#endif

	// -------------------------------
	// Initialize All Data 
	// -------------------------------
	for (int i = 0; i < N + 4; ++i) {
		// Initialize the arrays
		run_data->k[i]       = 0;
		run_data->u[i]       = 0.0 + 0.0 * I;
		RK_data->RK1_u[i]    = 0.0 + 0.0 * I;
		RK_data->RK2_u[i]    = 0.0 + 0.0 * I;
		RK_data->RK3_u[i]    = 0.0 + 0.0 * I;
		RK_data->RK4_u[i]    = 0.0 + 0.0 * I;
		RK_data->RK_u_tmp[i] = 0.0 + 0.0 * I;
		#if defined(AB4CN)
		RK_data->AB_tmp_u[i] = 0.0 + 0.0 * I;
		for (int j = 0; j < 3; ++j) {
			RK_data->AB_tmp_nonlin_u[j][i] = 0.0 + 0.0 * I;
		}
		#endif
		#if defined(PHASE_ONLY_FXD_AMP)
		run_data->a_n[i]   = 0.0;
		run_data->phi_n[i] = 0.0;
		#endif
		#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
		run_data->b[i]       = 0.0 + 0.0 * I;
		RK_data->RK1_b[i]    = 0.0 + 0.0 * I;
		RK_data->RK2_b[i]    = 0.0 + 0.0 * I;
		RK_data->RK3_b[i]    = 0.0 + 0.0 * I;
		RK_data->RK4_b[i]    = 0.0 + 0.0 * I;
		RK_data->RK_b_tmp[i] = 0.0 + 0.0 * I;
		#if defined(__ELSASSAR_MHD)
		run_data->z_plus[i]  = 0.0 + 0.0 * I;
		run_data->z_minus[i] = 0.0 + 0.0 * I;
		#endif
		#if defined(AB4CN)
		RK_data->AB_tmp_b[i] = 0.0 + 0.0 * I;
		for (int j = 0; j < 3; ++j) {
			RK_data->AB_tmp_nonlin_b[j][i] = 0.0 + 0.0 * I;
		}
		#endif
		#if defined(PHASE_ONLY_FXD_AMP)
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
	free(run_data->k);

	// Free system variables
	free(run_data->u);
	#if defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY)
	free(run_data->a_n);
	free(run_data->phi_n);
	#endif
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	free(run_data->b);
	#if defined(__ELSASSAR_MHD)
	free(run_data->z_plus);
	free(run_data->z_minus);
	#endif
	#if defined(PHASE_ONLY_FXD_AMP)
	free(run_data->b_n);
	free(run_data->psi_n);
	#endif
	#endif

	// ------------------------
	// Free System Msr Memory 
	// ------------------------
	#if defined(__SYS_MEASURES)
	free(run_data->tot_energy);
	free(run_data->tot_hel_u);
	free(run_data->tot_diss_u);
	#if defined(__MAGNETO)
	free(run_data->tot_hel_b);
	free(run_data->tot_cross_hel);
	free(run_data->tot_diss_b);
	#endif
	free(run_data->u_charact);
	free(run_data->int_scale);
	free(run_data->taylor_micro_scale);
	free(run_data->reynolds_no);
	free(run_data->kolmogorov_scale);
	#endif
	#if defined(__ENRG_SPECT) 
	free(run_data->energy_spect);
	#endif
	#if defined(__DISS_SPECT)
	free(run_data->diss_spect);
	#endif
	#if defined(__ENRG_SPECT_AVG) 
	free(run_data->energy_spect_t_avg);
	#endif
	#if defined(__DISS_SPECT_AVG)
	free(run_data->diss_spect_t_avg);
	#endif
	#if defined(__TIME)
	free(run_data->time);
	#endif
	#if defined(__ENRG_FLUX)
	free(run_data->energy_flux);
	free(run_data->energy_diss_u);
	free(run_data->energy_input_u);
	#if defined(__MAGNETO)
	free(run_data->energy_diss_b);
	free(run_data->energy_input_b);
	#endif
	#endif
	#if defined(__TOT_ENRG_FLUX)
	free(run_data->tot_energy_flux);
	free(run_data->tot_energy_diss_u);
	free(run_data->tot_energy_input_u);
	#if defined(__MAGNETO)
	free(run_data->tot_energy_diss_b);
	free(run_data->tot_energy_input_b);
	#endif
	#endif
	#if defined(__ENRG_FLUX_AVG)
	free(run_data->energy_flux_t_avg);
	#endif
	free(run_data->forcing_u);
	free(run_data->forcing_b);

	// ------------------------
	// Free Stats Objects
	// ------------------------
	// Free stats objects
	#if defined(STATS)
	for (int i = 0; i < NUM_POW; ++i) {
		#if defined(__STR_FUNC_VEL)
		free(stats_data->vel_str_func[i]);
		#endif
		#if defined(__STR_FUNC_VEL_FLUX)
		for (int j = 0; j < 2; ++j) {
			free(stats_data->vel_flux_str_func[j][i]);
			free(stats_data->vel_flux_str_func_abs[j][i]);
		}		
		#endif
		#if defined(__MAGNETO)
		#if defined(__STR_FUNC_MAG)
		free(stats_data->mag_str_func[i]);
		#endif
		#if defined(__STR_FUNC_MAG_FLUX)
		for (int j = 0; j < 2; ++j) {
			free(stats_data->mag_flux_str_func[j][i]);
			free(stats_data->mag_flux_str_func_abs[j][i]);
		}
		#endif
		#endif
	}
	for (int i = 0; i < sys_vars->N; ++i) {
		gsl_rstat_free(stats_data->real_vel_moments[i]);
		#if defined(__MAGNETO)
		gsl_rstat_free(stats_data->real_mag_moments[i]);
		#endif
		#if defined(__VEL_HIST)
		gsl_histogram_free(stats_data->real_vel_hist[i]);
		#endif
		#if defined(__MAG_HIST)
		gsl_histogram_free(stats_data->real_mag_hist[i]);
		#endif
	}
	free(stats_data->real_vel_moments);
	#if defined(__MAGNETO)
	free(stats_data->real_mag_moments);
	#endif
	#if defined(__VEL_HIST)
	free(stats_data->real_vel_hist);
	#endif
	#if defined(__MAG_HIST)
	free(stats_data->real_mag_hist);
	#endif
	#endif

	// ------------------------
	// Free Phase Sync Memory 
	// ------------------------
	#if defined(__CONSERVED_PHASES) || defined(__PHASE_SYNC) || defined(__PHASE_SYNC_STATS)
	FreePhaseSyncObjects();
	#endif

	// ----------------------------
	// Free Integration Variables
	// ----------------------------
	// Free integration variables
	free(RK_data->RK1_u);
	free(RK_data->RK2_u);
	free(RK_data->RK3_u);
	free(RK_data->RK4_u);
	free(RK_data->RK_u_tmp);
	#if defined(AB4CN)
	free(RK_data->AB_tmp_u);
	for (int i = 0; i < 3; ++i)	{
		free(RK_data->AB_tmp_nonlin_u[i]);
	}
	#endif
	#if defined(__MAGNETO)
	free(RK_data->RK1_b);
	free(RK_data->RK2_b);
	free(RK_data->RK3_b);
	free(RK_data->RK4_b);
	free(RK_data->RK_b_tmp);
	#if defined(AB4CN)
	free(RK_data->AB_tmp_b);
	for (int i = 0; i < 3; ++i)	{
		free(RK_data->AB_tmp_nonlin_b[i]);
	}
	#endif
	#endif
}
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------