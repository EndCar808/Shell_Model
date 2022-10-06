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
	// // If in testing / debug mode - create/ open test file
	// #if defined(DEBUG)
	// OpenTestingFile();
	// #endif

	// Initialize the collocation points and wavenumber space 
	InitializeShellWavenumbers(run_data->k, N);

	// Get initial conditions - seed for random number generator is set here
	InitialConditions(N);

	// Initialize the forcing
	InitializeForicing(N);


	// Initialize stats objects if required
	#if defined(STATS)
	InitializeStats();
	#endif

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
		// Write To File
		// -------------------------------
		if (iters % sys_vars->SAVE_EVERY == 0) {

			// Record System Measurables
			ComputeSystemMeasurables(t, save_data_indx, RK_data);

			// Compute stats
			#if defined(STATS)
			ComputeStats();
			#endif

			// If and when transient steps are complete write to file
			if (iters > trans_steps) {
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
#if defined(INT_FAC_RK4)
void IntFacRK4Step(const double dt, const long int N, RK_data_struct* RK_data) {

	// Initialize vairables
	int n;
	#if defined(PHASE_ONLY_DIRECT)
	double int_fac_u;
	double int_fac_u_1;
	#if defined(__MAGNETO)
	double int_fac_b;
	double int_fac_b_1;
	#endif
	#else
	fftw_complex int_fac_u;
	fftw_complex int_fac_u_1;
	#if defined(__MAGNETO)
	fftw_complex int_fac_b;
	fftw_complex int_fac_b_1;
	#endif
	#endif
	#if defined(PHASE_ONLY)
	double po_norm_fac_u;
	#if defined(__MAGNETO)
	double po_norm_fac_b;
	#endif
	#endif

	///----------------------- Forcing
	// Compute the forcing for the current iteration
	ComputeForicing(N);

	/////////////////////
	/// RK STAGES
	/////////////////////
	// ----------------------- Stage 1
	#if defined(PHASE_ONLY_DIRECT)
	NonlinearTerm(run_data->phi_n, run_data->psi_n, RK_data->RK1_u, RK_data->RK1_b, N);
	#else
	NonlinearTerm(run_data->u, run_data->b, RK_data->RK1_u, RK_data->RK1_b, N);
	#endif
	for (int i = 0; i < N; ++i) {
		// Get proper indx
		n = i + 2;


		#if defined(PHASE_ONLY_DIRECT)
		// Get the integrating factor
		int_fac_u = exp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] / 2.0 + cabs(run_data->forcing_u[n]) * sin(carg(run_data->forcing_u[n]) - run_data->phi_n[n]));

		// Update temorary velocity term
		RK_data->RK_u_tmp[n] = run_data->phi_n[n] * int_fac_u + dt * RK4_A21 * RK_data->RK1_u[n] * int_fac_u;
		#if defined(__MAGNETO)
		// Get the integrating factor
		int_fac_b = exp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] / 2.0 + cabs(run_data->forcing_b[n]) * sin(carg(run_data->forcing_b[n]) - run_data->psi_n[n]));

		// Update temporary magnetic term
		RK_data->RK_b_tmp[n] = run_data->psi_n[n] * int_fac_b + dt * RK4_A21 * RK_data->RK1_b[n] * int_fac_b;		
		#endif
		#else
		// Get the integrating factor
		int_fac_u = exp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] / 2.0 + run_data->forcing_u[n]);

		// Update temorary velocity term
		RK_data->RK_u_tmp[n] = run_data->u[n] * int_fac_u + dt * RK4_A21 * RK_data->RK1_u[n] * int_fac_u;
		#if defined(__MAGNETO)
		// Get the integrating factor
		int_fac_b = exp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] / 2.0 + run_data->forcing_b[n]);
		
		// Update temporary magnetic term
		RK_data->RK_b_tmp[n] = run_data->b[n] * int_fac_b + dt * RK4_A21 * RK_data->RK1_b[n] * int_fac_b;
		#endif
		#endif
	}

	// ----------------------- Stage 2
	NonlinearTerm(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK2_u, RK_data->RK2_b, N);
	for (int i = 0; i < N; ++i) {
		// Get proper indx
		n = i + 2;

		#if defined(PHASE_ONLY_DIRECT)
		// Get the integrating factor
		int_fac_u = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] / 2.0 + cabs(run_data->forcing_u[n]) * sin(carg(run_data->forcing_u[n]) - run_data->phi_n[n]));

		// Update temorary velocity term
		RK_data->RK_u_tmp[n] = run_data->phi_n[n] * int_fac_u + dt * RK4_A32 * RK_data->RK2_u[n];
		#if defined(__MAGNETO)
		// Get the integrating factor
		int_fac_b = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] / 2.0 + cabs(run_data->forcing_b[n]) * sin(carg(run_data->forcing_b[n]) - run_data->psi_n[n]));

		// Update temporary magnetic term
		RK_data->RK_b_tmp[n] = run_data->psi_n[n] * int_fac_b + dt * RK4_A32 * RK_data->RK2_b[n];		
		#endif
		#else
		// Get the integrating factor
		int_fac_u = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] / 2.0 + run_data->forcing_u[n]);

		// Update temorary velocity term
		RK_data->RK_u_tmp[n] = run_data->u[n] * int_fac_u + dt * RK4_A32 * RK_data->RK2_u[n];
		#if defined(__MAGNETO)
		// Get the integrating factor
		int_fac_b = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] / 2.0 + run_data->forcing_b[n]);

		// Update temporary magnetic term
		RK_data->RK_b_tmp[n] = run_data->b[n] * int_fac_b + dt * RK4_A32 * RK_data->RK2_b[n];		
		#endif
		#endif
	}

	// ----------------------- Stage 3
	NonlinearTerm(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK3_u, RK_data->RK3_b, N);
	for (int i = 0; i < N; ++i) {
		// Get proper indx
		n = i + 2;

		#if defined(PHASE_ONLY_DIRECT)
		// Get the integrating factor
		int_fac_u   = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] + cabs(run_data->forcing_u[n]) * sin(carg(run_data->forcing_u[n]) - run_data->phi_n[n]));
		int_fac_u_1 = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] / 2.0 + cabs(run_data->forcing_u[n]) * sin(carg(run_data->forcing_u[n]) - run_data->phi_n[n]));

		// Update temorary velocity term
		RK_data->RK_u_tmp[n] = run_data->phi_n[n] * int_fac_u + dt * RK4_A43 * RK_data->RK3_u[n] * int_fac_u_1;
		#if defined(__MAGNETO)
		// Get the integrating factor
		int_fac_b   = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] + cabs(run_data->forcing_b[n]) * sin(carg(run_data->forcing_b[n]) - run_data->psi_n[n]));
		int_fac_b_1 = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] / 2.0 + cabs(run_data->forcing_b[n]) * sin(carg(run_data->forcing_b[n]) - run_data->psi_n[n]));

		// Update temporary magnetic term
		RK_data->RK_b_tmp[n] = run_data->psi_n[n] * int_fac_b + dt * RK4_A43 * RK_data->RK3_b[n] * int_fac_b_1;		
		#endif
		#else
		// Get the integrating factor
		int_fac_u   = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] + run_data->forcing_u[n]);
		int_fac_u_1 = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] / 2.0 + run_data->forcing_u[n]);

		// Update temorary velocity term
		RK_data->RK_u_tmp[n] = run_data->u[n] * int_fac_u + dt * RK4_A43 * RK_data->RK3_u[n] * int_fac_u_1;
		#if defined(__MAGNETO)
		// Get the integrating factor
		int_fac_b   = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] + run_data->forcing_b[n]);
		int_fac_b_1 = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] / 2.0 + run_data->forcing_b[n]);
		
		// Update temporary magnetic term
		RK_data->RK_b_tmp[n] = run_data->b[n] * int_fac_b + dt * RK4_A43 * RK_data->RK3_b[n] * int_fac_b_1;		
		#endif
		#endif
	}

	// ----------------------- Stage 4
	NonlinearTerm(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK4_u, RK_data->RK4_b, N);


	/////////////////////
	/// UPDATE STEP
	/////////////////////
	for (int i = 0; i < N; ++i) {
		// Get temporary index
		n = i + 2;

		#if defined(PHASE_ONLY)
		// Pre-record the amplitudes for resetting after update step
		po_norm_fac_u = cabs(run_data->u[n]);
		#if defined(__MAGNETO)
		po_norm_fac_b = cabs(run_data->b[n]);
		#endif
		#endif

		///-------------------- Update Step
		#if defined(PHASE_ONLY_DIRECT)
		// Get the integrating factors
		int_fac_u   = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] + cabs(run_data->forcing_u[n]) * sin(carg(run_data->forcing_u[n]) - run_data->phi_n[n]));
		int_fac_u_1 = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] / 2.0 + cabs(run_data->forcing_u[n]) * sin(carg(run_data->forcing_u[n]) - run_data->phi_n[n]));

		// Update the new velocity field
		run_data->phi_n[n] = int_fac_u * run_data->phi_n[n] + dt * RK4_B1 * int_fac_u * RK_data->RK1_u[n] + dt * RK4_B2 * int_fac_u_1 * RK_data->RK2_u[n] + dt * RK4_B3 * int_fac_u_1 * RK_data->RK3_u[n] + dt * RK4_B4 * RK_data->RK4_u[n];
		#if defined(__MAGNETO)
		// Get the integrating factors
		int_fac_b   = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] + cabs(run_data->forcing_b[n]) * sin(carg(run_data->forcing_b[n]) - run_data->psi_n[n]));
		int_fac_b_1 = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] / 2.0 + cabs(run_data->forcing_b[n]) * sin(carg(run_data->forcing_b[n]) - run_data->psi_n[n]));

		// Update the new magnetic field
		run_data->psi_n[n] = int_fac_b * run_data->psi_n[n] + dt * RK4_B1 * int_fac_b * RK_data->RK1_b[n] + dt * RK4_B2 * int_fac_b_1 * RK_data->RK2_b[n] + dt * RK4_B3 * int_fac_b_1 * RK_data->RK3_b[n] + dt * RK4_B4 * RK_data->RK4_b[n];
		#endif
		#else
		// Get the integrating factors
		int_fac_u   = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] + run_data->forcing_u[n]);
		int_fac_u_1 = cexp(-sys_vars->NU * dt * run_data->k[n] * run_data->k[n] / 2.0 + run_data->forcing_u[n]);

		// Update the new velocity field
		run_data->u[n] = int_fac_u * run_data->u[n] + dt * RK4_B1 * int_fac_u * RK_data->RK1_u[n] + dt * RK4_B2 * int_fac_u_1 * RK_data->RK2_u[n] + dt * RK4_B3 * int_fac_u_1 * RK_data->RK3_u[n] + dt * RK4_B4 * RK_data->RK4_u[n];
		#if defined(__MAGNETO)
		// Get the integrating factors
		int_fac_b   = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] + run_data->forcing_b[n]);
		int_fac_b_1 = cexp(-sys_vars->ETA * dt * run_data->k[n] * run_data->k[n] / 2.0 + run_data->forcing_b[n]);

		// Update the new magnetic field
		run_data->b[n] = int_fac_b * run_data->b[n] + dt * RK4_B1 * int_fac_b * RK_data->RK1_b[n] + dt * RK4_B2 * int_fac_b_1 * RK_data->RK2_b[n] + dt * RK4_B3 * int_fac_b_1 * RK_data->RK3_b[n] + dt * RK4_B4 * RK_data->RK4_b[n];
		#endif
		#endif

		#if defined(PHASE_ONLY)
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
	#if defined(PHASE_ONLY)
	double po_norm_fac_u;
	#if defined(__MAGNETO)
	double po_norm_fac_b;
	#endif
	#endif

	///----------------------- Forcing
	// Compute the forcing for the current iteration
	ComputeForicing(N);

	/////////////////////
	/// RK STAGES
	/////////////////////
	// ----------------------- Stage 1
	#if defined(PHASE_ONLY_DIRECT)
	NonlinearTerm(run_data->phi_n, run_data->psi_n, RK_data->RK1_u, RK_data->RK1_b, N);
	#else
	NonlinearTerm(run_data->u, run_data->b, RK_data->RK1_u, RK_data->RK1_b, N);
	#endif
	for (int i = 0; i < N; ++i) {
		// Get proper index
		n = i + 2;

		#if defined(PHASE_ONLY_DIRECT)
		// Add forcing & Update temorary velocity term
		if (run_data->a_n[n] != 0) {
			RK_data->RK1_u[n] += cabs(run_data->forcing_u[n]) * sin(carg(run_data->forcing_u[n]) - run_data->phi_n[n]) / run_data->a_n[n];
		}
		RK_data->RK_u_tmp[n] = run_data->phi_n[n] + dt * RK4_A21 * RK_data->RK1_u[n];
		#if defined(__MAGNETO)
		// Add forcing & Update temporary magnetic term
		if (run_data->b_n[n] != 0) {
			RK_data->RK1_b[n] += cabs(run_data->forcing_b[n]) * sin(carg(run_data->forcing_b[n]) - run_data->psi_n[n]) / run_data->b_n[n];
		}
		RK_data->RK_b_tmp[n] = run_data->psi_n[n] + dt * RK4_A21 * RK_data->RK1_b[n];		
		#endif
		#else
		// Add dissipative and forcing terms & Update temorary velocity term
		RK_data->RK1_u[n] += - sys_vars->NU * run_data->k[n] * run_data->k[n] * run_data->u[n] + run_data->forcing_u[n];
		RK_data->RK_u_tmp[n] = run_data->u[n] + dt * RK4_A21 * RK_data->RK1_u[n];
		#if defined(__MAGNETO)
		// Add dissipative and forcing terms & Update temporary magnetic term
		RK_data->RK1_b[n] += - sys_vars->ETA * run_data->k[n] * run_data->k[n] * run_data->b[n] + run_data->forcing_b[n];
		RK_data->RK_b_tmp[n] = run_data->b[n] + dt * RK4_A21 * RK_data->RK1_b[n];
		#endif
		#endif
	}

	// ----------------------- Stage 2
	NonlinearTerm(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK2_u, RK_data->RK2_b, N);
	for (int i = 0; i < N; ++i) {
		// Get proper index
		n = i + 2;

		#if defined(PHASE_ONLY_DIRECT)
		// Add forcing & Update temorary velocity term
		if (run_data->a_n[n] != 0) {
			RK_data->RK2_u[n] += cabs(run_data->forcing_u[n]) * sin(carg(run_data->forcing_u[n]) - run_data->phi_n[n]) / run_data->a_n[n];
		}
		RK_data->RK_u_tmp[n] = run_data->phi_n[n] + dt * RK4_A32 * RK_data->RK2_u[n];
		#if defined(__MAGNETO)
		// Add forcing & Update temporary magnetic term
		if (run_data->b_n[n] != 0) {
			RK_data->RK2_b[n] += cabs(run_data->forcing_b[n]) * sin(carg(run_data->forcing_b[n]) - run_data->psi_n[n]) / run_data->b_n[n];
		}
		RK_data->RK_b_tmp[n] = run_data->psi_n[n] + dt * RK4_A32 * RK_data->RK2_b[n];		
		#endif
		#else
		// Add dissipative and forcing terms & Update temorary velocity term
		RK_data->RK2_u[n] += - sys_vars->NU * run_data->k[n] * run_data->k[n] * run_data->u[n] + run_data->forcing_u[n];
		RK_data->RK_u_tmp[n] = run_data->u[n] + dt * RK4_A32 * RK_data->RK2_u[n];
		#if defined(__MAGNETO)
		// Add dissipative and forcing terms & Update temporary magnetic term
		RK_data->RK2_b[n] += - sys_vars->ETA * run_data->k[n] * run_data->k[n] * run_data->b[n] + run_data->forcing_b[n];
		RK_data->RK_b_tmp[n] = run_data->b[n] + dt * RK4_A32 * RK_data->RK2_b[n];
		#endif
		#endif
	}

	// ----------------------- Stage 3
	NonlinearTerm(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK3_u, RK_data->RK3_b, N);
	for (int i = 0; i < N; ++i) {
		// Get proper index
		n = i + 2;

		#if defined(PHASE_ONLY_DIRECT)
		// Add forcing & Update temorary velocity term
		if (run_data->a_n[n] != 0) {
			RK_data->RK3_u[n] += cabs(run_data->forcing_u[n]) * sin(carg(run_data->forcing_u[n]) - run_data->phi_n[n]) / run_data->a_n[n];
		}
		RK_data->RK_u_tmp[n] = run_data->phi_n[n] + dt * RK4_A43 * RK_data->RK3_u[n];
		#if defined(__MAGNETO)
		// Add forcing & Update temporary magnetic term
		if (run_data->b_n[n] != 0) {
			RK_data->RK3_b[n] += cabs(run_data->forcing_b[n]) * sin(carg(run_data->forcing_b[n]) - run_data->psi_n[n]) / run_data->b_n[n];
		}
		RK_data->RK_b_tmp[n] = run_data->psi_n[n] + dt * RK4_A43 * RK_data->RK3_b[n];		
		#endif
		#else
		// Add dissipative and forcing terms & Update temorary velocity term
		RK_data->RK3_u[n] += - sys_vars->NU * run_data->k[n] * run_data->k[n] * run_data->u[n] + run_data->forcing_u[n];
		RK_data->RK_u_tmp[n] = run_data->u[n] + dt * RK4_A43 * RK_data->RK3_u[n];
		#if defined(__MAGNETO)
		// Add dissipative and forcing terms & Update temporary magnetic term
		RK_data->RK3_b[n] += - sys_vars->ETA * run_data->k[n] * run_data->k[n] * run_data->b[n] + run_data->forcing_b[n];
		RK_data->RK_b_tmp[n] = run_data->b[n] + dt * RK4_A43 * RK_data->RK3_b[n];
		#endif
		#endif
	}

	// ----------------------- Stage 4
	NonlinearTerm(RK_data->RK_u_tmp, RK_data->RK_b_tmp, RK_data->RK4_u, RK_data->RK4_b, N);
	// Add the forcing and the linear term
	for (int i = 0; i < N; ++i) {
		// Get proper index
		n = i + 2;

		#if defined(PHASE_ONLY_DIRECT)
		if (run_data->a_n[n] != 0) {
			RK_data->RK4_u[n] += cabs(run_data->forcing_u[n]) * sin(carg(run_data->forcing_u[n]) - run_data->phi_n[n]) / run_data->a_n[n];
		}
		#if defined(__MAGNETO)
		if (run_data->b_n[n] != 0) {
			RK_data->RK4_b[n] += cabs(run_data->forcing_b[n]) * sin(carg(run_data->forcing_b[n]) - run_data->psi_n[n]) / run_data->b_n[n];
		}
		#endif
		#else
		RK_data->RK4_u[n] += - sys_vars->NU * run_data->k[n] * run_data->k[n] * run_data->u[n] + run_data->forcing_u[n];
		#if defined(__MAGNETO)
		RK_data->RK4_b[n] += - sys_vars->ETA * run_data->k[n] * run_data->k[n] * run_data->b[n] + run_data->forcing_b[n];
		#endif
		#endif
	}

	/////////////////////
	/// UPDATE STEP
	/////////////////////
	for (int i = 0; i < N; ++i) {
		// Get tmp index
		n = i + 2;

		#if defined(PHASE_ONLY)
		// Pre-record the amplitudes for resetting after update step
		po_norm_fac_u = cabs(run_data->u[n]);
		#if defined(__MAGNETO)
		po_norm_fac_b = cabs(run_data->b[n]);
		#endif
		#endif

		///-------------------- Update Step
		#if defined(PHASE_ONLY_DIRECT)
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

		#if defined(PHASE_ONLY)
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
		memcpy(&(RK_data->AB_tmp_nonlin_u[iters - 1][0]), RK_data->RK1_u, sizeof(fftw_complex) * (N + 4));
		#if defined(__MAGNETO)
		memcpy(&(RK_data->AB_tmp_nonlin_b[iters - 1][0]), RK_data->RK1_b, sizeof(fftw_complex) * (N + 4));		
		#endif
	}
	else {
		// -----------------------------------
		// Compute Forcing
		// -----------------------------------
		// Compute the forcing for the current iteration
		ComputeForicing(N);

		// -----------------------------------
		// Compute RHS
		// -----------------------------------
		// Get the nonlinear term for the current step
		NonlinearTerm(run_data->u, run_data->b, RK_data->AB_tmp_u, RK_data->AB_tmp_b, N);
		// Add the forcing 
		for (int i = 0; i < N + 4; ++i) {
			RK_data->AB_tmp_u[i] += run_data->forcing_u[i];
			#if defined(__MAGNETO)
			RK_data->AB_tmp_b[i] += run_data->forcing_b[i];
			#endif
		}

		/////////////////////
		/// AB Update Step
		/////////////////////
		for (int i = 0; i < N; ++i) {
			// Get tmp index
			n = i + 2;

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
		}

		// -----------------------------------
		// Update Previous Nonlinear Terms
		// -----------------------------------
		// Update the previous Nonlinear term arrays for next iteration
		memcpy(&(RK_data->AB_tmp_nonlin_u[0][0]), &(RK_data->AB_tmp_nonlin_u[1][0]), sizeof(fftw_complex) * (N + 4));
		memcpy(&(RK_data->AB_tmp_nonlin_u[1][0]), &(RK_data->AB_tmp_nonlin_u[2][0]), sizeof(fftw_complex) * (N + 4));
		memcpy(&(RK_data->AB_tmp_nonlin_u[2][0]), RK_data->AB_tmp_u, sizeof(fftw_complex) * (N + 4));
		#if defined(__MAGNETO)
		memcpy(&(RK_data->AB_tmp_nonlin_b[0][0]), &(RK_data->AB_tmp_nonlin_b[1][0]), sizeof(fftw_complex) * (N + 4));
		memcpy(&(RK_data->AB_tmp_nonlin_b[1][0]), &(RK_data->AB_tmp_nonlin_b[2][0]), sizeof(fftw_complex) * (N + 4));
		memcpy(&(RK_data->AB_tmp_nonlin_b[2][0]), RK_data->AB_tmp_b, sizeof(fftw_complex) * (N + 4));
		#endif
	}
}
#endif
#if defined(PHASE_ONLY_DIRECT) && !defined(PHASE_ONLY)
/**
 * Function that performs the evluation of the nonlinear term
 * @param u        array containing the input velocity phases
 * @param b        array containing the input magnetic phases
 * @param u_nonlin array to hold the result of computing the nonlinear term for velocity field
 * @param b_nonlin array to hold the result of computing the nonlinear term for the magnetic field
 * @param N        int defining the number of shells
 */
void NonlinearTerm(double* u, double* b, double* u_nonlin, double* b_nonlin, const long int N) {

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
	double u_tmp_1, u_tmp_2, u_tmp_3;
	#if defined(__MAGNETO)
	double b_tmp_1, b_tmp_2, b_tmp_3;
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
		// Compute the velocity phase temporary terms in the nonlinear term
		u_tmp_1 = run_data->a_n[n + 1] * run_data->a_n[n + 2] * cos(u[n] + u[n + 1] + u[n + 2]);
		u_tmp_2 = run_data->a_n[n - 1] * run_data->a_n[n + 1] * cos(u[n - 1] + u[n] + u[n + 1]);
		u_tmp_3 = run_data->a_n[n - 2] * run_data->a_n[n - 1] * cos(u[n - 2] + u[n - 1] + u[n]);
		#if defined(__MAGNETO)
		// Update velocity temporary terms with magnetic field
		u_tmp_1 -= run_data->b_n[n + 1] * run_data->b_n[n + 2] * cos(u[n] + b[n + 1] + b[n + 2]);
		u_tmp_2 -= run_data->b_n[n - 1] * run_data->b_n[n + 1] * cos(b[n - 1] + u[n] + b[n + 1]);
		u_tmp_3 -= run_data->b_n[n - 2] * run_data->b_n[n - 1] * cos(b[n - 2] + b[n - 1] + u[n]);

		// Compute the magnetic temporary terms in the nonlinear term
		b_tmp_1 = run_data->a_n[n + 1] * run_data->b_n[n + 2] * cos(b[n] + u[n + 1] + b[n + 2]) - run_data->b_n[n + 1] * run_data->a_n[n + 2] * cos(b[n] + b[n + 1] + u[n + 2]);
		b_tmp_2 = run_data->a_n[n - 1] * run_data->b_n[n + 1] * cos(u[n - 1] + b[n] + b[n + 1]) - run_data->b_n[n - 1] * run_data->a_n[n + 1] * cos(b[n - 1] + b[n] + u[n + 1]);
		b_tmp_3 = run_data->a_n[n - 2] * run_data->b_n[n - 1] * cos(u[n - 2] + b[n - 1] + b[n]) - run_data->b_n[n - 2] * run_data->a_n[n - 1] * cos(b[n - 2] + u[n - 1] + b[n]);
		#endif
		
		// -----------------------------------
		// Compute Nonlinear Terms
		// -----------------------------------
		// Compute the nonlinear term for the velocity field
		u_nonlin[n] = (run_data->k[n] / run_data->a_n[n]) * (u_tmp_1 - interact_coeff_u_1 * u_tmp_2 - interact_coeff_u_2 *  u_tmp_3);
		#if defined(__MAGNETO)
		// Compute the nonlinear term for the magnetic field
		b_nonlin[n] = (run_data->k[n] / run_data->b_n[n]) * (interact_coeff_b_1 * b_tmp_1 + interact_coeff_b_2 * b_tmp_2 + interact_coeff_b_3 * b_tmp_3); 
		#endif
	}
}
#else
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
	const double lambda_pow         = sys_vars->Lambda * sys_vars->Lambda;
	const double interact_coeff_u_1 = sys_vars->EPS / sys_vars->Lambda;
	const double interact_coeff_u_2 = (1.0 - sys_vars->EPS) / lambda_pow;
	#if defined(__MAGNETO)
	const double interact_coeff_b_1 = 1.0 - sys_vars->EPS - sys_vars->EPS_M;
	const double interact_coeff_b_2 = sys_vars->EPS_M / sys_vars->Lambda;
	const double interact_coeff_b_3 = (1.0 - sys_vars->EPS_M) / lambda_pow;
	#endif
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
		u_nonlin[n] = I * run_data->k[n] * conj(u_tmp_1 - interact_coeff_u_1 * u_tmp_2 - interact_coeff_u_2 *  u_tmp_3);
		#if defined(__MAGNETO)
		// Compute the nonlinear term for the magnetic field
		b_nonlin[n] = I * run_data->k[n] * conj(interact_coeff_b_1 * b_tmp_1 + interact_coeff_b_2 * b_tmp_2 + interact_coeff_b_3 * b_tmp_3); 
		#endif
	}
}
#endif
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


	for (int i = 0; i < N + 4; ++i) {

		// Initialize the edges shells
		if (i < 2 || i > N + 1) {
			run_data->u[i] = 0.0 + 0.0 * I;

			#if defined(PHASE_ONLY)
			// Record the phases and amplitudes
			run_data->a_n[i]   = 0.0;
			run_data->phi_n[i] = 0.0;
			#endif
			#if defined(PHASE_ONLY_DIRECT)
			run_data->a_n[i]   = 0.0;
			run_data->phi_n[i] = 0.0;
			#endif
			#if defined(__MAGNETO)
			// Initialize the magnetic field
			run_data->b[i] = 0.0 + 0.0 * I;
			#if defined(PHASE_ONLY)
			// Record the phases and amplitudes
			run_data->b_n[i]   = 0.0;
			run_data->psi_n[i] = 0.0;
			#endif
			#if defined(PHASE_ONLY_DIRECT)
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
				run_data->u[i] = 1.0 / pow(run_data->k[i], sys_vars->ALPHA) * cexp(I * pow(i - 1, 2.0));
				#if defined(PHASE_ONLY)
				// Record the phases and amplitudes
				run_data->a_n[i]   = cabs(run_data->u[i]);
				run_data->phi_n[i] = carg(run_data->u[i]);
				#endif
				#if defined(PHASE_ONLY_DIRECT)
				run_data->a_n[i]   = 1.0 / pow(run_data->k[i], sys_vars->ALPHA);
				run_data->phi_n[i] = pow(i - 1, 2.0);
				#endif

				#if defined(__MAGNETO)
				// Initialize the magnetic field
				run_data->b[i] = 1.0 / pow(run_data->k[i], sys_vars->BETA) * cexp(I * pow(i - 1, 4.0)) * 1e-2;
				#if defined(PHASE_ONLY)
				// Record the phases and amplitudes
				run_data->b_n[i]   = cabs(run_data->b[i]);
				run_data->psi_n[i] = carg(run_data->b[i]);
				#endif
				#if defined(PHASE_ONLY_DIRECT)
				run_data->b_n[i]   = 1.0 / pow(run_data->k[i], sys_vars->BETA) * 1e-2;
				run_data->psi_n[i] = pow(i - 1, 4.0);
				#endif
				#endif
			}
			else if(!(strcmp(sys_vars->u0, "RANDOM"))) {
				// ------------------------------------------------
				// Default - Random Initial Conditions
				// ------------------------------------------------	
				// Get random uniform number
				r1 = (double)rand() / (double)RAND_MAX;
				#if defined(__MAGNETO)
				r2 = (double)rand() / (double)RAND_MAX;
				#endif

				// Initialize the velocity field
				run_data->u[i] = 1.0 / pow(run_data->k[i], sys_vars->ALPHA) * cexp(I * r1 * 2.0 * M_PI);

				#if defined(PHASE_ONLY)
				// Record the phases and amplitudes
				run_data->a_n[i]   = cabs(run_data->u[i]);
				run_data->phi_n[i] = carg(run_data->u[i]);
				#endif
				#if defined(PHASE_ONLY_DIRECT)
				run_data->a_n[i]   = 1.0 / pow(run_data->k[i], sys_vars->ALPHA);
				run_data->phi_n[i] = r1 * 2.0 * M_PI;
				#endif

				#if defined(__MAGNETO)
				// Initialize the magnetic field
				run_data->b[i] = 1.0 / pow(run_data->k[i], sys_vars->BETA) * cexp(I * r2 * 2.0 * M_PI) * 1e-2;
				#if defined(PHASE_ONLY)
				// Record the phases and amplitudes
				run_data->b_n[i]   = cabs(run_data->b[i]);
				run_data->psi_n[i] = carg(run_data->b[i]);
				#endif
				#if defined(PHASE_ONLY_DIRECT)
				run_data->b_n[i]   = 1.0 / pow(run_data->k[i], sys_vars->BETA) * 1e-2;
				run_data->psi_n[i] = r2 * 2.0 * M_PI;
				#endif
				#endif
			}
			else if(!(strcmp(sys_vars->u0, "ZERO"))) {
				// ------------------------------------------------
				// Zero Initial Conditions
				// ------------------------------------------------	
				// Initialize the velocity field
				run_data->u[i] = 1.0 / pow(run_data->k[i], sys_vars->ALPHA);

				#if defined(PHASE_ONLY)
				// Record the phases and amplitudes
				run_data->a_n[i]   = cabs(run_data->u[i]);
				run_data->phi_n[i] = carg(run_data->u[i]);
				#endif
				#if defined(PHASE_ONLY_DIRECT)
				run_data->a_n[i]   = 1.0 / pow(run_data->k[i], sys_vars->ALPHA);
				run_data->phi_n[i] = 0.0;
				#endif

				#if defined(__MAGNETO)
				// Initialize the magnetic field
				run_data->b[i] = 1.0 / pow(run_data->k[i], sys_vars->BETA) * 1e-2;
				#if defined(PHASE_ONLY)
				// Record the phases and amplitudes
				run_data->b_n[i]   = cabs(run_data->b[i]);
				run_data->psi_n[i] = carg(run_data->b[i]);
				#endif
				#if defined(PHASE_ONLY_DIRECT)
				run_data->b_n[i]   = 1.0 / pow(run_data->k[i], sys_vars->BETA) * 1e-2;
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
				#if defined(__MAGNETO)
				r2 = (double)rand() / (double)RAND_MAX;
				r4 = (double)rand() / (double)RAND_MAX;
				#endif

				// Initialize the velocity field
				run_data->u[i] = r1 + r3 * I;

				#if defined(PHASE_ONLY)
				// Record the phases and amplitudes
				run_data->a_n[i]   = cabs(run_data->u[i]);
				run_data->phi_n[i] = carg(run_data->u[i]);
				#endif
				#if defined(PHASE_ONLY_DIRECT)
				run_data->a_n[i]   = r1;
				run_data->phi_n[i] = r3;
				#endif

				#if defined(__MAGNETO)
				// Initialize the magnetic field
				run_data->b[i] = (r2 + r4) * 1e-2;
				#if defined(PHASE_ONLY)
				// Record the phases and amplitudes
				run_data->b_n[i]   = cabs(run_data->b[i]);
				run_data->psi_n[i] = carg(run_data->b[i]);
				#endif
				#if defined(PHASE_ONLY_DIRECT)
				run_data->b_n[i]   = r2;
				run_data->psi_n[i] = r4;
				#endif
				#endif
			}
		}
		// printf("u[%d]:\t%1.16lf\t%1.16lf i\tb[%d]:\t%1.16lf\t%1.16lf i\n", i - 1, creal(run_data->u[i]), cimag(run_data->u[i]),  i - 1, creal(run_data->b[i]), cimag(run_data->b[i]));		
		// printf("a_n[%d]:\t%1.16lf\tphi[%d]:\t%1.16lf\tb_n[%d]:\t%1.16lf\tpsi[%d]:\t%1.16lf\n", i, run_data->a_n[i], i, run_data->phi_n[i], i, run_data->b_n[i], i, run_data->psi_n[i]);
	}
	printf("\n");
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
			k[i] = sys_vars->k_0 * pow(sys_vars->Lambda, i - 2);
		}
		else {
			k[i] = 0.0;
		}
		// printf("k[%d]: %lf\n", i - 1, k[i]);
	}
	// printf("\n");

}
/**
 * Function to initialize the forcing 
 * @param N Number of shells
 */
void InitializeForicing(const long int N) {

	// ------------------------------------------------
	// Allocate Memory for Forcing Data
	// ------------------------------------------------
	run_data->forcing_u = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	run_data->forcing_b = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));

	// ------------------------------------------------
	// Initialize Forcing Data
	// ------------------------------------------------
	for (int i = 0; i < N + 4; ++i) {
		// Compute the forcing and intialize
		if(!(strcmp(sys_vars->forcing, "DELTA"))) {
			// ------------------------------------------------
			// Delta function on the zero mode
			// ------------------------------------------------
			// Record the forcing 
			if (i - 1 == sys_vars->force_k) {
				run_data->forcing_u[i] = sys_vars->force_scale_var * (1.0 + 1.0 * I) * my_delta(0.0, 0.0);
				#if defined(__MAGNETO)
				run_data->forcing_b[i] = sys_vars->force_scale_var * (1.0 + 1.0 * I) * my_delta(0.0, 0.0);
				#endif				
			}
			else {
				run_data->forcing_u[i] = 0.0 + 0.0 * I;
				#if defined(__MAGNETO)
				run_data->forcing_b[i] = 0.0 + 0.0 * I;
				#endif			
			}
			
		}
		else if(!(strcmp(sys_vars->forcing, "STOC"))) {
			// ------------------------------------------------
			// Stochastic Forcing
			// ------------------------------------------------
			
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
	}
	// printf("\n");
}
/**
 * Function compute the forcing for the current timestep
 * @param N Number of shells
 */
void ComputeForicing(const long int N) {

	// Initialize variables
	int n;

	// ------------------------------------------------
	// Initialize Forcing Data
	// ------------------------------------------------
	for (int i = 0; i < N + 4; ++i) {
		// Get temporary index
		n = i;

		// Compute the forcing and intialize
		 if(!(strcmp(sys_vars->forcing, "STOC"))) {
			// ------------------------------------------------
			// Stochastic Forcing
			// ------------------------------------------------
			
		}
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
	sys_vars->min_dt = MIN_STEP_SIZE;
	sys_vars->max_dt = 10;

	// -------------------------------
	// Integration Counters
	// -------------------------------
	// Number of time steps and saving steps
	sys_vars->num_t_steps = (int)round( ((*T) - (*t0)) / (*dt) );
	// printf("\n\n%1.16lf,%1.16lf,%1.16lf--%1.16lf\t%1.16lf\t%1.16lf\n\n", (*T), (*t0), (*dt), ((*T) - (*t0)) / (*dt), round(((*T) - (*t0)) / (*dt)));
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
	printf("Iter: %d/%ld\tt: %1.6lf/%1.3lf\tdt: %1.6g\tKE: %6.6g\tHEL_U: %6.6g\tHEL_B: %6.6g\tX-HEL: %6.6g\n", iters, sys_vars->num_t_steps, t, T, dt, run_data->tot_energy[save_data_indx], run_data->tot_hel_u[save_data_indx], run_data->tot_hel_b[save_data_indx], run_data->tot_cross_hel[save_data_indx]);
	#else
	printf("Iter: %d/%ld\tt: %1.6lf/%1.3lf\tdt: %1.6g\tKE: %6.6g\tHEL: %6.6g\t\n", iters, sys_vars->num_t_steps, t, T, dt, run_data->tot_energy[save_data_indx], run_data->tot_hel_u[save_data_indx]);
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
	run_data->k = (double* )fftw_malloc(sizeof(double) * (N + 4));  // k
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
	#if defined(PHASE_ONLY)	|| defined(PHASE_ONLY_DIRECT)
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
	#if defined(PHASE_ONLY)	|| defined(PHASE_ONLY_DIRECT)
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
	#if defined(PHASE_ONLY_DIRECT)
	RK_data->RK1_u       = (double* )fftw_malloc(sizeof(double) * (N + 4));
	RK_data->RK2_u       = (double* )fftw_malloc(sizeof(double) * (N + 4));
	RK_data->RK3_u       = (double* )fftw_malloc(sizeof(double) * (N + 4));
	RK_data->RK4_u       = (double* )fftw_malloc(sizeof(double) * (N + 4));
	RK_data->RK_u_tmp    = (double* )fftw_malloc(sizeof(double) * (N + 4));
	#if defined(AB4CN)
	RK_data->AB_tmp_u	 = (double* )fftw_malloc(sizeof(double) * (N + 4));
	for (int i = 0; i < 3; ++i) {
		RK_data->AB_tmp_nonlin_u[i] = (double* )fftw_malloc(sizeof(double) * (N + 4));
	}
	#endif
	#if defined(__MAGNETO)
	RK_data->RK1_b       = (double* )fftw_malloc(sizeof(double) * (N + 4));
	RK_data->RK2_b       = (double* )fftw_malloc(sizeof(double) * (N + 4));
	RK_data->RK3_b       = (double* )fftw_malloc(sizeof(double) * (N + 4));
	RK_data->RK4_b       = (double* )fftw_malloc(sizeof(double) * (N + 4));
	RK_data->RK_b_tmp    = (double* )fftw_malloc(sizeof(double) * (N + 4));
	#if defined(AB4CN)
	RK_data->AB_tmp_b	 = (double* )fftw_malloc(sizeof(double) * (N + 4));
	for (int i = 0; i < 3; ++i) {
		RK_data->AB_tmp_nonlin_b[i] = (double* )fftw_malloc(sizeof(double) * (N + 4));
	}
	#endif
	#endif
	#else
	RK_data->RK1_u       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	RK_data->RK2_u       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	RK_data->RK3_u       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	RK_data->RK4_u       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	RK_data->RK_u_tmp    = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	#if defined(AB4CN)
	RK_data->AB_tmp_u	 = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	for (int i = 0; i < 3; ++i) {
		RK_data->AB_tmp_nonlin_u[i] = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	}
	#endif
	#if defined(__MAGNETO)
	RK_data->RK1_b       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	RK_data->RK2_b       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	RK_data->RK3_b       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	RK_data->RK4_b       = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	RK_data->RK_b_tmp    = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	#if defined(AB4CN)
	RK_data->AB_tmp_b	 = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	for (int i = 0; i < 3; ++i) {
		RK_data->AB_tmp_nonlin_b[i] = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * (N + 4));
	}
	#endif
	#endif
	#endif
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
		#if defined(PHASE_ONLY)
		run_data->a_n[i]   = 0.0;
		run_data->phi_n[i] = 0.0;
		#endif
		#if defined(__MAGNETO)
		run_data->b[i]       = 0.0 + 0.0 * I;
		RK_data->RK1_b[i]    = 0.0 + 0.0 * I;
		RK_data->RK2_b[i]    = 0.0 + 0.0 * I;
		RK_data->RK3_b[i]    = 0.0 + 0.0 * I;
		RK_data->RK4_b[i]    = 0.0 + 0.0 * I;
		RK_data->RK_b_tmp[i] = 0.0 + 0.0 * I;
		#if defined(AB4CN)
		RK_data->AB_tmp_b[i] = 0.0 + 0.0 * I;
		for (int j = 0; j < 3; ++j) {
			RK_data->AB_tmp_nonlin_b[j][i] = 0.0 + 0.0 * I;
		}
		#endif
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
	#if defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT)
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
	fftw_free(run_data->tot_hel_u);
	#if defined(__MAGNETO)
	fftw_free(run_data->tot_hel_b);
	fftw_free(run_data->tot_cross_hel);
	#endif
	#endif
	#if defined(__TIME)
	fftw_free(run_data->time);
	#endif
	#if defined(__ENRG_FLUX)
	fftw_free(run_data->energy_flux);
	fftw_free(run_data->energy_diss_u);
	fftw_free(run_data->energy_input_u);
	#if defined(__MAGNETO)
	fftw_free(run_data->energy_diss_b);
	fftw_free(run_data->energy_input_b);
	#endif
	#endif
	fftw_free(run_data->forcing_u);
	fftw_free(run_data->forcing_b);

	// Free stats objects
	#if defined(STATS)
	for (int i = 0; i < NUM_POW - 2; ++i) {
		#if defined(__STR_FUNC_VEL)
		fftw_free(stats_data->vel_str_func[i]);
		#endif
		#if defined(__STR_FUNC_VEL_FLUX)
		fftw_free(stats_data->vel_flux_str_func[0][i]);
		fftw_free(stats_data->vel_flux_str_func[1][i]);
		#endif
		#if defined(__MAGNETO)
		#if defined(__STR_FUNC_MAG)
		fftw_free(stats_data->mag_str_func[i]);
		#endif
		#if defined(__STR_FUNC_MAG_FLUX)
		fftw_free(stats_data->mag_fluxstr_func[0][i]);
		fftw_free(stats_data->mag_fluxstr_func[1][i]);
		#endif
		#endif
	}
	for (int i = 0; i < sys_vars->N; ++i) {
		gsl_rstat_free(stats_data->vel_moments[i]);
		#if defined(__MAGNETO)
		gsl_rstat_free(stats_data->mag_moments[i]);
		#endif
	}
	#endif


	// Free integration variables
	fftw_free(RK_data->RK1_u);
	fftw_free(RK_data->RK2_u);
	fftw_free(RK_data->RK3_u);
	fftw_free(RK_data->RK4_u);
	fftw_free(RK_data->RK_u_tmp);
	#if defined(AB4CN)
	fftw_free(RK_data->AB_tmp_u);
	for (int i = 0; i < 3; ++i)	{
		fftw_free(RK_data->AB_tmp_nonlin_u[i]);
	}
	#endif
	#if defined(__MAGNETO)
	fftw_free(RK_data->RK1_b);
	fftw_free(RK_data->RK2_b);
	fftw_free(RK_data->RK3_b);
	fftw_free(RK_data->RK4_b);
	fftw_free(RK_data->RK_b_tmp);
	#if defined(AB4CN)
	fftw_free(RK_data->AB_tmp_b);
	for (int i = 0; i < 3; ++i)	{
		fftw_free(RK_data->AB_tmp_nonlin_b[i]);
	}
	#endif
	#endif
}
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------