/**
* @file phase_sync.c  
* @author Enda Carroll
* @date Sept 2022
* @brief File containing the phase sync functions for the solver
*/
// ---------------------------------------------------------------------
//  Standard Libraries and Headers
// ---------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h> 
#include <math.h>
#include <complex.h>

// ---------------------------------------------------------------------
//  User Libraries and Headers
// ---------------------------------------------------------------------
#include "data_types.h"
#include "utils.h"
// ---------------------------------------------------------------------
//  Function Definitions
// ---------------------------------------------------------------------
void ComputePhaseSyncData(const long int iter) {

	// Initialize variables
	int n;
	int N              = sys_vars->N;
	int num_triads     = phase_sync->num_triads;
	int num_phase_diff = phase_sync->num_phase_diff;
	int gsl_status;
	double phase_u;
	#if defined(__MAGNETO)
	double phase_b_1, phase_b_2, phase_b_3;	
	#endif

	// Update the phse sync counter
	phase_sync->num_phase_sync_steps++;

	// -------------------------------
	// Compute Phases & Phase Order
	// -------------------------------
	for (int i = 0; i < N; ++i)	{
		n = i + 2;

		///--------------------------------- Compute the triad phases and triad order parameters
		if (i < num_triads) {
			// Get the triad phase
			#if defined(PHASE_ONLY)
			phase_u = fmod(run_data->phi_n[n] + run_data->phi_n[n + 1] + run_data->phi_n[n + 2] + 6.0 * M_PI, 2.0 * M_PI);
			#if defined(__MAGNETO)
			phase_b_1 = fmod(run_data->phi_n[n] + run_data->psi_n[n + 1] + run_data->psi_n[n + 2] + 6.0 * M_PI, 2.0 * M_PI);
			phase_b_2 = fmod(run_data->psi_n[n] + run_data->phi_n[n + 1] + run_data->psi_n[n + 2] + 6.0 * M_PI, 2.0 * M_PI);
			phase_b_3 = fmod(run_data->psi_n[n] + run_data->psi_n[n + 1] + run_data->phi_n[n + 2] + 6.0 * M_PI, 2.0 * M_PI);
			#endif
			#else 
			phase_u = fmod(carg(run_data->u[n]) + carg(run_data->u[n + 1]) + carg(run_data->u[n + 2]) + 6.0 * M_PI, 2.0 * M_PI);
			#if defined(__MAGNETO)
			phase_b_1 = fmod(carg(run_data->u[n]) + carg(run_data->b[n + 1]) + carg(run_data->b[n + 2]) + 6.0 * M_PI, 2.0 * M_PI);
			phase_b_2 = fmod(carg(run_data->b[n]) + carg(run_data->u[n + 1]) + carg(run_data->b[n + 2]) + 6.0 * M_PI, 2.0 * M_PI);
			phase_b_3 = fmod(carg(run_data->b[n]) + carg(run_data->b[n + 1]) + carg(run_data->u[n + 2]) + 6.0 * M_PI, 2.0 * M_PI);
			#endif
			#endif
			
			///--------------------- Record the triad phases
			phase_sync->triads_u[i] = phase_u;
			#if defined(__MAGNETO)
			phase_sync->triads_b[0 * phase_sync->num_triads + i] = phase_b_1;
			phase_sync->triads_b[1 * phase_sync->num_triads + i] = phase_b_2;
			phase_sync->triads_b[2 * phase_sync->num_triads + i] = phase_b_3;
			#endif

			///--------------------- Record the triad phase order parameters
			#if defined(__PHASE_SYNC)
			phase_sync->triad_u_order[i] += cexp(I * phase_u);
			#if defined(__MAGNETO)
			phase_sync->triad_b_order[0 * phase_sync->num_triads + i] += cexp(I * phase_b_1);
			phase_sync->triad_b_order[1 * phase_sync->num_triads + i] += cexp(I * phase_b_2);
			phase_sync->triad_b_order[2 * phase_sync->num_triads + i] += cexp(I * phase_b_3);
			#endif
			#endif

			///--------------------- Record Phase Sync Stats
			#if defined(__PHASE_SYNC_STATS)
			gsl_status = gsl_histogram_increment(phase_sync->triad_u_hist[i], phase_u);
			if (gsl_status != 0) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update bin count for ["CYAN"%s"RESET"] for Iter ["CYAN"%ld"RESET"] -- GSL Exit Status [Err:"CYAN" %d"RESET" - Val:"CYAN" %lf"RESET"]\n-->> Exiting!!!\n", "Velocity Triads Histogram", iter, gsl_status, phase_u);
				exit(1);
			}
			#if defined(__MAGNETO)
			gsl_status = gsl_histogram_increment(phase_sync->triad_b_hist[0 * phase_sync->num_triads + i], phase_b_1);
			if (gsl_status != 0) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update bin count for ["CYAN"%s"RESET"] for Iter ["CYAN"%ld"RESET"] -- GSL Exit Status [Err:"CYAN" %d"RESET" - Val:"CYAN" %lf"RESET"]\n-->> Exiting!!!\n", "Magnetic Triads Histogram 1", iter, gsl_status, phase_b_1);
				exit(1);
			}
			gsl_status = gsl_histogram_increment(phase_sync->triad_b_hist[1 * phase_sync->num_triads + i], phase_b_2);
			if (gsl_status != 0) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update bin count for ["CYAN"%s"RESET"] for Iter ["CYAN"%ld"RESET"] -- GSL Exit Status [Err:"CYAN" %d"RESET" - Val:"CYAN" %lf"RESET"]\n-->> Exiting!!!\n", "Magnetic Triads Histogram 2", iter, gsl_status, phase_b_2);
				exit(1);
			}
			gsl_status = gsl_histogram_increment(phase_sync->triad_b_hist[2 * phase_sync->num_triads + i], phase_b_3);
			if (gsl_status != 0) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update bin count for ["CYAN"%s"RESET"] for Iter ["CYAN"%ld"RESET"] -- GSL Exit Status [Err:"CYAN" %d"RESET" - Val:"CYAN" %lf"RESET"]\n-->> Exiting!!!\n", "Magnetic Triads Histogram 3", iter, gsl_status, phase_b_3);
				exit(1);
			}
			#endif
			#endif
		}

		///--------------------------------- Compute the phase differences and phase difference order parameters
		if (i < num_phase_diff) {
			// Get the phase differences
			#if defined(PHASE_ONLY)
			phase_u = fmod(run_data->phi_n[n] - run_data->phi_n[n + 3] + 4.0 * M_PI, 2.0 * M_PI);
			#if defined(__MAGNETO)
			phase_b_1 = fmod(run_data->psi_n[n] - run_data->psi_n[n + 3] + 4.0 * M_PI, 2.0 * M_PI);
			#endif
			#else 
			phase_u = fmod(carg(run_data->u[n]) - carg(run_data->u[n + 3]) + 4.0 * M_PI, 2.0 * M_PI);
			#if defined(__MAGNETO)
			phase_b_1 = fmod(carg(run_data->b[n]) - carg(run_data->b[n + 3]) + 4.0 * M_PI, 2.0 * M_PI);
			#endif
			#endif
			
			///--------------------- Record the phase differences
			phase_sync->phase_diff_u[i] = phase_u;
			#if defined(__MAGNETO)
			phase_sync->phase_diff_b[i] = phase_b_1;
			#endif

			///--------------------- Record the phase difference order parameters
			#if defined(__PHASE_SYNC)
			phase_sync->phase_diff_u_order[i] += cexp(I * phase_u);
			#if defined(__MAGNETO)
			phase_sync->phase_diff_b_order[i] += cexp(I * phase_b_1);
			#endif
			#endif

			///--------------------- Record Phase Sync Stats
			#if defined(__PHASE_SYNC_STATS)
			gsl_status = gsl_histogram_increment(phase_sync->phase_diff_u_hist[i], phase_u);
			if (gsl_status != 0) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update bin count for ["CYAN"%s"RESET"] for Iter ["CYAN"%ld"RESET"] -- GSL Exit Status [Err:"CYAN" %d"RESET" - Val:"CYAN" %lf"RESET"]\n-->> Exiting!!!\n", "Velocity Phase Difference Histogram", iter, gsl_status, phase_u);
				exit(1);
			}
			#if defined(__MAGNETO)
			gsl_status = gsl_histogram_increment(phase_sync->phase_diff_b_hist[i], phase_b_1);
			if (gsl_status != 0) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update bin count for ["CYAN"%s"RESET"] for Iter ["CYAN"%ld"RESET"] -- GSL Exit Status [Err:"CYAN" %d"RESET" - Val:"CYAN" %lf"RESET"]\n-->> Exiting!!!\n", "Magnetic Phase Difference Histogram", iter, gsl_status, phase_b_1);
				exit(1);
			}
			#endif
			#endif
		}
	}	
}
/**
 * Wrapper funciton to allocate the phase sync memory and phase sync stats
 */
void InitializePhaseSyncObjects(void) {

	// Initialize variables
	phase_sync->num_triads     = sys_vars->N - 2;
	phase_sync->num_phase_diff = sys_vars->N - 3;
	int num_triads             = phase_sync->num_triads;
	int num_phase_diff         = phase_sync->num_phase_diff;


	// -------------------------------
	// Allocate Triads Memory
	// -------------------------------
	///---------------- Allocate Velocity Triads
	phase_sync->triads_u = (double* )malloc(sizeof(double) * num_triads);
	if (phase_sync->triads_u == NULL) {
	    fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Velocity Triads");
	    exit(1);
	}
	
	///---------------- Allocate Magentic Triads
	#if defined(__MAGNETO)
	phase_sync->triads_b = (double* )malloc(sizeof(double) * num_triads * NUM_MAG_TRIAD_TYPES);
	if (phase_sync->triads_b == NULL) {
	    fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Triads");
	    exit(1);
	}
	#endif

	// -------------------------------
	// Allocate Phase Difference Memory
	// -------------------------------
	///---------------- Allocate Velocity Phase Difference
	phase_sync->phase_diff_u = (double* )malloc(sizeof(double) * num_phase_diff);
	if (phase_sync->phase_diff_u == NULL) {
	    fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Velocity Phase Difference");
	    exit(1);
	}
	
	///---------------- Allocate Magentic Phase Difference
	#if defined(__MAGNETO)
	phase_sync->phase_diff_b = (double* )malloc(sizeof(double) * num_phase_diff);
	if (phase_sync->phase_diff_b == NULL) {
	    fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Phase Difference");
	    exit(1);
	}
	#endif


	// -------------------------------
	// Allocate Kuramoto Order Params
	// -------------------------------
	#if defined(__PHASE_SYNC)
	///----------------- Velocity Triad Order Parameter
	phase_sync->triad_u_order = (complex double* )malloc(sizeof(complex double) * num_triads);
	if (phase_sync->triad_u_order == NULL) {
	    fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Velocity Triad Order Parameter");
	    exit(1);
	}

	///----------------- Magentic Triad Order Parameter
	#if defined(__MAGNETO)
	phase_sync->triad_b_order = (complex double* )malloc(sizeof(complex double) * num_triads * NUM_MAG_TRIAD_TYPES);
	if (phase_sync->triad_b_order == NULL) {
	    fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Triad Order Parameter");
	    exit(1);
	}
	#endif

	///----------------- Velocity Phase Difference Order Parameter
	phase_sync->phase_diff_u_order = (complex double* )malloc(sizeof(complex double) * num_phase_diff);
	if (phase_sync->phase_diff_u_order == NULL) {
	    fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Velocity Phase Difference Order Parameter");
	    exit(1);
	}

	///----------------- Magentic Phase Difference Order Parameter
	#if defined(__MAGNETO)
	phase_sync->phase_diff_b_order = (complex double* )malloc(sizeof(complex double) * num_phase_diff);
	if (phase_sync->phase_diff_b_order == NULL) {
	    fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Phase Difference Order Parameter");
	    exit(1);
	}
	#endif
	#endif


	// -------------------------------
	// Allocate Phase Sync Stats
	// -------------------------------
	#if defined(__PHASE_SYNC_STATS)
	int gsl_status; 

	///---------------- Allocate Velocity Histogram memory
	phase_sync->triad_u_hist = (gsl_histogram** )malloc(sizeof(gsl_histogram*) * num_triads);
	if (phase_sync->triad_u_hist == NULL) {
	    fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Velocity Triads Histograms");
	    exit(1);
	}
	phase_sync->phase_diff_u_hist = (gsl_histogram** )malloc(sizeof(gsl_histogram*) * num_phase_diff);
	if (phase_sync->phase_diff_u_hist == NULL) {
	    fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Velocity Phase Difference Histograms");
	    exit(1);
	}

	// Initialize Histograms
	for (int i = 0; i < num_triads; ++i) {
		// Allocate  number of bins
		phase_sync->triad_u_hist[i] = gsl_histogram_alloc(NUM_PHASE_SYNC_HIST_BINS);
		
		// Set the bin ranges
		gsl_status = gsl_histogram_set_ranges_uniform(phase_sync->triad_u_hist[i], -0.0 - 0.05, 2.0 * M_PI + 0.05);
		if (gsl_status != 0) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to set bin ranges for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Velocity Triads Histogram");
			exit(1);
		}

		if (i < num_phase_diff) {
			// Allocate  number of bins
			phase_sync->phase_diff_u_hist[i] = gsl_histogram_alloc(NUM_PHASE_SYNC_HIST_BINS);

			// Set the bin ranges
			gsl_status = gsl_histogram_set_ranges_uniform(phase_sync->phase_diff_u_hist[i], -0.0 - 0.05, 2.0 * M_PI + 0.05);
			if (gsl_status != 0) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to set bin ranges for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Velocity Phase Difference Histogram");
				exit(1);
			}
		}
	}

	///---------------- Allocate Magentic Field Phase Histogram memory
	#if defined(__MAGNETO)
	phase_sync->triad_b_hist = (gsl_histogram** )malloc(sizeof(gsl_histogram*) * num_triads * NUM_MAG_TRIAD_TYPES);
	if (phase_sync->triad_b_hist == NULL) {
	    fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Triads Histograms");
	    exit(1);
	}
	phase_sync->phase_diff_b_hist = (gsl_histogram** )malloc(sizeof(gsl_histogram*) * num_phase_diff);
	if (phase_sync->phase_diff_b_hist == NULL) {
	    fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Phase Difference Histograms");
	    exit(1);
	}

	// Initialize Histograms
	for (int i = 0; i < num_triads; ++i) {
		// Allocate  number of bins
		phase_sync->triad_b_hist[0 * phase_sync->num_triads + i] = gsl_histogram_alloc(NUM_PHASE_SYNC_HIST_BINS);
		phase_sync->triad_b_hist[1 * phase_sync->num_triads + i] = gsl_histogram_alloc(NUM_PHASE_SYNC_HIST_BINS);
		phase_sync->triad_b_hist[2 * phase_sync->num_triads + i] = gsl_histogram_alloc(NUM_PHASE_SYNC_HIST_BINS);
		
		// Set the bin ranges
		gsl_status = gsl_histogram_set_ranges_uniform(phase_sync->triad_b_hist[0 * NUM_MAG_TRIAD_TYPES + i], -0.0 - 0.05, 2.0 * M_PI + 0.05);
		if (gsl_status != 0) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to set bin ranges for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Triads Histogram");
			exit(1);
		}
		gsl_status = gsl_histogram_set_ranges_uniform(phase_sync->triad_b_hist[1 * NUM_MAG_TRIAD_TYPES + i], -0.0 - 0.05, 2.0 * M_PI + 0.05);
		if (gsl_status != 0) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to set bin ranges for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Triads Histogram");
			exit(1);
		}
		gsl_status = gsl_histogram_set_ranges_uniform(phase_sync->triad_b_hist[2 * NUM_MAG_TRIAD_TYPES + i], -0.0 - 0.05, 2.0 * M_PI + 0.05);
		if (gsl_status != 0) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to set bin ranges for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Triads Histogram");
			exit(1);
		}

		if (i < num_phase_diff) {
			// Allocate  number of bins
			phase_sync->phase_diff_b_hist[i] = gsl_histogram_alloc(NUM_PHASE_SYNC_HIST_BINS);

			// Set the bin ranges
			gsl_status = gsl_histogram_set_ranges_uniform(phase_sync->phase_diff_b_hist[i], -0.0 - 0.05, 2.0 * M_PI + 0.05);
			if (gsl_status != 0) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to set bin ranges for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Phase Difference Histogram");
				exit(1);
			}
		}
	}
	#endif
	#endif
}
/**
 * Wrapper function to free phase sync memory and clean up phase sync stats objects
 */
void FreePhaseSyncObjects(void) {

	// -------------------------------
	// Free Phase Memory
	// -------------------------------
	free(phase_sync->triads_u);
	free(phase_sync->triads_b);
	free(phase_sync->phase_diff_u);
	free(phase_sync->phase_diff_b);

	// -------------------------------
	// Free Order Parameter Memory
	// -------------------------------
	#if defined(__PHASE_SYNC)
	free(phase_sync->triad_u_order);
	free(phase_sync->triad_b_order);
	free(phase_sync->phase_diff_u_order);
	free(phase_sync->phase_diff_b_order);
	#endif

	// -------------------------------
	// Free Order Stats objects
	// -------------------------------
	#if defined(__PHASE_SYNC_STATS)
	for (int i = 0; i < phase_sync->num_triads; ++i) {
		gsl_histogram_free(phase_sync->triad_u_hist[i]);
		gsl_histogram_free(phase_sync->triad_b_hist[0 * phase_sync->num_triads + i]);
		gsl_histogram_free(phase_sync->triad_b_hist[1 * phase_sync->num_triads + i]);
		gsl_histogram_free(phase_sync->triad_b_hist[2 * phase_sync->num_triads + i]);
		if (i < phase_sync->num_phase_diff) {
			gsl_histogram_free(phase_sync->phase_diff_u_hist[i]);
			gsl_histogram_free(phase_sync->phase_diff_b_hist[i]);
		}
	}
	#endif
}
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------