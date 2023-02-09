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
	double phase_u_1, phase_u_2, phase_u_3;
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	double phase_b_type1, phase_b_type2, phase_b_type3;	
	double phase_b_1, phase_b_2, phase_b_3;	
	#endif

	// Update the phse sync counter
	phase_sync->num_phase_sync_steps++;

	// -------------------------------
	// Get the Phases
	// -------------------------------
	#if !defined(PHASE_ONLY)
	for (int i = 0; i < N; ++i)	{
		n = i + 2;

		#if defined(__ELSASSAR_MHD)
		run_data->u[n] = (run_data->z_plus[n] + run_data->z_minus[n]) / 2.0;
		run_data->b[n] = (run_data->z_plus[n] - run_data->z_minus[n]) / 2.0;
		#endif

		// Get the phases
		run_data->phi_n[n] = fmod(carg(run_data->u[n]) + 2.0 * M_PI, 2.0 * M_PI);	
		#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
		run_data->psi_n[n] = fmod(carg(run_data->b[n]) + 2.0 * M_PI, 2.0 * M_PI);
		#endif
	}
	#endif


	// -------------------------------
	// Compute Phases & Phase Order
	// -------------------------------
	for (int i = 0; i < N; ++i)	{
		n = i + 2;

		///--------------------------------- Compute the triad phases and triad order parameters
		if (i < num_triads) {

			// Get the triad phase
			phase_u_1 = fmod(run_data->phi_n[n] + 2.0 * M_PI, 2.0 * M_PI);
			phase_u_2 = fmod(run_data->phi_n[n + 1] + 2.0 * M_PI, 2.0 * M_PI);
			phase_u_3 = fmod(run_data->phi_n[n + 2] + 2.0 * M_PI, 2.0 * M_PI);
			phase_u   = fmod(phase_u_1 + phase_u_2 + phase_u_3 + 2.0 * M_PI, 2.0 * M_PI);
			#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
			phase_b_type1 = fmod(run_data->phi_n[n] + run_data->psi_n[n + 1] + run_data->psi_n[n + 2] + 2.0 * M_PI, 2.0 * M_PI);
			phase_b_type2 = fmod(run_data->psi_n[n] + run_data->phi_n[n + 1] + run_data->psi_n[n + 2] + 2.0 * M_PI, 2.0 * M_PI);
			phase_b_type3 = fmod(run_data->psi_n[n] + run_data->psi_n[n + 1] + run_data->phi_n[n + 2] + 2.0 * M_PI, 2.0 * M_PI);
			#endif
			
			
			///--------------------- Record the triad phases
			phase_sync->triads_u[i] = phase_u;
			#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
			phase_sync->triads_b[0 * phase_sync->num_triads + i] = phase_b_type1;
			phase_sync->triads_b[1 * phase_sync->num_triads + i] = phase_b_type2;
			phase_sync->triads_b[2 * phase_sync->num_triads + i] = phase_b_type3;
			#endif

			///--------------------- Record the triad phase order parameters
			#if defined(__PHASE_SYNC)
			phase_sync->triad_u_order[i] += cexp(I * phase_u);
			#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
			phase_sync->triad_b_order[0 * phase_sync->num_triads + i] += cexp(I * phase_b_type1);
			phase_sync->triad_b_order[1 * phase_sync->num_triads + i] += cexp(I * phase_b_type2);
			phase_sync->triad_b_order[2 * phase_sync->num_triads + i] += cexp(I * phase_b_type3);
			#endif
			#endif

			///--------------------- Record Phase Sync Stats
			#if defined(__PHASE_SYNC_STATS)
			gsl_status = gsl_histogram_increment(phase_sync->triad_u_hist[i], phase_u);
			if (gsl_status != 0) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update bin count for ["CYAN"%s"RESET"] for Iter ["CYAN"%ld"RESET"] -- GSL Exit Status [Err:"CYAN" %d"RESET" - Val:"CYAN" %lf"RESET"]\n-->> Exiting!!!\n", "Velocity Triads Histogram", iter, gsl_status, phase_u);
				exit(1);
			}
			#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
			gsl_status = gsl_histogram_increment(phase_sync->triad_b_hist[0 * phase_sync->num_triads + i], phase_b_type1);
			if (gsl_status != 0) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update bin count for ["CYAN"%s"RESET"] for Iter ["CYAN"%ld"RESET"] -- GSL Exit Status [Err:"CYAN" %d"RESET" - Val:"CYAN" %lf"RESET"]\n-->> Exiting!!!\n", "Magnetic Triads Histogram 1", iter, gsl_status, phase_b_type1);
				exit(1);
			}
			gsl_status = gsl_histogram_increment(phase_sync->triad_b_hist[1 * phase_sync->num_triads + i], phase_b_type2);
			if (gsl_status != 0) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update bin count for ["CYAN"%s"RESET"] for Iter ["CYAN"%ld"RESET"] -- GSL Exit Status [Err:"CYAN" %d"RESET" - Val:"CYAN" %lf"RESET"]\n-->> Exiting!!!\n", "Magnetic Triads Histogram 2", iter, gsl_status, phase_b_type2);
				exit(1);
			}
			gsl_status = gsl_histogram_increment(phase_sync->triad_b_hist[2 * phase_sync->num_triads + i], phase_b_type3);
			if (gsl_status != 0) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update bin count for ["CYAN"%s"RESET"] for Iter ["CYAN"%ld"RESET"] -- GSL Exit Status [Err:"CYAN" %d"RESET" - Val:"CYAN" %lf"RESET"]\n-->> Exiting!!!\n", "Magnetic Triads Histogram 3", iter, gsl_status, phase_b_type3);
				exit(1);
			}
			#endif
			#endif
		}

		///--------------------------------- Compute the phase differences and phase difference order parameters
		if (i < num_phase_diff) {

			// Get the phase differences
			phase_u_1 = fmod(run_data->phi_n[n] + 2.0 * M_PI, 2.0 * M_PI);
			phase_u_2 = fmod(run_data->phi_n[n + 3] + 2.0 * M_PI, 2.0 * M_PI);
			phase_u   = fmod((phase_u_1 - phase_u_2) + 2.0 * M_PI, 2.0 * M_PI);			
			#if defined(__MAGNETO)
			phase_b_1     = fmod(run_data->psi_n[n] + 2.0 * M_PI, 2.0 * M_PI);
			phase_b_2     = fmod(run_data->psi_n[n + 3] + 2.0 * M_PI, 2.0 * M_PI);
			phase_b_type1 = fmod((phase_b_1 - phase_b_2) + 2.0 * M_PI, 2.0 * M_PI);
			#endif
			
			
			///--------------------- Record the phase differences
			phase_sync->phase_diff_u[i] = phase_u;
			#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
			phase_sync->phase_diff_b[i] = phase_b_type1;
			#endif

			///--------------------- Record the phase difference order parameters
			#if defined(__PHASE_SYNC)
			phase_sync->phase_diff_u_order[i] += cexp(I * phase_u);
			#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
			phase_sync->phase_diff_b_order[i] += cexp(I * phase_b_type1);
			#endif
			#endif

			///--------------------- Record Phase Sync Stats
			#if defined(__PHASE_SYNC_STATS)
			gsl_status = gsl_histogram_increment(phase_sync->phase_diff_u_hist[i], phase_u);
			if (gsl_status != 0) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update bin count for ["CYAN"%s"RESET"] for Iter ["CYAN"%ld"RESET"] -- GSL Exit Status [Err:"CYAN" %d"RESET" - Val:"CYAN" %lf"RESET"]\n-->> Exiting!!!\n", "Velocity Phase Difference Histogram", iter, gsl_status, phase_u);
				exit(1);
			}
			#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
			gsl_status = gsl_histogram_increment(phase_sync->phase_diff_b_hist[i], phase_b_type1);
			if (gsl_status != 0) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update bin count for ["CYAN"%s"RESET"] for Iter ["CYAN"%ld"RESET"] -- GSL Exit Status [Err:"CYAN" %d"RESET" - Val:"CYAN" %lf"RESET"]\n-->> Exiting!!!\n", "Magnetic Phase Difference Histogram", iter, gsl_status, phase_b_type1);
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
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
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
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
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
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
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
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
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
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
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
		gsl_status = gsl_histogram_set_ranges_uniform(phase_sync->triad_b_hist[0 * phase_sync->num_triads + i], -0.0 - 0.05, 2.0 * M_PI + 0.05);
		if (gsl_status != 0) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to set bin ranges for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Triads Histogram");
			exit(1);
		}
		gsl_status = gsl_histogram_set_ranges_uniform(phase_sync->triad_b_hist[1 * phase_sync->num_triads + i], -0.0 - 0.05, 2.0 * M_PI + 0.05);
		if (gsl_status != 0) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to set bin ranges for ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Magnetic Triads Histogram");
			exit(1);
		}
		gsl_status = gsl_histogram_set_ranges_uniform(phase_sync->triad_b_hist[2 * phase_sync->num_triads + i], -0.0 - 0.05, 2.0 * M_PI + 0.05);
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
 * Wrapper Function That Writes the Phase Sync Stats Data at the end of the simulation to file
 */
void WritePhaseSyncStatsToFile(void) {

	const hsize_t D1 = 1;
	hsize_t dims1D[D1];
	const hsize_t D2 = 2;
	hsize_t dims2D[D2];
	#if (defined(__MAGNETO) || defined(__ELSASSAR_MHD))
	const hsize_t D3 = 3;
	hsize_t dims3D[D3];
	#endif

	// -------------------------------
	// Write Velocity Triads
	// -------------------------------
	///--------------------- Write Velocity triads histogram
	double* vel_triads_counts = (double* )malloc(sizeof(double) * phase_sync->num_triads * NUM_PHASE_SYNC_HIST_BINS);
	double* vel_triads_ranges = (double* )malloc(sizeof(double) * (NUM_PHASE_SYNC_HIST_BINS + 1));
	for (int i = 0; i < NUM_PHASE_SYNC_HIST_BINS + 1; ++i) {
		vel_triads_ranges[i] = phase_sync->triad_u_hist[0]->range[i];
		if (i < NUM_PHASE_SYNC_HIST_BINS) {
			for (int n = 0; n < phase_sync->num_triads; ++n) {
				vel_triads_counts[n * NUM_PHASE_SYNC_HIST_BINS + i] = phase_sync->triad_u_hist[n]->bin[i];
			}
		}
	}


	// Write Bin Count data 
	dims2D[0] = phase_sync->num_triads;
	dims2D[1] = NUM_PHASE_SYNC_HIST_BINS;
	if ( (H5LTmake_dataset(file_info->phase_sync_file_handle, "VelTriads_Counts", D2, dims2D, H5T_NATIVE_DOUBLE, vel_triads_counts)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelTriads_Counts");
	}
	// Write Bin Range data 
	dims1D[0] = NUM_PHASE_SYNC_HIST_BINS + 1;
	if ( (H5LTmake_dataset(file_info->phase_sync_file_handle, "VelTriads_Ranges", D1, dims1D, H5T_NATIVE_DOUBLE, vel_triads_ranges)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelTriads_Ranges");
	}

	// Free temporary memory
	free(vel_triads_counts);
	free(vel_triads_ranges);


	// -------------------------------
	// Write Velocity Phase Diff
	// -------------------------------
	///--------------------- Write Velocity phase difference histogram
	double* vel_phase_diff_counts = (double* )malloc(sizeof(double) * phase_sync->num_phase_diff * NUM_PHASE_SYNC_HIST_BINS);
	double* vel_phase_diff_ranges = (double* )malloc(sizeof(double) * (NUM_PHASE_SYNC_HIST_BINS + 1));
	for (int i = 0; i < NUM_PHASE_SYNC_HIST_BINS + 1; ++i) {
		vel_triads_ranges[i] = phase_sync->phase_diff_u_hist[0]->range[i];
		if (i < NUM_PHASE_SYNC_HIST_BINS) {
			for (int n = 0; n < phase_sync->num_phase_diff; ++n) {
				vel_phase_diff_counts[n * NUM_PHASE_SYNC_HIST_BINS + i] = phase_sync->phase_diff_u_hist[n]->bin[i];
			}
		}
	}

	printf("\n\nBEFROE HERERERERE\n");
	
	// Write Bin Count data 
	dims2D[0] = phase_sync->num_phase_diff;
	dims2D[1] = NUM_PHASE_SYNC_HIST_BINS;
	if ( (H5LTmake_dataset(file_info->phase_sync_file_handle, "VelPhaseDifference_Counts", D2, dims2D, H5T_NATIVE_DOUBLE, vel_phase_diff_counts)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelPhaseDifference_Counts");
	}
	// Write Bin Range data 
	dims1D[0] = NUM_PHASE_SYNC_HIST_BINS + 1;
	if ( (H5LTmake_dataset(file_info->phase_sync_file_handle, "VelPhaseDifference_Ranges", D1, dims1D, H5T_NATIVE_DOUBLE, vel_phase_diff_ranges)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelPhaseDifference_Ranges");
	}

	// Free temporary memory
	free(vel_phase_diff_counts);
	free(vel_phase_diff_ranges);

	printf("\n\nHEREERE\n");


	// -------------------------------
	// Write Magnetic Phase Diff
	// -------------------------------
	#if (defined(__MAGNETO) || defined(__ELSASSAR_MHD))
	///--------------------- Write Mag triads histogram
	double* mag_triads_counts = (double* )malloc(sizeof(double) * phase_sync->num_triads * NUM_PHASE_SYNC_HIST_BINS * NUM_MAG_TRIAD_TYPES);
	double* mag_triads_ranges = (double* )malloc(sizeof(double) * (NUM_PHASE_SYNC_HIST_BINS + 1));
	for (int i = 0; i < NUM_PHASE_SYNC_HIST_BINS + 1; ++i) {
		mag_triads_ranges[i] = phase_sync->triad_b_hist[0]->range[i];
		if (i < NUM_PHASE_SYNC_HIST_BINS) {
			for (int n = 0; n < phase_sync->num_triads; ++n) {
				for (int type = 0; type < NUM_MAG_TRIAD_TYPES; ++type)	{
					mag_triads_counts[NUM_MAG_TRIAD_TYPES * (n * NUM_PHASE_SYNC_HIST_BINS + i) + type] = phase_sync->triad_b_hist[type * phase_sync->num_triads + n]->bin[i];
				}
			}
		}
	}

	// Write Bin Count data 
	dims3D[0] = phase_sync->num_triads;
	dims3D[1] = NUM_PHASE_SYNC_HIST_BINS;
	dims3D[2] = NUM_MAG_TRIAD_TYPES;
	if ( (H5LTmake_dataset(file_info->phase_sync_file_handle, "MagTriads_Counts", D3, dims3D, H5T_NATIVE_DOUBLE, mag_triads_counts)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagTriads_Counts");
	}
	// Write Bin Range data 
	dims1D[0] = NUM_PHASE_SYNC_HIST_BINS + 1;
	if ( (H5LTmake_dataset(file_info->phase_sync_file_handle, "MagTriads_Ranges", D1, dims1D, H5T_NATIVE_DOUBLE, mag_triads_ranges)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagTriads_Ranges");
	}

	// Free temporary memory
	free(mag_triads_counts);
	free(mag_triads_ranges);


	// -------------------------------
	// Write Magnetic Phase Diff
	// -------------------------------
	///--------------------- Write Mag phase difference histogram
	double* mag_phase_diff_counts = (double* )malloc(sizeof(double) * phase_sync->num_phase_diff * NUM_PHASE_SYNC_HIST_BINS);
	double* mag_phase_diff_ranges = (double* )malloc(sizeof(double) * (NUM_PHASE_SYNC_HIST_BINS + 1));
	for (int i = 0; i < NUM_PHASE_SYNC_HIST_BINS + 1; ++i) {
		mag_triads_ranges[i] = phase_sync->phase_diff_b_hist[0]->range[i];
		if (i < NUM_PHASE_SYNC_HIST_BINS) {
			for (int n = 0; n < phase_sync->num_phase_diff; ++n) {
				mag_phase_diff_counts[n * NUM_PHASE_SYNC_HIST_BINS + i] = phase_sync->phase_diff_b_hist[n]->bin[i];
			}
		}
	}

	// Write Bin Count data 
	dims2D[0] = phase_sync->num_phase_diff;
	dims2D[1] = NUM_PHASE_SYNC_HIST_BINS;
	if ( (H5LTmake_dataset(file_info->phase_sync_file_handle, "MagPhaseDifference_Counts", D2, dims2D, H5T_NATIVE_DOUBLE, mag_phase_diff_counts)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagPhaseDifference_Counts");
	}
	// Write Bin Range data 
	dims1D[0] = NUM_PHASE_SYNC_HIST_BINS + 1;
	if ( (H5LTmake_dataset(file_info->phase_sync_file_handle, "MagPhaseDifference_Ranges", D1, dims1D, H5T_NATIVE_DOUBLE, mag_phase_diff_ranges)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagPhaseDifference_Ranges");
	}

	// Free temporary memory
	free(mag_phase_diff_counts);
	free(mag_phase_diff_ranges);
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
		#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
		gsl_histogram_free(phase_sync->triad_b_hist[0 * phase_sync->num_triads + i]);
		gsl_histogram_free(phase_sync->triad_b_hist[1 * phase_sync->num_triads + i]);
		gsl_histogram_free(phase_sync->triad_b_hist[2 * phase_sync->num_triads + i]);
		#endif
		if (i < phase_sync->num_phase_diff) {
			gsl_histogram_free(phase_sync->phase_diff_u_hist[i]);
			#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
			gsl_histogram_free(phase_sync->phase_diff_b_hist[i]);
			#endif
		}
	}
	#endif
}
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------
