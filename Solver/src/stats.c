/**
* @file utils.c  
* @author Enda Carroll
* @date Sept 2022
* @brief File containing the stats functions for the solver
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
/**
 * Function to compute the statistics for the current iteration
 * @param iters          The current iteration of the simulation
 * @param save_data_indx The current index for saving data to file
 */
void ComputeStats(const long int iters, const long int save_data_indx) {


	// Initialize variables
	int n;
	int gsl_status;
	#if defined(__VEL_HIST)
	double std = 10;
	double min = -10;
	double max = 10;
	#endif
	const long int N = sys_vars->N;
	double vel_enrg_flux_term, vel_hel_flux_term, vel_trip_prod_term;
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	double mag_enrg_flux_term, mag_hel_flux_term, mag_trip_prod_term;
	#endif

	// ------------------------------------
    // Get Fourier Fields
    // ------------------------------------	
	#if defined(PHASE_ONLY) || defined(PHASE_ONLY_FXD_AMP) || defined(__ELSASSAR_MHDL)
	for (int i = 0; i < N; ++i) {
		// Get temp indx
		n = i + 2;

		// Get the Fourier fields
		#if defined(__ELSASSAR_MHD)
		run_data->u[n] = (run_data->z_plus[n] + run_data->z_minus[n]) / 2.0;
		run_data->b[n] = (run_data->z_plus[n] - run_data->z_minus[n]) / 2.0;
		#elif defined(PHASE_ONLY) || defined(PHASE_ONLY_FXD_AMP)
		run_data->u[n] = run_data->a_n[n] * cexp(I * run_data->phi_n[n]);
		#if defined(__MAGNETO) 
		run_data->b[n] = run_data->b_n[n] * cexp(I * run_data->psi_n[n]);
		#endif
		#endif
	}
	#endif

	// ------------------------------------
    // Check System for Stationarity
    // ------------------------------------
    ///---------------------------------- System is not stationary yet -> Compute histogram limits
	if (iters < sys_vars->trans_iters) {
    	for (int i = 0; i < N; ++i) {
			// Get temp indx
			n = i + 2;

			// Update the running sums for the field stats
	    	gsl_status = gsl_rstat_add(creal(run_data->u[n]), stats_data->real_vel_moments[i]);
	    	if (gsl_status != 0) {
	    		// Print error message to error stream
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update running stats counter ["CYAN"%s"RESET"] - Error Value: ["CYAN"%d"RESET"] Field Value: ["CYAN"%lf"RESET"] \n-->> Exiting!!!\n", "Velocity Stats", gsl_status, creal(run_data->u[n]));

				// Save current state of system to file before exiting
				FinalWriteAndCloseOutputFile(N, iters, save_data_indx);
				
				// Exit programme
				exit(1);
			}
			
	    	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	    	gsl_status = gsl_rstat_add(creal(run_data->b[n]), stats_data->real_mag_moments[i]);
		    if (gsl_status != 0) {
		    	// Print error message to error stream
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update running stats counter ["CYAN"%s"RESET"] - Error Value: ["CYAN"%d"RESET"] Field Value: ["CYAN"%lf"RESET"] \n-->> Exiting!!!\n", "Magnetic Stats", gsl_status, creal(run_data->b[n]));

				// Save current state of system to file before exiting
				FinalWriteAndCloseOutputFile(N, iters, save_data_indx);
				
				// Exit programme
				exit(1);
			}
	    	#endif
	   	}
	}
    ///---------------------------------- System is stationary -> calculate stats
	else {
		// ---------------------------------------------
	    // Update Stats Counter & Reset Stats Objects
	    // --------------------------------------------
	    // Update counter
		stats_data->num_stats_steps++;

		// Set histogram bin limits then reset running stats counters
		if (stats_data->set_stats_flag) {
			for (int i = 0; i < N; ++i) {

				///-------------------------------------------- Set Velocity Histogram Bin ranges
				#if defined(__VEL_HIST)
				if (sys_vars->TRANS_ITERS_FLAG == TRANSIENT_ITERS) {
					// Get the standard deviation / rms
					std = gsl_rstat_sd(stats_data->real_vel_moments[i]);

					// Set bin ranges for the histograms -> Set (in units) of standard deviations
					gsl_status = gsl_histogram_set_ranges_uniform(stats_data->real_vel_hist[i], -VEL_BIN_LIM * std, VEL_BIN_LIM * std);
					if (gsl_status != 0) {
						fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to set bin ranges for: ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Real Velocity Histogram");
						exit(1);
					}
				}
				else {
					// Set bin ranges for the histograms -> Set (in units) of standard deviations
					gsl_status = gsl_histogram_set_ranges_uniform(stats_data->real_vel_hist[i], min - 0.05 * fabs(min), max + 0.05 * fabs(max));
					if (gsl_status != 0) {
						fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to set bin ranges for: ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Real Velocity Histogram");
						exit(1);
					}
				}
				#endif



				///-------------------------------------------- Set Magnetic Histogram Bin ranges				
				#if defined(__MAG_HIST) && (defined(__MAGNETO) || defined(__ELSASSAR_MHD))
				if (sys_vars->TRANS_ITERS_FLAG == TRANSIENT_ITERS) {
					// Get the standard deviation / rms
					std = gsl_rstat_sd(stats_data->real_mag_moments[i]);

					// Set bin ranges for the histograms -> Set (in units) of standard deviations
					gsl_status = gsl_histogram_set_ranges_uniform(stats_data->real_mag_hist[i], -VEL_BIN_LIM * std, VEL_BIN_LIM * std);
					if (gsl_status != 0) {
						fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to set bin ranges for: ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Real Magnetic Histogram");
						exit(1);
					}
				}
				else {
					// Set bin ranges for the histograms -> Set (in units) of standard deviations
					gsl_status = gsl_histogram_set_ranges_uniform(stats_data->real_mag_hist[i], min - 0.05 * fabs(min), max + 0.05 * fabs(max));
					if (gsl_status != 0) {
						fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to set bin ranges for: ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Real Magnetic Histogram");
						exit(1);
					}
				}
				#endif


				// Reset stats counters
				gsl_status = gsl_rstat_reset(stats_data->real_vel_moments[i]);
				#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
				gsl_status = gsl_rstat_reset(stats_data->real_mag_moments[i]);
				#endif
			}

			// Reset set stats flag
			stats_data->set_stats_flag = 0;
		}

		// ------------------------------------
	    // Compute Stats
	    // ------------------------------------
	    for (int i = 0; i < N; ++i) {
	    	n = i + 2;

	    	///---------------------------------- Velocity Histogram
	    	#if defined(__VEL_HIST)
	    	// Update Real velocity histogram
	    	gsl_status = gsl_histogram_increment(stats_data->real_vel_hist[i], creal(run_data->u[n]));
	    	if (gsl_status != 0) {
	    		// Print error message to error stream
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to increment histogram ["CYAN"%s"RESET"] - Error Value: ["CYAN"%d"RESET"] Field Value: ["CYAN"%lf"RESET"] Bin Limits: ["CYAN"%lf"RESET", "CYAN"%lf"RESET"]\n-->> Exiting!!!\n", "Real Velocity", gsl_status, creal(run_data->u[n]), gsl_histogram_min(stats_data->real_vel_hist[i]), gsl_histogram_max(stats_data->real_vel_hist[i]));

				// Save current state of system to file before exiting
				FinalWriteAndCloseOutputFile(N, iters, save_data_indx);
				
				// Exit programme
				exit(1);
			}
			#endif

			///---------------------------------- Magnetic Histogram
	    	#if defined(__MAG_HIST) && (defined(__MAGNETO) || defined(__ELSASSAR_MHD))
	    	// Update Real velocity histogram
	    	gsl_status = gsl_histogram_increment(stats_data->real_mag_hist[i], creal(run_data->b[n]));
	    	if (gsl_status != 0) {
	    		// Print error message to error stream
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to increment histogram ["CYAN"%s"RESET"] - Error Value: ["CYAN"%d"RESET"] Field Value: ["CYAN"%lf"RESET"] Bin Limits: ["CYAN"%lf"RESET", "CYAN"%lf"RESET"]\n-->> Exiting!!!\n", "Real Velocity", gsl_status, creal(run_data->u[n]), gsl_histogram_min(stats_data->real_mag_hist[i]), gsl_histogram_max(stats_data->real_mag_hist[i]));

				// Save current state of system to file before exiting
				FinalWriteAndCloseOutputFile(N, iters, save_data_indx);
				
				// Exit programme
				exit(1);
			}
			#endif

	    	///---------------------------------- Running Stats
	    	// Update the running sums for the field stats
	    	gsl_status = gsl_rstat_add(creal(run_data->u[n]), stats_data->real_vel_moments[i]);
	    	if (gsl_status != 0) {
	    		// Print error message to error stream
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update running stats counter ["CYAN"%s"RESET"] - Error Value: ["CYAN"%d"RESET"] Field Value: ["CYAN"%lf"RESET"] \n-->> Exiting!!!\n", "Velocity Stats", gsl_status, creal(run_data->u[n]));

				// Save current state of system to file before exiting
				FinalWriteAndCloseOutputFile(N, iters, save_data_indx);
				
				// Exit programme
				exit(1);
			}
	    	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	    	gsl_status = gsl_rstat_add(creal(run_data->b[n]), stats_data->real_mag_moments[i]);
		    if (gsl_status != 0) {
		    	// Print error message to error stream
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update running stats counter ["CYAN"%s"RESET"] - Error Value: ["CYAN"%d"RESET"] Field Value: ["CYAN"%lf"RESET"] \n-->> Exiting!!!\n", "Magnetic Stats", gsl_status, creal(run_data->b[n]));

				// Save current state of system to file before exiting
				FinalWriteAndCloseOutputFile(N, iters, save_data_indx);
				
				// Exit programme
				exit(1);
			}
	    	#endif
	    	
	    	///---------------------------------- Structure Functions
	    	// Compute the flux terms
	    	#if defined(__STR_FUNC_VEL_FLUX) || defined(__STR_FUNC_MAG_FLUX) || defined(__STR_FUNC_TRIP_PROD_VEL) || defined(__STR_FUNC_TRIP_PROD_MAG)
	    	// Tripple product term
	    	vel_trip_prod_term = cimag(run_data->u[n + 1] * run_data->u[n] * run_data->u[n - 1]);

	    	// Flux terms
			vel_enrg_flux_term = cimag(run_data->u[n + 2] * run_data->u[n + 1] * run_data->u[n] + (1.0 - sys_vars->EPS) / sys_vars->Lambda * run_data->u[n + 1] * run_data->u[n] * run_data->u[n - 1]);
			vel_hel_flux_term  = cimag(run_data->u[n + 2] * run_data->u[n + 1] * run_data->u[n] - (sys_vars->EPS * sys_vars->Lambda + 1.0) / pow(sys_vars->Lambda, 2.0) * run_data->u[n + 1] * run_data->u[n] * run_data->u[n - 1]);
			#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
			// Tripple product term
			mag_trip_prod_term = cimag(run_data->b[n + 1] * run_data->b[n] * run_data->b[n - 1]);

			// Flux product term
			mag_enrg_flux_term = cimag(run_data->b[n + 2] * run_data->b[n + 1] * run_data->b[n] + (1.0 - sys_vars->EPS) / sys_vars->Lambda * run_data->b[n + 1] * run_data->b[n] * run_data->b[n - 1]);
			mag_hel_flux_term  = cimag(run_data->b[n + 2] * run_data->b[n + 2] * run_data->b[n + 2]);
			#endif    	
			#endif

	    	// Compute the moments for the structure functions
	    	for (int p = 1; p <= NUM_POW; ++p) {
	    		
	    		// Compute the moments of the fields
	    		#if defined(__STR_FUNC_MAG) || defined(__STR_FUNC_VEL)
				stats_data->vel_str_func[p - 1][i] += pow(cabs(run_data->u[n]), p);
				#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
				stats_data->mag_str_func[p - 1][i] += pow(cabs(run_data->b[n]), p); 
	    		#endif
	    		#endif

				// Compute the moments of the fluxes
				#if defined(__STR_FUNC_VEL_FLUX) || defined(__STR_FUNC_MAG_FLUX)
				// Tripple product str func
				stats_data->vel_trip_prod_str_func[p - 1][i] += pow(sgn(vel_trip_prod_term), p) * pow(fabs(vel_trip_prod_term), p / 3.0);
				stats_data->vel_trip_prod_str_func_abs[p - 1][i] += pow(fabs(vel_trip_prod_term), p / 3.0);

				// Flux term str func
				stats_data->vel_flux_str_func[0][p - 1][i] += pow(sgn(vel_enrg_flux_term), p) * pow(fabs(vel_enrg_flux_term), p / 3.0);
				stats_data->vel_flux_str_func[1][p - 1][i] += pow(sgn(vel_hel_flux_term), p) * pow(fabs(vel_hel_flux_term), p / 3.0);
				stats_data->vel_flux_str_func_abs[0][p - 1][i] += pow(fabs(vel_enrg_flux_term), p / 3.0);
				stats_data->vel_flux_str_func_abs[1][p - 1][i] += pow(fabs(vel_hel_flux_term), p / 3.0);
	    		#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	    		// Tripple product str func
				stats_data->mag_trip_prod_str_func[p - 1][i] += pow(sgn(mag_trip_prod_term), p) * pow(fabs(mag_trip_prod_term), p / 3.0);
				stats_data->mag_trip_prod_str_func_abs[p - 1][i] += pow(fabs(mag_trip_prod_term), p / 3.0);

				// Flux term str func
				stats_data->mag_flux_str_func[0][p - 1][i] += pow(sgn(mag_enrg_flux_term), p) * pow(fabs(mag_enrg_flux_term), p / 3.0);
				stats_data->mag_flux_str_func[1][p - 1][i] += pow(sgn(mag_hel_flux_term), p) * pow(fabs(mag_hel_flux_term), p / 3.0);
				stats_data->mag_flux_str_func_abs[0][p - 1][i] += pow(fabs(mag_enrg_flux_term), p / 3.0);
				stats_data->mag_flux_str_func_abs[1][p - 1][i] += pow(fabs(mag_hel_flux_term), p / 3.0);
	    		#endif
	    		#endif
	    	}
	    }
	}
	
}
/**
 * Function to allocate and initialize the stats objects
 */
void InitializeStats(void) {

	// Initialize variables
	const long int N = sys_vars->N;
	int gsl_status;

	// Initialize the stats counter
	stats_data->num_stats_steps = 0;

	// Set stats setting flag
	stats_data->set_stats_flag = 1;


	// ------------------------------------
    // Allocate & Initialize Stats Objects
    // ------------------------------------
	///--------------------------------- Structure Functions
	// Allocate memory for the structure functions
	for (int i = 0; i < NUM_POW; ++i) {
		#if defined(__STR_FUNC_VEL)
		stats_data->vel_str_func[i] = (double* )malloc(sizeof(double) * N);
		if (stats_data->vel_str_func[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Velocity Structure Function");
			exit(1);
		}
		#endif
		#if defined(__STR_FUNC_TRIP_PROD_VEL)
		stats_data->vel_trip_prod_str_func[i] = (double* )malloc(sizeof(double) * N);
		if (stats_data->vel_trip_prod_str_func[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Velocity Tripple Product Structure Function");
			exit(1);
		}
		stats_data->vel_trip_prod_str_func_abs[i] = (double* )malloc(sizeof(double) * N);
		if (stats_data->vel_trip_prod_str_func_abs[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Absolute Velocity Tripple Product Structure Function");
			exit(1);
		}		
		#endif
		#if defined(__STR_FUNC_VEL_FLUX)
		for (int j = 0; j < 2; ++j) {
			stats_data->vel_flux_str_func[j][i] = (double* )malloc(sizeof(double) * N);
			if (stats_data->vel_flux_str_func[j][i] == NULL) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Velocity Flux Structure Function");
				exit(1);
			}
			stats_data->vel_flux_str_func_abs[j][i] = (double* )malloc(sizeof(double) * N);
			if (stats_data->vel_flux_str_func_abs[j][i] == NULL) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Absolute Velocity Flux Structure Function");
				exit(1);
			}		
		}
		#endif
		#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
		#if defined(__STR_FUNC_MAG)
		stats_data->mag_str_func[i] = (double* )malloc(sizeof(double) * N);
		if (stats_data->mag_str_func[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Magnetic Structure Function");
			exit(1);
		}
		#endif
		#if defined(__STR_FUNC_TRIP_PROD_MAG)
		stats_data->mag_trip_prod_str_func[i] = (double* )malloc(sizeof(double) * N);
		if (stats_data->mag_trip_prod_str_func[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Magnetic Tripple Product Structure Function");
			exit(1);
		}
		stats_data->mag_trip_prod_str_func_abs[i] = (double* )malloc(sizeof(double) * N);
		if (stats_data->mag_trip_prod_str_func_abs[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Absolute Magnetic Tripple Product Structure Function");
			exit(1);
		}		
		#endif
		#if defined(__STR_FUNC_MAG_FLUX)
		for (int j = 0; j < 2; ++j) {
			stats_data->mag_flux_str_func[j][i] = (double* )malloc(sizeof(double) * N);
			if (stats_data->mag_flux_str_func[j][i] == NULL) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Magnetic Flux Sturcure Function");
				exit(1);
			}
			stats_data->mag_flux_str_func_abs[j][i] = (double* )malloc(sizeof(double) * N);
			if (stats_data->mag_flux_str_func_abs[j][i] == NULL) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Absolute Magnetic Flux Sturcure Function");
				exit(1);
			}
		}		
		#endif
		#endif
	}


	///--------------------------------- Running Stats Objects
	// Allocate memory for the running stats objects
	stats_data->real_vel_moments = (double** )malloc(sizeof(double* ) * N);
	if (stats_data->real_vel_moments == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Velocity Stats");
		exit(1);
	}
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	stats_data->real_mag_moments = (double** )malloc(sizeof(double* ) * N);
	if (stats_data->real_mag_moments == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Magnetic Stats");
		exit(1);
	}
	#endif

	// Initialize the running stats objects
	for (int i = 0; i < N; ++i) {
		stats_data->real_vel_moments[i] = gsl_rstat_alloc();
		if (stats_data->real_vel_moments[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Running stats object: ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Velocity Stats");
			exit(1);
		}
		#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
		stats_data->real_mag_moments[i] = gsl_rstat_alloc();
		if (stats_data->real_mag_moments[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Running stats object: ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Magnetic Stats");
			exit(1);
		}
		#endif
	}

	///--------------------------------- Histograms
	#if defined(__VEL_HIST)
	// Allocate memory for the array of histogram structs
	stats_data->real_vel_hist = (double** )malloc(sizeof(double* ) * N);
	for (int i = 0; i < N; ++i) {
		// Allocate	the histogram structs
		stats_data->real_vel_hist[i] = gsl_histogram_alloc(VEL_NUM_BINS);
		if (stats_data->real_vel_hist[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for histogram struct for: ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Real Velocity Histogram");
			exit(1);
		}
	}
	#endif

	#if defined(__MAG_HIST) && (defined(__MAGNETO) || defined(__ELSASSAR_MHD))
	// Allocate memory for the array of histogram structs
	stats_data->real_mag_hist = (double** )malloc(sizeof(double* ) * N);
	for (int i = 0; i < N; ++i) {
		// Allocate	the histogram structs
		stats_data->real_mag_hist[i] = gsl_histogram_alloc(VEL_NUM_BINS);
		if (stats_data->real_vel_hist[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for histogram struct for: ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Real Magnetic Histogram");
			exit(1);
		}
	}
	#endif
}
/**
 * Wrapper Function for writing all stats data at the end of the simulation to file
 */
void WriteStatsToFile(void) {

	// Initialize variables
	herr_t status;
	static const hsize_t D1 = 1;
	hsize_t dims1D[D1];
	const hsize_t D2 = 2;
	hsize_t dims2D[D2];
	const hsize_t D3 = 3;
	hsize_t dims3D[D3];

	// ------------------------------------
    // Write Number of Stats Steps To File
    // ------------------------------------
	dims1D[0] = 1;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "NumStatsSteps", D1, dims1D, H5T_NATIVE_LONG, &(stats_data->num_stats_steps))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "NumStatsSteps");
	}

	// ------------------------------------
    // Write Velocity Field Data
    // ------------------------------------
	///-------------------------- Write velocity field 
	dims1D[0] = sys_vars->N;
	#if defined(PHASE_ONLY)
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "VelAmps", D1, dims1D, H5T_NATIVE_DOUBLE, &(run_data->a_n[2]))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelAmps");
	}
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "VelPhases", D1, dims1D, H5T_NATIVE_DOUBLE, &(run_data->phi_n[2]))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelPhases");
	}
	#elif defined(__ELSASSAR_MHD)
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "ZPlus", D1, dims1D, file_info->COMPLEX_DTYPE, &(run_data->z_plus[2]))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "ZPlus");
	}
	#else
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "VelModes", D1, dims1D, file_info->COMPLEX_DTYPE, &(run_data->u[2]))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelModes");
	}
	#endif


	// ------------------------------------
    // Write Velocity Field Stats Data
    // ------------------------------------
	///-------------------------- Velocity field stats
	double* tmp_vel_stats = (double* )malloc(sizeof(double) * (NUM_RUN_STATS) * sys_vars->N);
	for (int i = 0; i < sys_vars->N; ++i) {
		tmp_vel_stats[i * (NUM_RUN_STATS) + 0] = gsl_rstat_mean(stats_data->real_vel_moments[i]);
		tmp_vel_stats[i * (NUM_RUN_STATS) + 1] = gsl_rstat_variance(stats_data->real_vel_moments[i]);
		tmp_vel_stats[i * (NUM_RUN_STATS) + 2] = gsl_rstat_skew(stats_data->real_vel_moments[i]);
		tmp_vel_stats[i * (NUM_RUN_STATS) + 3] = gsl_rstat_kurtosis(stats_data->real_vel_moments[i]);		
		tmp_vel_stats[i * (NUM_RUN_STATS) + 4] = gsl_rstat_rms(stats_data->real_vel_moments[i]);
		tmp_vel_stats[i * (NUM_RUN_STATS) + 5] = gsl_rstat_min(stats_data->real_vel_moments[i]);
		tmp_vel_stats[i * (NUM_RUN_STATS) + 6] = gsl_rstat_max(stats_data->real_vel_moments[i]);
	}

	// Write data 
	dims2D[0] = sys_vars->N;
	dims2D[1] = NUM_RUN_STATS;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "VelStats", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_vel_stats)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelStats");
	}

	// Free temp memory
	free(tmp_vel_stats);



	// ------------------------------------
    // Write Velocity Field Hist Data
    // ------------------------------------
	///-------------------------- Real Velocity Histogram
	#if defined(__VEL_HIST)
	double* tmp_vel_hist_bin    = (double* )malloc(sizeof(double) * (VEL_NUM_BINS) * sys_vars->N);
	double* tmp_vel_hist_ranges = (double* )malloc(sizeof(double) * (VEL_NUM_BINS + 1) * sys_vars->N);
	for (int i = 0; i < sys_vars->N; ++i) {
		for (int j = 0; j < VEL_NUM_BINS + 1; ++j) {
			if (j < VEL_NUM_BINS) {
				tmp_vel_hist_bin[i * VEL_NUM_BINS + j] = stats_data->real_vel_hist[i]->bin[j];
			}
			tmp_vel_hist_ranges[i * (VEL_NUM_BINS + 1) + j] = stats_data->real_vel_hist[i]->range[j];
		}
	}

	// Write Bin Count data 
	dims2D[0] = sys_vars->N;
	dims2D[1] = VEL_NUM_BINS;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "RealVelHist_Counts", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_vel_hist_bin)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "RealVelHist_Counts");
	}
	// Write Bin Range data 
	dims2D[0] = sys_vars->N;
	dims2D[1] = VEL_NUM_BINS + 1;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "RealVelHist_Ranges", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_vel_hist_ranges)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "RealVelHist_Ranges");
	}

	// Free temporary memory
	free(tmp_vel_hist_bin);
	free(tmp_vel_hist_ranges);
	#endif


	// ------------------------------------
    // Write Velocity Str Func Data
    // ------------------------------------
	///-------------------------- Velocity structure function
	#if defined(__STR_FUNC_VEL)
	// Allocate temporary contiguous array
	double* tmp_vel_str = (double* )malloc(sizeof(double) * NUM_POW * sys_vars->N);
	for (int i = 0; i < sys_vars->N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_vel_str[i * NUM_POW + p] = stats_data->vel_str_func[p][i];
		}
	}

	// Write data 
	dims2D[0] = sys_vars->N;
	dims2D[1] = NUM_POW;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionVel", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_vel_str)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionVel");
	}

	// Free temp memory
	free(tmp_vel_str);
	#endif

	///-------------------------- Tripple product structure function
	#if defined(__STR_FUNC_TRIP_PROD_VEL) 
	// Allocate temporary contiguous array
	double* tmp_vel_trip_prod_str = (double* )malloc(sizeof(double) * NUM_POW * sys_vars->N);
	for (int i = 0; i < sys_vars->N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_vel_trip_prod_str[i * NUM_POW + p] = stats_data->vel_trip_prod_str_func[p][i];
		}
	}

	// Write data 
	dims2D[0] = sys_vars->N;
	dims2D[1] = NUM_POW;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionTrippleProdVel", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_vel_trip_prod_str)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionTrippleProdVel");
	}

	for (int i = 0; i < sys_vars->N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_vel_trip_prod_str[i * NUM_POW + p] = stats_data->vel_trip_prod_str_func_abs[p][i];
		}
	}

	// Write data 
	dims2D[0] = sys_vars->N;
	dims2D[1] = NUM_POW;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionTrippleProdVelAbs", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_vel_trip_prod_str)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionTrippleProdVelAbs");
	}

	// Free temp memory
	free(tmp_vel_trip_prod_str);
	#endif

	///-------------------------- Velocity Flux structure function
	#if defined(__STR_FUNC_VEL_FLUX)
	// Allocate temporary contiguous array
	double* tmp_vel_str_flux = (double* )malloc(sizeof(double) * 2 * NUM_POW * sys_vars->N);
	for (int i = 0; i < sys_vars->N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_vel_str_flux[2 * (i * NUM_POW + p) + 0] = stats_data->vel_flux_str_func[0][p][i];
			tmp_vel_str_flux[2 * (i * NUM_POW + p) + 1] = stats_data->vel_flux_str_func[1][p][i];
		}
	}

	// Write data 
	dims3D[0] = sys_vars->N;
	dims3D[1] = NUM_POW;
	dims3D[2] = 2;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionVelFlux", D3, dims3D, H5T_NATIVE_DOUBLE, tmp_vel_str_flux)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionVelFlux");
	}

	for (int i = 0; i < sys_vars->N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_vel_str_flux[2 * (i * NUM_POW + p) + 0] = stats_data->vel_flux_str_func_abs[0][p][i];
			tmp_vel_str_flux[2 * (i * NUM_POW + p) + 1] = stats_data->vel_flux_str_func_abs[1][p][i];
		}
	}

	// Write data 
	dims3D[0] = sys_vars->N;
	dims3D[1] = NUM_POW;
	dims3D[2] = 2;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionVelFluxAbs", D3, dims3D, H5T_NATIVE_DOUBLE, tmp_vel_str_flux)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionVelFluxAbs");
	}

	// Free temp memory
	free(tmp_vel_str_flux);
	#endif



	// ------------------------------------
    // Write Magnetic Field Data
    // ------------------------------------
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	///-------------------------- Write Magnetic field 
	dims1D[0] = sys_vars->N;
	#if defined(PHASE_ONLY)
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "MagAmps", D1, dims1D, H5T_NATIVE_DOUBLE, &(run_data->b_n[2]))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagAmps");
	}
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "MagPhases", D1, dims1D, H5T_NATIVE_DOUBLE, &(run_data->psi_n[2]))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagPhases");
	}
	#elif defined(__ELSASSAR_MHD)
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "ZMinus", D1, dims1D, file_info->COMPLEX_DTYPE, &(run_data->z_minus[2]))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "ZMinus");
	}
	#else
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "MagModes", D1, dims1D, file_info->COMPLEX_DTYPE, &(run_data->b[2]))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagModes");
	}
	#endif


	// ------------------------------------
    // Write Magnetic Field Stats Data
    // ------------------------------------
	///-------------------------- Magnetic field stats
	double* tmp_mag_stat = (double* )malloc(sizeof(double) * (NUM_RUN_STATS) * sys_vars->N);
	for (int i = 0; i < sys_vars->N; ++i) {
		tmp_mag_stat[i * (NUM_RUN_STATS) + 0] = gsl_rstat_mean(stats_data->real_mag_moments[i]);
		tmp_mag_stat[i * (NUM_RUN_STATS) + 1] = gsl_rstat_variance(stats_data->real_mag_moments[i]);
		tmp_mag_stat[i * (NUM_RUN_STATS) + 2] = gsl_rstat_skew(stats_data->real_mag_moments[i]);
		tmp_mag_stat[i * (NUM_RUN_STATS) + 3] = gsl_rstat_kurtosis(stats_data->real_mag_moments[i]);		
		tmp_mag_stat[i * (NUM_RUN_STATS) + 4] = gsl_rstat_rms(stats_data->real_mag_moments[i]);
		tmp_mag_stat[i * (NUM_RUN_STATS) + 5] = gsl_rstat_min(stats_data->real_mag_moments[i]);
		tmp_mag_stat[i * (NUM_RUN_STATS) + 6] = gsl_rstat_max(stats_data->real_mag_moments[i]);
	}

	// Write data 
	dims2D[0] = sys_vars->N;
	dims2D[1] = NUM_RUN_STATS;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "MagStats", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_mag_stat)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagStats");
	}

	// Free temp memory
	free(tmp_mag_stat);


	// ------------------------------------
    // Write Magnetic Field Hist Data
    // ------------------------------------
	///-------------------------- Real Magnetic Histogram
	#if defined(__MAG_HIST)
	double* tmp_mag_hist_bin    = (double* )malloc(sizeof(double) * (VEL_NUM_BINS) * sys_vars->N);
	double* tmp_mag_hist_ranges = (double* )malloc(sizeof(double) * (VEL_NUM_BINS + 1) * sys_vars->N);
	for (int i = 0; i < sys_vars->N; ++i) {
		for (int j = 0; j < VEL_NUM_BINS + 1; ++j) {
			if (j < VEL_NUM_BINS) {
				tmp_mag_hist_bin[i * VEL_NUM_BINS + j] = stats_data->real_mag_hist[i]->bin[j];
			}
			tmp_mag_hist_ranges[i * (VEL_NUM_BINS + 1) + j] = stats_data->real_mag_hist[i]->range[j];
		}
	}

	// Write Bin Count data 
	dims2D[0] = sys_vars->N;
	dims2D[1] = VEL_NUM_BINS;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "RealMagHist_Counts", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_mag_hist_bin)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "RealMagHist_Counts");
	}
	// Write Bin Range data 
	dims2D[0] = sys_vars->N;
	dims2D[1] = VEL_NUM_BINS + 1;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "RealMagHist_Ranges", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_mag_hist_ranges)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "RealMagHist_Ranges");
	}

	// Free temporary memory
	free(tmp_mag_hist_bin);
	free(tmp_mag_hist_ranges);
	#endif


	// ------------------------------------
    // Write Magnetic Str Func Data
    // ------------------------------------
	///-------------------------- Magnetic structure function
	#if defined(__STR_FUNC_MAG) 
	// Allocate temporary contiguous array
	double* tmp_mag_str = (double* )malloc(sizeof(double) * NUM_POW * sys_vars->N);
	for (int i = 0; i < sys_vars->N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_mag_str[i * NUM_POW + p] = stats_data->mag_str_func[p][i];
		}
	}

	// Write data 
	dims2D[0] = sys_vars->N;
	dims2D[1] = NUM_POW;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionMag", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_mag_str)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionMag");
	}

	// Free temp memory
	free(tmp_mag_str);
	#endif

	///-------------------------- Tripple product structure function
	#if defined(__STR_FUNC_TRIP_PROD_MAG) 
	// Allocate temporary contiguous array
	double* tmp_mag_trip_prod_str = (double* )malloc(sizeof(double) * NUM_POW * sys_vars->N);
	for (int i = 0; i < sys_vars->N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_mag_trip_prod_str[i * NUM_POW + p] = stats_data->mag_trip_prod_str_func[p][i];
		}
	}

	// Write data 
	dims2D[0] = sys_vars->N;
	dims2D[1] = NUM_POW;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionTrippleProdMag", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_mag_trip_prod_str)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionTrippleProdMag");
	}

	for (int i = 0; i < sys_vars->N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_mag_trip_prod_str[i * NUM_POW + p] = stats_data->mag_trip_prod_str_func_abs[p][i];
		}
	}

	// Write data 
	dims2D[0] = sys_vars->N;
	dims2D[1] = NUM_POW;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionTrippleProdMagAbs", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_mag_trip_prod_str)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionTrippleProdMagAbs");
	}

	// Free temp memory
	free(tmp_mag_trip_prod_str);
	#endif

	///-------------------------- Magnetic flux structure function
	#if defined(__STR_FUNC_MAG_FLUX)
	// Allocate temporary contiguous array
	double* tmp_mag_str_flux = (double* )malloc(sizeof(double) * 2 * NUM_POW * sys_vars->N);
	for (int i = 0; i < sys_vars->N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_mag_str_flux[2 * (i * NUM_POW + p) + 0] = stats_data->mag_flux_str_func[0][p][i];
			tmp_mag_str_flux[2 * (i * NUM_POW + p) + 1] = stats_data->mag_flux_str_func[1][p][i];
		}
	}

	// Write data 
	dims3D[0] = sys_vars->N;
	dims3D[1] = NUM_POW;
	dims3D[2] = 2;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionMagFlux", D3, dims3D, H5T_NATIVE_DOUBLE, tmp_mag_str_flux)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionMagFlux");
	}

	for (int i = 0; i < sys_vars->N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_mag_str_flux[2 * (i * NUM_POW + p) + 0] = stats_data->mag_flux_str_func_abs[0][p][i];
			tmp_mag_str_flux[2 * (i * NUM_POW + p) + 1] = stats_data->mag_flux_str_func_abs[1][p][i];
		}
	}

	// Write data 
	dims3D[0] = sys_vars->N;
	dims3D[1] = NUM_POW;
	dims3D[2] = 2;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionMagFluxAbs", D3, dims3D, H5T_NATIVE_DOUBLE, tmp_mag_str_flux)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionMagFluxAbs");
	}

	// Free temp memory
	free(tmp_mag_str_flux);
	#endif
	#endif
}
/**
 * Wrapper function to free all of the stats objects
 */
void FreeStatsObjects(void) {

	for (int i = 0; i < NUM_POW; ++i) {
		#if defined(__STR_FUNC_VEL)
		free(stats_data->vel_str_func[i]);
		#endif
		#if defined(__STR_FUNC_TRIP_PROD_VEL)
		free(stats_data->vel_trip_prod_str_func[i]);
		free(stats_data->vel_trip_prod_str_func_abs[i]);
		#endif
		#if defined(__STR_FUNC_VEL_FLUX)
		for (int j = 0; j < 2; ++j) {
			free(stats_data->vel_flux_str_func[j][i]);
			free(stats_data->vel_flux_str_func_abs[j][i]);
		}		
		#endif
		#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
		#if defined(__STR_FUNC_MAG)
		free(stats_data->mag_str_func[i]);
		#endif
		#if defined(__STR_FUNC_TRIP_PROD_MAG)
		free(stats_data->mag_trip_prod_str_func[i]);
		free(stats_data->mag_trip_prod_str_func_abs[i]);
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
		#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
		gsl_rstat_free(stats_data->real_mag_moments[i]);
		#endif
		#if defined(__VEL_HIST)
		gsl_histogram_free(stats_data->real_vel_hist[i]);
		#endif
		#if (defined(__MAGNETO) || defined(__ELSASSAR_MHD)) && defined(__MAG_HIST)
		gsl_histogram_free(stats_data->real_mag_hist[i]);
		#endif
	}
	free(stats_data->real_vel_moments);
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	free(stats_data->real_mag_moments);
	#endif
	#if defined(__VEL_HIST)
	free(stats_data->real_vel_hist);
	#endif
	#if (defined(__MAGNETO) || defined(__ELSASSAR_MHD)) && defined(__MAG_HIST)
	free(stats_data->real_mag_hist);
	#endif
}
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------
