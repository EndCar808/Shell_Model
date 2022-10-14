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
	double std, min, max;
	#endif
	const long int N = sys_vars->N;
	double vel_enrg_flux_term, vel_hel_flux_term;
	#if defined(__MAGNETO)
	double mag_enrg_flux_term, mag_hel_flux_term;
	#endif

	// ------------------------------------
    // Check System for Stationary
    // ------------------------------------
    ///---------------------------------- System is not stationary yet -> Compute histogram limits
	if (iters < sys_vars->trans_iters) {
    	for (int i = 0; i < N; ++i) {
			n = i + 2;

			// Update the running sums for the field stats
	    	gsl_status = gsl_rstat_add(creal(run_data->u[n]), stats_data->vel_moments[i]);
	    	if (gsl_status != 0) {
	    		// Print error message to error stream
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update running stats counter ["CYAN"%s"RESET"] - Error Value: ["CYAN"%d"RESET"] Field Value: ["CYAN"%lf"RESET"] \n-->> Exiting!!!\n", "Velocity Stats", gsl_status, creal(run_data->u[n]));

				// Save current state of system to file before exiting
				FinalWriteAndCloseOutputFile(N, iters, save_data_indx);
				
				// Exit programme
				exit(1);
			}
	    	#if defined(__MAGNETO)
	    	gsl_status = gsl_rstat_add(creal(run_data->b[n]), stats_data->mag_moments[i]);
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
				#if defined(__VEL_HIST)
				// Get the standard deviation / rms
				std = gsl_rstat_sd(stats_data->vel_moments[i]);

				// Set bin ranges for the histograms -> Set (in units) of standard deviations
				gsl_status = gsl_histogram_set_ranges_uniform(stats_data->real_vel_hist[i], -VEL_BIN_LIM * std, VEL_BIN_LIM * std);
				if (gsl_status != 0) {
					fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to set bin ranges for: ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Real Velocity Histogram");
					exit(1);
				}
				#endif

				// Reset stats counters
				gsl_status = gsl_rstat_reset(stats_data->vel_moments[i]);
				#if defined(__MAGNETO)
				gsl_status = gsl_rstat_reset(stats_data->mag_moments[i]);
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

	    	///---------------------------------- Histograms
	    	#if defined(__VEL_HIST)
	    	// Update Real velocity histogram
	    	gsl_status = gsl_histogram_increment(stats_data->real_vel_hist[i], creal(run_data->u[n]));
	    	if (gsl_status != 0) {
	    		// Print error message to error stream
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to increment histogram ["CYAN"%s"RESET"] - Error Value: ["CYAN"%d"RESET"] Field Value: ["CYAN"%lf"RESET"] \n-->> Exiting!!!\n", "Real Velocity", gsl_status, creal(run_data->u[n]));

				// Save current state of system to file before exiting
				FinalWriteAndCloseOutputFile(N, iters, save_data_indx);
				
				// Exit programme
				exit(1);
			}
			#endif

	    	///---------------------------------- Running Stats
	    	// Update the running sums for the field stats
	    	gsl_status = gsl_rstat_add(creal(run_data->u[n]), stats_data->vel_moments[i]);
	    	if (gsl_status != 0) {
	    		// Print error message to error stream
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to update running stats counter ["CYAN"%s"RESET"] - Error Value: ["CYAN"%d"RESET"] Field Value: ["CYAN"%lf"RESET"] \n-->> Exiting!!!\n", "Velocity Stats", gsl_status, creal(run_data->u[n]));

				// Save current state of system to file before exiting
				FinalWriteAndCloseOutputFile(N, iters, save_data_indx);
				
				// Exit programme
				exit(1);
			}
	    	#if defined(__MAGNETO)
	    	gsl_status = gsl_rstat_add(creal(run_data->b[n]), stats_data->mag_moments[i]);
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
	    	#if defined(__STR_FUNC_VEL_FLUX) || defined(__STR_FUNC_MAG_FLUX)
			vel_enrg_flux_term = cimag(run_data->u[n + 2] * run_data->u[n + 1] * run_data->u[n] + (1.0 - sys_vars->EPS) / sys_vars->Lambda * run_data->u[n + 1] * run_data->u[n] * run_data->u[n - 1]);
			vel_hel_flux_term  = cimag(run_data->u[n + 2] * run_data->u[n + 1] * run_data->u[n] - (sys_vars->EPS * sys_vars->Lambda + 1.0) / pow(sys_vars->Lambda, 2.0) * run_data->u[n + 1] * run_data->u[n] * run_data->u[n - 1]);
			#if defined(__MAGNETO) 
			// mag_enrg_flux_term = cimag(run_data->u[n + 2] * run_data->u[n + 2] * run_data->u[n + 2])
			// mag_hel_flux_term  = cimag(run_data->u[n + 2] * run_data->u[n + 2] * run_data->u[n + 2])
			#endif    	
			#endif

	    	// Compute the moments for the structure functions
	    	for (int p = 1; p <= NUM_POW; ++p) {
	    		
	    		// Compute the moments of the fields
	    		#if defined(__STR_FUNC_MAG) || defined(__STR_FUNC_VEL)
				stats_data->vel_str_func[p - 1][i] += pow(cabs(run_data->u[n]), p);
				#if defined(__MAGNETO)
				stats_data->mag_str_func[p - 1][i] += pow(cabs(run_data->b[n]), p);
	    		#endif
	    		#endif

				// Compute the moments of the fluxes
				#if defined(__STR_FUNC_VEL_FLUX) || defined(__STR_FUNC_MAG_FLUX)
				stats_data->vel_flux_str_func[0][p - 1][i] += pow(pow(sgn(vel_enrg_flux_term), p) * fabs(vel_enrg_flux_term), p / 3.0);
				stats_data->vel_flux_str_func[1][p - 1][i] += pow(pow(sgn(vel_enrg_flux_term), p) * fabs(vel_hel_flux_term), p / 3.0);
				stats_data->vel_flux_str_func_abs[0][p - 1][i] += pow(fabs(vel_enrg_flux_term), p / 3.0);
				stats_data->vel_flux_str_func_abs[1][p - 1][i] += pow(fabs(vel_hel_flux_term), p / 3.0);
	    		#if defined(__MAGNETO)
				// stats_data->mag_flux_str_func[0][p - 1][i] += pow(mag_enrg_flux_term, p);
				// stats_data->mag_flux_str_func[1][p - 1][i] += pow(mag_hel_flux_term, p);
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
		stats_data->vel_str_func[i]      	= (double* )fftw_malloc(sizeof(double) * N);
		if (stats_data->vel_str_func[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Velocity Structure Function");
			exit(1);
		}
		#endif
		#if defined(__STR_FUNC_VEL_FLUX)
		for (int j = 0; j < 2; ++j) {
			stats_data->vel_flux_str_func[j][i] = (double* )fftw_malloc(sizeof(double) * N);
			if (stats_data->vel_flux_str_func[j][i] == NULL) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Velocity Flux Structure Function");
				exit(1);
			}
			stats_data->vel_flux_str_func_abs[j][i] = (double* )fftw_malloc(sizeof(double) * N);
			if (stats_data->vel_flux_str_func_abs[j][i] == NULL) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Absolute Velocity Flux Structure Function");
				exit(1);
			}		
		}
		#endif
		#if defined(__MAGNETO)
		#if defined(__STR_FUNC_MAG)
		stats_data->mag_str_func[i]         = (double* )fftw_malloc(sizeof(double) * N);
		if (stats_data->mag_str_func[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Magnetic Structure Function");
			exit(1);
		}
		#endif
		#if defined(__STR_FUNC_MAG_FLUX)
		for (int j = 0; j < 2; ++j) {
			stats_data->mag_flux_str_func[j][i] = (double* )fftw_malloc(sizeof(double) * N);
			if (stats_data->mag_flux_str_func[j][i] == NULL) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Magnetic Flux Sturcure Function");
				exit(1);
			}
			stats_data->mag_flux_str_func_abs[j][i] = (double* )fftw_malloc(sizeof(double) * N);
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
	stats_data->vel_moments = (double** )fftw_malloc(sizeof(double* ) * N);
	if (stats_data->vel_moments == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Velocity Stats");
		exit(1);
	}
	#if defined(__MAGNETO)
	stats_data->mag_moments = (double** )fftw_malloc(sizeof(double* ) * N);
	if (stats_data->mag_moments == NULL) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Magnetic Stats");
		exit(1);
	}
	#endif

	// Initialize the running stats objects
	for (int i = 0; i < N; ++i) {
		stats_data->vel_moments[i] = gsl_rstat_alloc();
		if (stats_data->vel_moments[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Running stats object: ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Velocity Stats");
			exit(1);
		}
		#if defined(__MAGNETO)
		stats_data->mag_moments[i] = gsl_rstat_alloc();
		if (stats_data->mag_moments[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Running stats object: ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Magnetic Stats");
			exit(1);
		}
		#endif
	}



	///--------------------------------- Histograms
	#if defined(__VEL_HIST)
	// Allocate memory for the array of histogram structs
	stats_data->real_vel_hist = (double** )fftw_malloc(sizeof(double* ) * N);
	for (int i = 0; i < N; ++i) {
		// Allocate	the histogram structs
		stats_data->real_vel_hist[i] = gsl_histogram_alloc(VEL_NUM_BINS);
		if (stats_data->real_vel_hist[i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for histogram struct for: ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Real Velocity Histogram");
			exit(1);
		}
	}
	#endif
}
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------