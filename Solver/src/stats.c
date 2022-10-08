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
// ---------------------------------------------------------------------
//  Function Definitions
// ---------------------------------------------------------------------
/**
 * Function to compute the statistics for the current iteration
 * @param iters          The current iteration of the simulation
 * @param save_data_indx The current index for saving data to file
 */
void ComputeStats(const int iters, const int save_data_indx) {

	// Initialize variables
	int n;
	int gsl_status;
	const long int N = sys_vars->N;
	double vel_enrg_flux_term, vel_hel_flux_term;
	#if defined(__MAGNETO)
	double mag_enrg_flux_term, mag_hel_flux_term;
	#endif

	// ------------------------------------
    // Update Stats Counter
    // ------------------------------------
	stats_data->num_stats_steps++;

	// ------------------------------------
    // Compute Stats
    // ------------------------------------
    for (int i = 0; i < N; ++i) {
    	n = i + 2;

    	// Update the running sums for the field stats
    	gsl_status = gsl_rstat_add(creal(run_data->u[n]), stats_data->vel_moments[i]);
    	if (gsl_status != 0) {
    		// Print error message to error stream
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to running stats counter ["CYAN"%s"RESET"] - Error Value: ["CYAN"%d"RESET"] Field Value: ["CYAN"%lf"RESET"] \n-->> Exiting!!!\n", "Velocity Stats", gsl_status, creal(run_data->u[n]));

			// Save current state of system to file before exiting
			FinalWriteAndCloseOutputFile(N, iters, save_data_indx);
			
			// Exit programme
			exit(1);
		}
    	#if defined(__MAGNETO)
    	gsl_status = gsl_rstat_add(creal(run_data->b[n]), stats_data->mag_moments[i]);
	    if (gsl_status != 0) {
	    	// Print error message to error stream
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to running stats counter ["CYAN"%s"RESET"] - Error Value: ["CYAN"%d"RESET"] Field Value: ["CYAN"%lf"RESET"] \n-->> Exiting!!!\n", "Magnetic Stats", gsl_status, creal(run_data->b[n]));

			// Save current state of system to file before exiting
			FinalWriteAndCloseOutputFile(N, iters, save_data_indx);
			
			// Exit programme
			exit(1);
		}
    	#endif
    	
    	// Compute the flux terms
    	#if defined(__STR_FUNC_VEL_FLUX) || defined(__STR_FUNC_MAG_FLUX)
		vel_enrg_flux_term = fabs(cimag(run_data->u[n + 2] * run_data->u[n + 1] * run_data->u[n] + (1.0 - sys_vars->EPS) / sys_vars->Lambda * run_data->u[n + 1] * run_data->u[n] * run_data->u[n - 1]));
		vel_hel_flux_term  = fabs(cimag(run_data->u[n + 2] * run_data->u[n + 1] * run_data->u[n] - (sys_vars->EPS * sys_vars->Lambda + 1.0) / pow(sys_vars->Lambda, 2.0) * run_data->u[n + 1] * run_data->u[n] * run_data->u[n - 1]));
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
			stats_data->vel_flux_str_func[0][p - 1][i] += pow(vel_enrg_flux_term, p / 3.0);
			stats_data->vel_flux_str_func[1][p - 1][i] += pow(vel_hel_flux_term, p / 3.0);
    		#if defined(__MAGNETO)
			// stats_data->mag_flux_str_func[0][p - 1][i] += pow(mag_enrg_flux_term, p);
			// stats_data->mag_flux_str_func[1][p - 1][i] += pow(mag_hel_flux_term, p);
    		#endif
    		#endif
    	}
    }
}
/**
 * Function to allocate and initialize the stats objects
 */
void InitializeStats(void) {

	// Initialize variables
	const long int N = sys_vars->N;

	// Initialize the stats counter
	stats_data->num_stats_steps = 0;


	// ------------------------------------
    // Allocate & Initialize Stats Objects
    // ------------------------------------
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
		stats_data->vel_flux_str_func[0][i] = (double* )fftw_malloc(sizeof(double) * N);
		if (stats_data->vel_flux_str_func[0][i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Velocity Flux Structure Function");
			exit(1);
		}
		stats_data->vel_flux_str_func[1][i] = (double* )fftw_malloc(sizeof(double) * N);
		if (stats_data->vel_flux_str_func[1][i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Velocity Flux Structure Function");
			exit(1);
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
		stats_data->mag_flux_str_func[0][i] = (double* )fftw_malloc(sizeof(double) * N);
		if (stats_data->mag_flux_str_func[0][i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Magnetic Flux Sturcure Function");
			exit(1);
		}
		stats_data->mag_flux_str_func[1][i] = (double* )fftw_malloc(sizeof(double) * N);
		if (stats_data->mag_flux_str_func[1][i] == NULL) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for Stats Array ["CYAN"%s"RESET"] \n-->> Exiting!!!\n", "Magnetic Flux Sturcure Function");
			exit(1);
		}
		#endif
		#endif
	}

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
}
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------