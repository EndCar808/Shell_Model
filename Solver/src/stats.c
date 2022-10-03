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
void ComputeStats(void) {

	// Initialize variables
	int n;
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
    	gsl_rstat_add(creal(run_data->u[n]), stats_data->vel_moments[i]);
    	#if defined(__MAGNETO)
    	gsl_rstat_add(creal(run_data->b[n]), stats_data->mag_moments[i]);
    	#endif
    	
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
    	for (int p = 2; p < NUM_POW; ++p) {
    		
    		// Compute the moments of the fields
    		#if defined(__STR_FUNC_MAG) || defined(__STR_FUNC_VEL)
			stats_data->vel_str_func[p - 2][i] += pow(cabs(run_data->u[n]), p);
			#if defined(__MAGNETO)
			stats_data->mag_str_func[p - 2][i] += pow(cabs(run_data->b[n]), p);
    		#endif
    		#endif

			// Compute the moments of the fluxes
			#if defined(__STR_FUNC_VEL_FLUX) || defined(__STR_FUNC_MAG_FLUX)
			stats_data->vel_flux_str_func[0][p - 2][i] += pow(vel_enrg_flux_term, p / 3.0);
			stats_data->vel_flux_str_func[1][p - 2][i] += pow(vel_hel_flux_term, p / 3.0);
    		#if defined(__MAGNETO)
			// stats_data->mag_flux_str_func[0][p - 2][i] += pow(mag_enrg_flux_term, p);
			// stats_data->mag_flux_str_func[1][p - 2][i] += pow(mag_hel_flux_term, p);
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
	for (int i = 0; i < NUM_POW - 2; ++i) {
		#if defined(__STR_FUNC_VEL)
		stats_data->vel_str_func[i]      	= (double* )fftw_malloc(sizeof(double) * N);
		#endif
		#if defined(__STR_FUNC_VEL_FLUX)
		stats_data->vel_flux_str_func[0][i] = (double* )fftw_malloc(sizeof(double) * N);
		stats_data->vel_flux_str_func[1][i] = (double* )fftw_malloc(sizeof(double) * N);
		#endif
		#if defined(__MAGNETO)
		#if defined(__STR_FUNC_MAG)
		stats_data->mag_str_func[i]         = (double* )fftw_malloc(sizeof(double) * N);
		#endif
		#if defined(__STR_FUNC_MAG_FLUX)
		stats_data->mag_flux_str_func[0][i] = (double* )fftw_malloc(sizeof(double) * N);
		stats_data->mag_flux_str_func[1][i] = (double* )fftw_malloc(sizeof(double) * N);
		#endif
		#endif
	}

	// Allocate memory for the running stats objects
	stats_data->vel_moments = (double** )fftw_malloc(sizeof(double* ) * N);
	#if defined(__MAGNETO)
	stats_data->mag_moments = (double** )fftw_malloc(sizeof(double* ) * N);
	#endif

	// Initialize the running stats objects
	for (int i = 0; i < N; ++i) {
		stats_data->vel_moments[i] = gsl_rstat_alloc();
		#if defined(__MAGNETO)
		stats_data->mag_moments[i] = gsl_rstat_alloc();
		#endif
	}
}
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------