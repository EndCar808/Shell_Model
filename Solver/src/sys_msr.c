/**
* @file utils.c  
* @author Enda Carroll
* @date Jun 2021
* @brief File containing the system measurables functions for the pseudospectral solver
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
#include "sys_msr.h"
#include "solver.h"
#include "utils.h"
// ---------------------------------------------------------------------
//  Function Definitions
// ---------------------------------------------------------------------
/**
 * Function to compute the system measurables such as energy, enstrophy, palinstrophy, helicity, energy and enstrophy dissipation rates, and spectra at once on the local processes for the current timestep
 * @param t 		The current time in the simulation
 * @param iter 		The index in the system arrays for the current timestep
 * @param RK_data 	Struct containing the integration varaiables needed for the nonlinear term function
 */
void ComputeSystemMeasurables(double t, int iter, RK_data_struct* RK_data) {

	// Initialize variables
    const long int N = sys_vars->N; 
    int n, l;
    double k_fac;
    const double lambda_pow         = sys_vars->Lambda * sys_vars->Lambda;
    const double interact_coeff_u_1 = sys_vars->EPS / sys_vars->Lambda;
    const double interact_coeff_u_2 = (1.0 - sys_vars->EPS) / lambda_pow;
    const double interact_coeff_b_1 = 1.0 - sys_vars->EPS - sys_vars->EPS_M;
    const double interact_coeff_b_2 = sys_vars->EPS_M / sys_vars->Lambda;
    const double interact_coeff_b_3 = (1.0 - sys_vars->EPS_M) / lambda_pow;
    double k_pre_fac_1, k_pre_fac_2, k_pre_fac_3, k_pre_fac_4;

    // Record the initial time
    #if defined(__TIME)
    if (sys_vars->TRANS_ITERS_FLAG != TRANSIENT_ITERS) {
        run_data->time[iter] = t;
    }
    #endif

    // If adaptive stepping check if within memory limits
    if ((iter >= sys_vars->num_print_steps) && (iter % 100 == 0)) {
        // Print warning to screen if we have exceeded the memory limits for the system measurables arrays
        printf("\n["MAGENTA"WARNING"RESET"] --- Unable to write system measures at Indx: [%d] t: [%lf] ---- Number of intergration steps is now greater then memory allocated\n", iter, t);
    }

    // ------------------------------------
    // Initialize Measurables
    // ------------------------------------
    #if defined(__SYS_MEASURES)
    if (iter < sys_vars->num_print_steps) {
        // Initialize totals
        run_data->tot_energy[iter]    = 0.0;
        run_data->tot_hel_u[iter]     = 0.0;
        #if defined(__MAGNETO)
        run_data->tot_hel_b[iter]     = 0.0;
        run_data->tot_cross_hel[iter] = 0.0;
        #endif
    }
    #endif
    #if defined(__ENRG_FLUX)
    for (int i = 0; i < N; ++i) {
        // Initialize the energy dissipation
        run_data->energy_flux[i]  = 0.0;
        run_data->energy_diss[i]  = 0.0;
        run_data->energy_input[i] = 0.0;
    }
    #endif

    // ------------------------------------
    // Remove Forcing
    // ------------------------------------
    // Remove forcing to compute the energy flux balance terms
    for (int i = 0; i < (N + 4); ++i) {
        run_data->u[n] -= run_data->forcing_u[n];
        #if defined(__MAGNETO)
        run_data->b[n] -= run_data->forcing_b[n];
        #endif
    }
    
    // -------------------------------------
    // Compute Measurables
    // -------------------------------------
    for (int i = 0; i < N; ++i) {
        // Get tmp index
        n = i + 2;

        if (iter < sys_vars->num_print_steps) {
            //-------------- System Measures
            #if defined(__SYS_MEASURES)
            // Compute Velocity helicity k factors
            if (-log_lambda(fabs(sys_vars->EPS - 1.0) / 2.0) < 0) { 
                k_fac = 1.0 / pow(run_data->k[i], log_lambda(fabs(sys_vars->EPS - 1.0) / 2.0));
            } 
            else {
              k_fac = pow(run_data->k[i], -log_lambda(fabs(sys_vars->EPS - 1.0) / 2.0));  
            }
            // Update sum for totals
            #if defined(PHASE_ONLY_DIRECT)
            run_data->tot_energy[iter]    += run_data->a_n[n] * run_data->a_n[n];
            run_data->tot_hel_u[iter]     += pow(sgn(sys_vars->EPS - 1.0), i) * (run_data->a_n[n] * run_data->a_n[n]) * k_fac;
            #else
            run_data->tot_energy[iter]    += cabs(run_data->u[n] * conj(run_data->u[n]));
            run_data->tot_hel_u[iter]     += pow(sgn(sys_vars->EPS - 1.0), i) * cabs(run_data->u[n] * conj(run_data->u[n])) * k_fac;
            #endif
            #if defined(__MAGNETO)
            #if defined(PHASE_ONLY_DIRECT)
            run_data->tot_energy[iter]    += run_data->b_n[n] * run_data->b_n[n];
            run_data->tot_hel_b[iter]     += pow(sgn(sys_vars->EPS - 1.0), i) * (run_data->b_n[n] * run_data->b_n[n]) / run_data->k[i];
            run_data->tot_cross_hel[iter] += creal(run_data->a_n[n] * run_data->b_n[n]);
            #else
            run_data->tot_energy[iter]    += cabs(run_data->b[n] * conj(run_data->b[n]));
            run_data->tot_hel_b[iter]     += pow(sgn(sys_vars->EPS - 1.0), i) * cabs(run_data->b[n] * conj(run_data->b[n])) / run_data->k[i];
            run_data->tot_cross_hel[iter] += creal(run_data->u[n] * conj(run_data->b[n]));
            #endif
            #endif
            #endif

            //-------------- Energy Flux and Dissipation
            #if defined(__ENRG_FLUX)
            // Compute the energy dissipation
            for (int j = 0; j < i; ++j) {
                // Get temp indx
                l = j + 2;
                #if defined(PHASE_ONLY_DIRECT)
                run_data->energy_diss[i]  += run_data->k[j] * run_data->k[j] * sys_vars->NU * run_data->a_n[n] * run_data->a_n[n];
                run_data->energy_force[i] += run_data->a_n[l] * cabs(run_data->forcing_u[l]) * cos(run_data->phi_n[l] - carg(run_data->forcing_u[l]));
                #if defined(__MAGNETO) 
                run_data->energy_diss[i]  += run_data->k[j] * run_data->k[j] * sys_vars->ETA * run_data->b_n[n] * run_data->b_n[n];
                run_data->energy_force[i] += run_data->b_n[l] * cabs(run_data->forcing_b[l]) * cos(run_data->psi_n[l] - carg(run_data->forcing_b[l]));
                #endif
                #else
                run_data->energy_diss[i]  += run_data->k[j] * run_data->k[j] * sys_vars->NU * cabs(run_data->u[l] * conj(run_data->u[l]));
                run_data->energy_force[i] += creal(run_data->u[l] * conj(run_data->forcing_u[l]));
                #if defined(__MAGNETO) 
                run_data->energy_diss[i]  += run_data->k[j] * run_data->k[j] * sys_vars->ETA * cabs(run_data->b[l] * conj(run_data->b[l])); 
                run_data->energy_force[i] += creal(run_data->b[l] * conj(run_data->forcing_b[l]));
                #endif
                #endif                  
            }

            // Get the correct k prefactor terms for the nonlinear flux term 
            if (i == 0) {
                k_pre_fac_1 =  - run_data->k[i] * interact_coeff_u_1;
                #if defined(__MAGNETO)
                k_pre_fac_2 = run_data->k[i] * interact_coeff_b_2;
                k_pre_fac_3 = run_data->k[i] * interact_coeff_u_1;
                k_pre_fac_4 = run_data->k[i] * interact_coeff_b_2;
                #endif
            }
            else {
                k_pre_fac_1 = run_data->k[i - 1] - run_data->k[i] * interact_coeff_u_1;
                #if defined(__MAGNETO)
                k_pre_fac_2 = - run_data->k[i - 1] + run_data->k[i] * interact_coeff_b_2;
                k_pre_fac_3 = run_data->k[i - 1] * interact_coeff_b_1 + run_data->k[i] * interact_coeff_u_1;
                k_pre_fac_4 = run_data->k[i - 1] * interact_coeff_b_1 + run_data->k[i] * interact_coeff_b_2;
                #endif
            }

            // Compute the energy flux
            #if defined(PHASE_ONLY_DIRECT)
            // First term
            run_data->energy_flux[i] = k_pre_fac_1 * run_data->a_n[n - 1] * run_data->a_n[n] * run_data->a_n[n + 1] * sin(run_data->phi_n[n - 1] + run_data->phi_n[n] + run_data->phi_n[n + 1]);
            run_data->energy_flux[i] += run_data->k[i] * run_data->a_n[n] * run_data->a_n[n + 1] * run_data->a_n[n + 2] * sin(run_data->phi_n[n] + run_data->phi_n[n + 1] + run_data->phi_n[n + 2]);
            #if defined(__MAGNETO)
            // Second term
            run_data->energy_flux[i] += k_pre_fac_2 * run_data->a_n[n - 1] * run_data->b_n[n] * run_data->b_n[n + 1] * sin(run_data->phi_n[n - 1] + run_data->psi_n[n] + run_data->psi_n[n + 1]);
            run_data->energy_flux[i] += -run_data->k[i] * run_data->a_n[n] * run_data->b_n[n + 1] * run_data->b_n[n + 2] * sin(run_data->phi_n[n] + run_data->psi_n[n + 1] + run_data->psi_n[n + 2]);
            // Third term
            run_data->energy_flux[i] += k_pre_fac_3 * run_data->b_n[n - 1] * run_data->a_n[n] * run_data->b_n[n + 1] * sin(run_data->psi_n[n - 1] + run_data->phi_n[n] + run_data->psi_n[n + 1]);
            run_data->energy_flux[i] += interact_coeff_b_1 * run_data->k[i] * run_data->b_n[n] * run_data->a_n[n + 1] * run_data->b_n[n + 2] * sin(run_data->psi_n[n] + run_data->phi_n[n + 1] + run_data->psi_n[n + 2]);
            // Fourth term
            run_data->energy_flux[i] += -k_pre_fac_4 * run_data->b_n[n - 1] * run_data->a_n[n] * run_data->b_n[n + 1] * sin(run_data->psi_n[n - 1] + run_data->psi_n[n] + run_data->psi_n[n + 1]);
            run_data->energy_flux[i] += - interact_coeff_b_1 * run_data->k[i] * run_data->b_n[n] * run_data->a_n[n + 1] * run_data->b_n[n + 2] * sin(run_data->psi_n[n] + run_data->psi_n[n + 1] + run_data->psi_n[n + 2]);
            #endif
            #else
            // First term
            run_data->energy_flux[i] = cimag(k_pre_fac_1 * run_data->u[n - 1] * run_data->u[n] * run_data->u[n + 1] + run_data->k[i] * run_data->u[n] * run_data->u[n + 1] * run_data->u[n + 2]);
            #if defined(__MAGNETO)
            // Second term
            run_data->energy_flux[i] += cimag(k_pre_fac_2 * run_data->u[n - 1] * run_data->b[n] * run_data->b[n + 1] - run_data->k[i] * run_data->u[n] * run_data->b[n + 1] * run_data->b[n + 2]);
            // Third term
            run_data->energy_flux[i] += cimag(k_pre_fac_3 * run_data->b[n - 1] * run_data->u[n] * run_data->b[n + 1] + interact_coeff_b_1 * run_data->k[i] * run_data->b[n] * run_data->u[n + 1] * run_data->b[n + 2]);
            // Fourth term
            run_data->energy_flux[i] += cimag(-k_pre_fac_4 * run_data->b[n - 1] * run_data->b[n] * run_data->u[n + 1] - interact_coeff_b_1 * run_data->k[i] * run_data->b[n] * run_data->b[n + 1] * run_data->u[n + 2]);
            #endif
            #endif
            #endif
        }
    }    

    // ------------------------------------
    // Normalize Measureables 
    // ------------------------------------ 
    #if defined(__SYS_MEASURES)
    if (iter < sys_vars->num_print_steps) {
        run_data->tot_energy[iter]    *= 0.5;
        run_data->tot_hel_u[iter]     *= 0.5;
        #if defined(__MAGNETO)
        run_data->tot_hel_b[iter]     *= 0.5;
        run_data->tot_cross_hel[iter] *= 0.5;
        #endif
    }
    #endif

    // ------------------------------------
    // Add Back Forcing
    // ------------------------------------
    // Add back forcing after computing the energy flux balance terms
    for (int i = 0; i < (N + 4); ++i) {
        run_data->u[n] += run_data->forcing_u[n];
        #if defined(__MAGNETO)
        run_data->b[n] += run_data->forcing_b[n];
        #endif
    }
}
/**
 * Function to initialize and compute the system measurables and spectra of the initial conditions
 * @param RK_data The struct containing the Runge-Kutta arrays to compute the nonlinear term for the fluxes
 */
void InitializeSystemMeasurables(RK_data_struct* RK_data) {

    // Set the size of the arrays to twice the number of printing steps to account for extra steps due to adaptive stepping
    if (sys_vars->ADAPT_STEP_FLAG == ADAPTIVE_STEP) {
        sys_vars->num_print_steps = 2 * sys_vars->num_print_steps;
    }
    else {
        sys_vars->num_print_steps = sys_vars->num_print_steps;
    }
    int print_steps = sys_vars->num_print_steps;

    // --------------------------------
    // Allocate System Totals Memory
    // --------------------------------
    #if defined(__SYS_MEASURES)
    // Total Energy in the system
    run_data->tot_energy = (double* )fftw_malloc(sizeof(double) * print_steps);
    if (run_data->tot_energy == NULL) {
        fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Total Energy");
        exit(1);
    }
    #if defined(__MAGNETO)
    // Total Velocity Helicity
    run_data->tot_hel_u = (double* )fftw_malloc(sizeof(double) * print_steps);
    if (run_data->tot_hel_u == NULL) {
        fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Total Velocity Helicity");
        exit(1);
    }   
    // Total Magnetic Helicity
    run_data->tot_hel_b = (double* )fftw_malloc(sizeof(double) * print_steps);
    if (run_data->tot_hel_b == NULL) {
        fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Total Magnetic Helicity");
        exit(1);
    }   

    // Total Enstrophy
    run_data->tot_cross_hel = (double* )fftw_malloc(sizeof(double) * print_steps);
    if (run_data->tot_cross_hel == NULL) {
        fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Total Cross Helicity");
        exit(1);
    }   
    #endif
    #endif

    // Time
    #if defined(__TIME)
    run_data->time = (double* )fftw_malloc(sizeof(double) * print_steps);
    if (run_data->time == NULL) {
        fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Time");
        exit(1);
    }
    #endif

    #if defined(__ENRG_FLUX)
    // Allocate energy flux
    run_data->energy_flux = (double* )fftw_malloc(sizeof(double) * sys_vars->N);
    if (run_data->energy_flux == NULL) {
        fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Energy Flux");
        exit(1);
    }
    // Allocate energy dissipation
    run_data->energy_diss = (double* )fftw_malloc(sizeof(double) * sys_vars->N);
    if (run_data->energy_diss == NULL) {
        fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Energy Dissipation");
        exit(1);
    }
    // Allocate energy input
    run_data->energy_force = (double* )fftw_malloc(sizeof(double) * sys_vars->N);
    if (run_data->energy_force == NULL) {
        fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for the ["CYAN"%s"RESET"]\n-->> Exiting!!!\n", "Energy Input");
        exit(1);
    }
    #endif

    // ----------------------------
    // Get Measurables of the ICs
    // ----------------------------
    if (sys_vars->TRANS_ITERS_FLAG != TRANSIENT_ITERS) {
        ComputeSystemMeasurables(0.0, 0, RK_data);
    }
}
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------