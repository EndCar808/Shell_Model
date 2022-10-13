/**
* @file utils.c  
* @author Enda Carroll
* @date Sept 2022
* @brief File containing the utilities functions for the solver
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
 * Function to read in arguements given at the command line upon execution of the solver
 * @param  argc The number of arguments given
 * @param  argv Array containg the arguments specified
 * @return      Returns 0 if the parsing of arguments has been successful
 */
int GetCMLArgs(int argc, char** argv) {

	// Initialize Variables
	int c;
	int force_flag      = 0;
	int cfl_flag        = 0;
	int interact_flag   = 0;
	int shell_wave_flag = 0;
	int visc_flag       = 0;
	int mag_diff_flag   = 0;
	int output_dir_flag = 0;
	int trans_iter_flag = 0;
	int time_step_flag  = 0;

	// -------------------------------
	// Initialize Default Values
	// -------------------------------
	// Input file path
	strncpy(file_info->input_file_name, "NONE", 512);
	strncpy(file_info->input_dir, "./Data/Solver/InitialConditions/", 512);  // Set default to Initial Conditions folder
	// Output file directory
	strncpy(file_info->output_dir, "./Data/Tmp/", 512);  // Set default output directory to the Tmp folder
	strncpy(file_info->output_tag, "No-Tag", 64);
	file_info->file_only = 0; // used to indicate if output file should be file only i.e., not output folder
	// System dimensions
	sys_vars->N = 19;
	// Integration time 
	sys_vars->t0              = 0.0;
	sys_vars->dt              = 1e-3;
	sys_vars->T               = 12.0;
	sys_vars->CFL_CONST       = 0.9;
	sys_vars->ADAPT_STEP_FLAG = 0;
	sys_vars->CFL_COND_FLAG   = 1;
	// Initial conditions
	strncpy(sys_vars->u0, "N_SCALING", 64);
	// Energy Spectra scaling
	sys_vars->ALPHA = 1.5;
	sys_vars->BETA  = 1.5;
	// Shell wavenumber coefficients
	sys_vars->k_0    = K_0;
	sys_vars->Lambda = LAMBDA;	
	// Interaction coefficients
	sys_vars->EPS   = EPSILON;
	sys_vars->EPS_M = EPSILON_M;	
	// Forcing
	strncpy(sys_vars->forcing, "NONE", 64);	
	sys_vars->force_k         = 0;
	sys_vars->force_scale_var = 1.0;
	// Transient iterations
	sys_vars->TRANS_ITERS_FLAG = 0;
	sys_vars->TRANS_ITERS_FRAC = TRANSIENT_FRAC;
	// System viscosity / diffusivity
	sys_vars->NU  = 1e-7;
	sys_vars->HYPER_VISC_FLAG = 0;
	sys_vars->HYPER_VISC_POW  = VISC_POW;
	// // System ekman dray / hypodiffusivity
	sys_vars->ETA = 1e-7;
	sys_vars->HYPO_MAG_DIFF_FLAG = 0;
	sys_vars->HYPO_MAG_DIFF_POW  = HYPO_DIFF_POW;
	// Write to file every 
	sys_vars->SAVE_EVERY = (long int)100;

	// -------------------------------
	// Parse CML Arguments
	// -------------------------------
	while ((c = getopt(argc, argv, "o:h:n:a:b:d:s:e:t:v:i:c:p:f:z:y:w:T:")) != -1) {
		switch(c) {
			case 'o':
				if (output_dir_flag == 0) {
					// Read in location of output directory
					strncpy(file_info->output_dir, optarg, 512);
					output_dir_flag++;
				}
				else if (output_dir_flag == 1) {
					// Output file only indicated
					file_info->file_only = 1;
				}
				break;
			case 'n':
				// Read in the dimensions of the system
				sys_vars->N = atoi(optarg);

				// Check dimensions satisfy requirements
				if (sys_vars->N <= 2) {
					fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The minimum dimension size of [%ld] must be greater than 2\n-->> Exiting!\n\n", sys_vars->N);		
					exit(1);
				}
				break;
			case 's':
				// Read in intial time
				sys_vars->t0 = atof(optarg);
				if (sys_vars->t0 < 0) {
				    fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The integration start time [%lf] must be a positive\n-->> Exiting!\n", sys_vars->t0);      
				    exit(1);
				}	
				break;
			case 'e':
				// Read in final time
				sys_vars->T = atof(optarg);	
				if (sys_vars->T < 0) {
				    fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The integration end time [%lf] must be a positive\n-->> Exiting!\n", sys_vars->T);      
				    exit(1);
				}
				else if (sys_vars->T < sys_vars->t0) {
					fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided end time: [%lf] must be greater than the initial time: [%lf]\n-->> Exiting!\n\n", sys_vars->T, sys_vars->t0);		
					exit(1);
				}
				break;
			case 'h':
				if (time_step_flag == 0) {
					// Read in initial timestep
					sys_vars->dt = atof(optarg);

					if (sys_vars->dt <= 0) {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided timestep: [%lf] must be strictly positive\n-->> Exiting!\n\n", sys_vars->dt);		
						exit(1);
					}
					time_step_flag = 1;
				}
				else if (time_step_flag == 1) {
					// Read in adaptive step indicator
					sys_vars->ADAPT_STEP_FLAG = atoi(optarg);
					if ((sys_vars->ADAPT_STEP_FLAG == 0) || (sys_vars->ADAPT_STEP_FLAG == 1)) {
					}
					else {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: Incorrect Adaptive stepping flag: [%d] Must be either 0 or 1 -- Set to 0 (no adaptive stepping) by default\n-->> Exiting!\n\n", sys_vars->ADAPT_STEP_FLAG);		
						exit(1);
					}
					// time_step_flag = 2;
					break;
				}
				break;			
			case 'c':
				if (cfl_flag == 0) {
					// Read in value of the CFL -> this can be used to control the timestep
					sys_vars->CFL_COND_FLAG = atoi(optarg);
					if ((sys_vars->CFL_COND_FLAG == 0) || (sys_vars->CFL_COND_FLAG == 1)) {
					}
					else {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: Incorrect CFL condition flag: [%d] Must be either 0 or 1 -- Set to 1 (use CFL condition) by default\n-->> Exiting!\n\n", sys_vars->CFL_COND_FLAG);		
						exit(1);
					}
					cfl_flag = 1;
					break;	
				}
				else if (cfl_flag == 1) {
					sys_vars->CFL_CONST = atof(optarg);
					if (sys_vars->CFL_CONST <= 0) {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided CFL Constant: [%lf] must be strictly positive\n-->> Exiting!\n\n", sys_vars->CFL_CONST);		
						exit(1);
					}
					break;
				}
				break;
			case 'a':
				// Read in scalining exponent of the velocity energy spectra
				sys_vars->ALPHA = atof(optarg);

				if (sys_vars->ALPHA < 0) {
				    fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The velocity energy spectra slope [%lf] must be a positive\n-->> Exiting!\n", sys_vars->ALPHA);      
				    exit(1);
				}	
				break;
			case 'b':
				// Read in the scalining exponent of the velocity energy spectra
				sys_vars->BETA = atof(optarg);

				if (sys_vars->BETA < 0) {
				    fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The magnetic energy spectra slope [%lf] must be a positive\n-->> Exiting!\n", sys_vars->BETA);      
				    exit(1);
				}	
				break;				
			case 'v':
				if (visc_flag == 0) {
					// Read in the viscosity
					sys_vars->NU = atof(optarg);
					if (sys_vars->NU < 0) {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided viscosity: [%lf] must be positive\n-->> Exiting!\n\n", sys_vars->NU);		
						exit(1);
					}
					visc_flag = 1;
					break;
				}
				else if (visc_flag == 1) {
					// Read in hyperviscosity flag
					sys_vars->HYPER_VISC_FLAG = atoi(optarg);
					if ((sys_vars->HYPER_VISC_FLAG == 0) || (sys_vars->HYPER_VISC_FLAG == 1)) {
					}
					else {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: Incorrect Hyperviscosity flag: [%d] Must be either 0 or 1 -- Set to 0 (no Hyperviscosity) by default\n-->> Exiting!\n\n", sys_vars->HYPER_VISC_FLAG);		
						exit(1);
					}
					visc_flag = 2;
					break;	
				}
				else if (visc_flag == 2) {
					// Read in the hyperviscosity power
					sys_vars->HYPER_VISC_POW = atof(optarg);
					if (sys_vars->HYPER_VISC_POW <= 0.0) {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided hyperviscosity power: [%lf] must be strictly positive\n-->> Exiting!\n\n", sys_vars->HYPER_VISC_POW);		
						exit(1);
					}
					break;
				}
				break;
			case 'd':
				if (mag_diff_flag == 0) {
					// Read in the magnetic diffusivity
					sys_vars->ETA = atof(optarg);
					if (sys_vars->ETA < 0) {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided Magnetic Diffusivity: [%lf] must be positive\n-->> Exiting!\n\n", sys_vars->ETA);		
						exit(1);
					}
					mag_diff_flag = 1;
					break;
				}
				else if (mag_diff_flag == 1) {
					// Read in hypodiffusivity flag
					sys_vars->HYPO_MAG_DIFF_FLAG = atoi(optarg);
					if ((sys_vars->HYPO_MAG_DIFF_FLAG == 0) || (sys_vars->HYPO_MAG_DIFF_FLAG == 1)) {
					}
					else {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: Incorrect Hypodiffusivity flag: [%d] Must be either 0 or 1 -- Set to 0 (no drag) by default\n-->> Exiting!\n\n", sys_vars->HYPO_MAG_DIFF_FLAG);		
						exit(1);
					}
					mag_diff_flag = 2;
					break;	
				}
				else if (mag_diff_flag == 2) {
					// Read in the hypodiffusivity power
					sys_vars->HYPO_MAG_DIFF_POW = atof(optarg);
					if (sys_vars->HYPO_MAG_DIFF_POW >= 0.0) {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided Hypodiffusivity power: [%lf] must be strictly negative\n-->> Exiting!\n\n", sys_vars->HYPO_MAG_DIFF_POW);		
						exit(1);
					}
					break;
				}
				break;
			case 'w':
				if (shell_wave_flag == 0) {
					// Read in the shell wavenumber prefactor
					sys_vars->k_0 = atof(optarg);
					if (sys_vars->k_0 < 0) {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided shell wavenumber factor: [%lf] must be positive\n-->> Exiting!\n\n", sys_vars->k_0);		
						exit(1);
					}
					shell_wave_flag = 1;
					break;
				}
				else if (shell_wave_flag == 1) {
					// Read in the intershell ratio for the shell radius
					sys_vars->Lambda = atof(optarg);
					if (sys_vars->Lambda < 0) {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided intershell ratio lambda: [%lf] must be positive\n-->> Exiting!\n\n", sys_vars->Lambda);		
						exit(1);
					}
					break;	
				}
				break;
			case 'y':
				if (interact_flag == 0) {
					// Read in the interaction coefficient for the velocity equation
					sys_vars->EPS = atof(optarg);
					if (sys_vars->EPS < 0) {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The velocity interaction coefficient: [%lf] must be positive\n-->> Exiting!\n\n", sys_vars->EPS);		
						exit(1);
					}
					interact_flag = 1;
					break;
				}
				else if (interact_flag == 1) {
					// Read in the interaction coefficient for the magnetic equation
					sys_vars->EPS_M = atof(optarg);
					if (sys_vars->EPS_M < 0) {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The magnetic interaction coefficient: [%lf] must be positive\n-->> Exiting!\n\n", sys_vars->EPS_M);		
						exit(1);
					}
					break;	
				}
				break;
			case 'i':
				// Read in the initial conditions
				if (!(strcmp(optarg,"N_SCALING"))) {
					// Scaled Initial Conditions
					strncpy(sys_vars->u0, "N_SCALING", 64);
					break;
				}
				else if (!(strcmp(optarg,"RANDOM"))) {
					// Random Initial Conditions
					strncpy(sys_vars->u0, "RANDOM", 64);
					break;
				}
				else if (!(strcmp(optarg,"ZERO"))) {
					// Zero Phases Initial Conditions
					strncpy(sys_vars->u0, "ZERO", 64);
					break;
				}
				else {
					// No initial conditions specified -> this will default to random initial conditions
					strncpy(sys_vars->u0, "NONE", 64);
					break;
				}
				break;	
			case 't':
				// Read in output directory tag
				strncpy(file_info->output_tag, optarg, 64);	
				break;
			case 'T':
				if (trans_iter_flag == 0) {
					// Read in transient iterations flag
					sys_vars->TRANS_ITERS_FLAG = atoi(optarg);
					if ((sys_vars->TRANS_ITERS_FLAG == 0) || (sys_vars->TRANS_ITERS_FLAG == 1)) {
					}
					else {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: Incorrect transient iterations flag: [%d] Must be either 0 or 1 -- Set to 0 (no transient iterations) by default\n-->> Exiting!\n\n", sys_vars->TRANS_ITERS_FLAG);		
						exit(1);
					}
					trans_iter_flag = 1;
					break;	
				}
				else if (trans_iter_flag == 1) {
					// Read in transient iterations fraction of total iterations
					sys_vars->TRANS_ITERS_FRAC = atof(optarg);
					if ((sys_vars->TRANS_ITERS_FRAC >= 1.0) || (sys_vars->TRANS_ITERS_FRAC <= 0.0)) {
						fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: Incorrect transient iterations fraction: [%lf] Must be between 0 and 1 -- Set to 0.2 by default\n-->> Exiting!\n\n", sys_vars->TRANS_ITERS_FRAC);		
						// exit(1);
						sys_vars->TRANS_ITERS_FRAC = 0.2;
					}
					trans_iter_flag = 2;
					break;
				}
				break;				
			case 'p':
				// Read in how often to print to file
				sys_vars->SAVE_EVERY = (long int)atoi(optarg);
				break;
			case 'f':
				// Read in the forcing type
				if (!(strcmp(optarg,"DELTA")) && (force_flag == 0)) {
					// Delta function on mode 0
					strncpy(sys_vars->forcing, "DELTA", 64);
					force_flag = 1;
					break;
				}
				else if (!(strcmp(optarg,"NONE"))  && (force_flag == 0)) {
					// No forcing
					strncpy(sys_vars->forcing, "NONE", 64);
					force_flag = 1;
					break;
				}
				else if (!(strcmp(optarg,"STOC"))  && (force_flag == 0)) {
					// Stochastic forcing
					strncpy(sys_vars->forcing, "STOC", 64);
					force_flag = 1;
					break;
				}
				else if ((force_flag == 1)) {
					// Get the forcing wavenumber
					sys_vars->force_k = atoi(optarg);
					force_flag = 2;
					break;
				}
				else if ((force_flag == 2)) {
					// Get the force scaling variable
					sys_vars->force_scale_var = atof(optarg);
					break;
				}
				break;
			case 'z':
				strncpy((file_info->input_file_name), optarg, 512);	// copy the input file name given as a command line argument
				if ( access((file_info->input_file_name), F_OK) != 0) {
					fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The input file [%s] cannot be found, please ensure correct path to file is specified.\n", (file_info->input_file_name));		
					exit(1);					
				}
				break;
			default:
				fprintf(stderr, "\n["RED"ERROR"RESET"] Incorrect command line flag encountered\n");		
				fprintf(stderr, "Use"YELLOW" -o"RESET" to specify the output directory\n");
				fprintf(stderr, "Use"YELLOW" -n"RESET" to specify the size of each dimension in the system\n");
				fprintf(stderr, "Use"YELLOW" -s"RESET" to specify the start time of the simulation\n");
				fprintf(stderr, "Use"YELLOW" -e"RESET" to specify the end time of the simulation\n");
				fprintf(stderr, "Use"YELLOW" -h"RESET" to specify the timestep\n");
				fprintf(stderr, "Use"YELLOW" -c"RESET" to specify the CFL constant for the adaptive stepping\n");
				fprintf(stderr, "Use"YELLOW" -v"RESET" to specify the system viscosity\n");
				fprintf(stderr, "Use"YELLOW" -i"RESET" to specify the initial condition\n");
				fprintf(stderr, "Use"YELLOW" -t"RESET" to specify the tag name to be used in the output file directory\n");
				fprintf(stderr, "Use"YELLOW" -f"RESET" to specify the forcing type\n");
				fprintf(stderr, "Use"YELLOW" -T"RESET" to specify if transient iterations are to be performed\n");
				fprintf(stderr, "Use"YELLOW" -p"RESET" to specify how often to print to file\n");
				fprintf(stderr, "Use"YELLOW" -z"RESET" to specify an input file to read parameters from\n");
				fprintf(stderr, "\nExample usage:\n"CYAN"\t./bin/main -o \"../Data/Tmp\" -n 64 -n 64 -s 0.0 -e 1.0 -h 0.0001 -v 1.0 -i \"TG_VORT\" -t \"TEMP_RUN\" \n"RESET);
				fprintf(stderr, "-->> Now Exiting!\n\n");
				exit(1);
		}
	}

	return 0;
}
/**
 * Function that prints the summary details of the simulation to a .txt
 * @param sim_time The execution time of the simulation 
 * @param argc     The number of command line arguments
 * @param argv     Array conatining the command line arguments
 */
void PrintSimulationDetails(int argc, char** argv, double sim_time) {

	// Initialize variables
	FILE *sim_file;
	char sys_type[64];
	char solv_type[64];
	char model_type[64];
	char file_path[512];

	// -------------------------------
	// Open File
	// -------------------------------
	strcpy(file_path, file_info->output_dir);
	strcat(file_path, "SimulationDetails.txt");
	sim_file = fopen(file_path, "w");

	// -------------------------------
	// Print Executing Command
	// -------------------------------
	fprintf(sim_file, "Executing Command:"); 
	fprintf(sim_file, "\n\n\t");
	for (int i = 0; i < argc; ++i) {
		fprintf(sim_file, "%s ", argv[i]);
	}
	fprintf(sim_file, "\n\n");

	// -------------------------------
	// Print Simulation Details
	// -------------------------------
	// Simulation Mode
	#if defined(__MAGNETO)
	sprintf(sys_type, "%s", "MAG_HYDRO");
	#elif ! defined(__MAGNETO)
	sprintf(sys_type, "%s", "HYDRO");
	#else
	sprintf(sys_type, "%s", "SYS_UNKN");
	#endif
	#if defined(__INT_FAC_RK4)
	sprintf(solv_type, "%s", "INT_FAC_RK4");
	#elif defined(__RK5)
	sprintf(solv_type, "%s", "RK5");
	#elif defined(__DPRK5)
	sprintf(solv_type, "%s", "DP5");
	#else 
	sprintf(solv_type, "%s", "SOLV_UKN");
	#endif
	#if defined(PHASE_ONLY)
	sprintf(model_type, "%s", "PHAEONLY");
	#else
	sprintf(model_type, "%s", "FULL");
	#endif
	fprintf(sim_file, "Systen Type: %s\nSolver Type: %s\nModel Type: %s\n", sys_type, solv_type, model_type);

	// System Params
	fprintf(sim_file, "Viscosity: %1.10g\n", sys_vars->NU);
	fprintf(sim_file, "Re: %5.1lf\n", 1.0 / sys_vars->NU);
	fprintf(sim_file, "Magnetic Diffusivity: %1.10g\n", sys_vars->ETA);
	if (sys_vars->HYPO_MAG_DIFF_FLAG == HYPO_DIFF) {
		fprintf(sim_file, "Hypodiffusivity: YES\n");
		fprintf(sim_file, "Hypodiff Power: %1.1lf\n", HYPO_DIFF_POW);		
	}
	else {
		fprintf(sim_file, "Hypodiffusivity: NO\n");
	}
	if (sys_vars->HYPER_VISC_FLAG == HYPER_VISC){
		fprintf(sim_file, "Hyperviscosity: YES\n");
		fprintf(sim_file, "Hyperviscosity Power: %1.1lf\n", VISC_POW);	
	}
	else {
		fprintf(sim_file, "Hyperviscosity: NO\n\n");
	}
	// Equation Interaction coefficients
	fprintf(sim_file, "Velocity Interaction Coefficient: %1.3lf\n", sys_vars->EPS);
	fprintf(sim_file, "Magnetic Interaction Coefficient: %1.3lf\n\n", sys_vars->EPS_M);


	// Shell Wavenumber variables
	fprintf(sim_file, "Shell Wavenumber Prefactor: %1.3lf\n", sys_vars->k_0);
	fprintf(sim_file, "Intershell Ratio: %1.3lf\n\n", sys_vars->Lambda);

	// Spectra Slopes
	fprintf(sim_file, "Velocity Spect Slope: %1.3lf\n", sys_vars->ALPHA);
	fprintf(sim_file, "Magnetic Spect Slope: %1.3lf\n\n", sys_vars->BETA);
	
	// Spatial details
	fprintf(sim_file, "Fourier Modes: %ld\n\n", sys_vars->N);

	// Initial Conditions
	fprintf(sim_file, "Initial Conditions: %s\n\n", sys_vars->u0);

	// Forcing
	fprintf(sim_file, "Forcing Type: %s\n", sys_vars->forcing);
	fprintf(sim_file, "Forcing Wavenumber: %d\n", sys_vars->force_k);
	fprintf(sim_file, "Forcing Scale Val: %lf\n\n", sys_vars->force_scale_var);

	// Time details
	fprintf(sim_file, "Time Range: [%1.1lf - %1.1lf]\n", sys_vars->t0, sys_vars->T);
	fprintf(sim_file, "Finishing Timestep: %1.10lf\n", sys_vars->dt);
	fprintf(sim_file, "CFL No.: %1.5lf\n", sys_vars->CFL_CONST);
	if (sys_vars->ADAPT_STEP_FLAG == ADAPTIVE_STEP) {
		fprintf(sim_file, "Adaptive Stepping: YES\n");
		if (sys_vars->CFL_COND_FLAG == CFL_STEP) {
			fprintf(sim_file, "CFL Stepping Mode: YES\n");
		}
		else {
			fprintf(sim_file, "CFL Stepping Mode: NO\n");	
		}	
		fprintf(sim_file, "Total Timesteps: %ld\n", sys_vars->num_t_steps);
		fprintf(sim_file, "Total Timesteps Executed: %ld\n", sys_vars->tot_iters);
		fprintf(sim_file, "Timestep Range [min - max]: [%1.10lf - %1.10lf]\n", sys_vars->min_dt, sys_vars->max_dt);
	}
	else {
		fprintf(sim_file, "Adaptive Stepping: NO\n");
		fprintf(sim_file, "Total Timesteps: %ld\n", sys_vars->num_t_steps);
	}
	
	// Printing
	fprintf(sim_file, "Data Saved Every: %d\n", sys_vars->print_every);
	fprintf(sim_file, "Total Saving Steps: %ld\n", sys_vars->tot_save_steps);
	
	// -------------------------------
	// Print Execution Time to File
	// -------------------------------
	int hh = (int) sim_time / 3600;
	int mm = ((int )sim_time - hh * 3600) / 60;
	int ss = sim_time - hh * 3600 - mm * 60;
	fprintf(sim_file, "\n\nTotal Execution Time: %5.10lf --> %d hrs : %d mins : %d secs\n\n", sim_time, hh, mm, ss);

	// -------------------------------
	// Close File
	// -------------------------------
	fclose(sim_file);
}
/**
 * Function to carry out the signum function 
 * @param  x double Input value
 * @return   Output of performing the signum function on the input
 */
double sgn(double x) {

	if (x < 0.0) return -1.0;
	if (x > 0.0) return 1.0;
	return 0.0;
}
/**
 * Function to compute the log of base lambda 
 * @param  x double input to compute the log of
 * @return   result of computing log of the base lambda of the input
 */
double log_lambda(double x) {
	return log(x) / log(sys_vars->Lambda);
}
/**
 * Function to perform the Kronecker delta
 * @param  i Input
 * @param  j Input
 * @return   Returns the Kronecker delta of the input
 */
double my_delta(double i, double j) {

	double ans;

	if (i == j) {
		ans = 1.0;
	}
	else {
		ans = 0.0;
	}

	return ans;
}
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------