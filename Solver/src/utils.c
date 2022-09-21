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
	int dim_flag        = 0;
	int force_flag      = 0;
	int cfl_flag        = 0;
	int visc_flag       = 0;
	int drag_flag       = 0;
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
	sys_vars->dt              = 1e-4;
	sys_vars->T               = 1.0;
	sys_vars->CFL_CONST       = 0.9;
	sys_vars->ADAPT_STEP_FLAG = 0;
	sys_vars->CFL_COND_FLAG   = 1;
	// Initial conditions
	strncpy(sys_vars->u0, "N_SCALING", 64);
	sys_vars->ALPHA = 1.5;
	sys_vars->BETA  = 1.5;
	// Forcing
	strncpy(sys_vars->forcing, "NONE", 64);	
	sys_vars->force_k         = 0;
	sys_vars->force_scale_var = 1.0;
	// Transient iterations
	sys_vars->TRANS_ITERS_FLAG = 0;
	sys_vars->TRANS_ITERS_FRAC = TRANSIENT_FRAC;
	// System viscosity / hyperviscosity
	sys_vars->NU              = 0.0001;
	// sys_vars->HYPER_VISC_FLAG = 0;
	// sys_vars->HYPER_VISC_POW  = VISC_POW;
	// // System ekman dray / hypodiffusivity
	// sys_vars->EKMN_ALPHA     = 0.0;
	// sys_vars->EKMN_DRAG_FLAG = 0;
	// sys_vars->EKMN_DRAG_POW  = EKMN_POW;
	// Write to file every 
	sys_vars->SAVE_EVERY = 100;

	// -------------------------------
	// Parse CML Arguments
	// -------------------------------
	// while ((c = getopt(argc, argv, "o:h:n:d:s:e:t:v:i:c:p:f:z:T:")) != -1) {
	// 	switch(c) {
	// 		case 'o':
	// 			if (output_dir_flag == 0) {
	// 				// Read in location of output directory
	// 				strncpy(file_info->output_dir, optarg, 512);
	// 				output_dir_flag++;
	// 			}
	// 			else if (output_dir_flag == 1) {
	// 				// Output file only indicated
	// 				file_info->file_only = 1;
	// 			}
	// 			break;
	// 		case 'n':
	// 			// Read in the dimensions of the system
	// 			sys_vars->N[dim_flag] = atoi(optarg);

	// 			// Check dimensions satisfy requirements
	// 			if (sys_vars->N[dim_flag] <= 2) {
	// 				fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The minimum dimension size of [%ld] must be greater than 2\n-->> Exiting!\n\n", sys_vars->N[dim_flag]);		
	// 				exit(1);
	// 			}
	// 			else if (sys_vars->N[dim_flag] < sys_vars->rank) {
	// 				fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The system dimension of [%ld] cannot be less than the number of MPI processes of [%d]\n-->> Exiting!\n\n", sys_vars->N[dim_flag], sys_vars->rank);		
	// 				exit(1);
	// 			}
	// 			else if (sys_vars->N[dim_flag] % 2 != 0) {
	// 				fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The system dimension of [%ld] must be a multiple of 2\n-->> Exiting!\n", sys_vars->N[dim_flag]);		
	// 				exit(1);
	// 			}
	// 			dim_flag++;
	// 			break;
	// 		case 's':
	// 			// Read in intial time
	// 			sys_vars->t0 = atof(optarg);

	// 			if (sys_vars->t0 < 0) {
	// 			    fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The integration start time [%lf] must be a positive\n-->> Exiting!\n", sys_vars->t0);      
	// 			    exit(1);
	// 			}	
	// 			break;
	// 		case 'e':
	// 			// Read in final time
	// 			sys_vars->T = atof(optarg);	
	// 			if (sys_vars->T < 0) {
	// 			    fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The integration end time [%lf] must be a positive\n-->> Exiting!\n", sys_vars->T);      
	// 			    exit(1);
	// 			}
	// 			else if (sys_vars->T < sys_vars->t0) {
	// 				fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided end time: [%lf] must be greater than the initial time: [%lf]\n-->> Exiting!\n\n", sys_vars->T, sys_vars->t0);		
	// 				exit(1);
	// 			}
	// 			break;
	// 		case 'h':
	// 			if (time_step_flag == 0) {
	// 				// Read in initial timestep
	// 				sys_vars->dt = atof(optarg);
	// 				if (sys_vars->dt <= 0) {
	// 					fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided timestep: [%lf] must be strictly positive\n-->> Exiting!\n\n", sys_vars->dt);		
	// 					exit(1);
	// 				}
	// 				time_step_flag = 1;
	// 			}
	// 			else if (time_step_flag == 1) {
	// 				// Read in adaptive step indicator
	// 				sys_vars->ADAPT_STEP_FLAG = atoi(optarg);
	// 				if ((sys_vars->ADAPT_STEP_FLAG == 0) || (sys_vars->ADAPT_STEP_FLAG == 1)) {
	// 				}
	// 				else {
	// 					fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: Incorrect Adaptive stepping flag: [%d] Must be either 0 or 1 -- Set to 0 (no adaptive stepping) by default\n-->> Exiting!\n\n", sys_vars->ADAPT_STEP_FLAG);		
	// 					exit(1);
	// 				}
	// 				// time_step_flag = 2;
	// 				break;
	// 			}
	// 			break;			
	// 		case 'c':
	// 			if (cfl_flag == 0) {
	// 				// Read in value of the CFL -> this can be used to control the timestep
	// 				sys_vars->CFL_COND_FLAG = atoi(optarg);
	// 				if ((sys_vars->CFL_COND_FLAG == 0) || (sys_vars->CFL_COND_FLAG == 1)) {
	// 				}
	// 				else {
	// 					fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: Incorrect CFL condition flag: [%d] Must be either 0 or 1 -- Set to 1 (use CFL condition) by default\n-->> Exiting!\n\n", sys_vars->CFL_COND_FLAG);		
	// 					exit(1);
	// 				}
	// 				cfl_flag = 1;
	// 				break;	
	// 			}
	// 			else if (cfl_flag == 1) {
	// 				sys_vars->CFL_CONST = atof(optarg);
	// 				if (sys_vars->CFL_CONST <= 0) {
	// 					fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided CFL Constant: [%lf] must be strictly positive\n-->> Exiting!\n\n", sys_vars->CFL_CONST);		
	// 					exit(1);
	// 				}
	// 				// cfl_flag = 2;
	// 				break;
	// 			}
	// 			break;				
	// 		case 'v':
	// 			if (visc_flag == 0) {
	// 				// Read in the viscosity
	// 				sys_vars->NU = atof(optarg);
	// 				if (sys_vars->NU <= 0) {
	// 					fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided viscosity: [%lf] must be strictly positive\n-->> Exiting!\n\n", sys_vars->NU);		
	// 					exit(1);
	// 				}
	// 				visc_flag = 1;
	// 				break;
	// 			}
	// 			else if (visc_flag == 1) {
	// 				// Read in hyperviscosity flag
	// 				sys_vars->HYPER_VISC_FLAG = atoi(optarg);
	// 				if ((sys_vars->HYPER_VISC_FLAG == 0) || (sys_vars->HYPER_VISC_FLAG == 1)) {
	// 				}
	// 				else {
	// 					fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: Incorrect Hyperviscosity flag: [%d] Must be either 0 or 1 -- Set to 0 (no Hyperviscosity) by default\n-->> Exiting!\n\n", sys_vars->HYPER_VISC_FLAG);		
	// 					exit(1);
	// 				}
	// 				visc_flag = 2;
	// 				break;	
	// 			}
	// 			else if (visc_flag == 2) {
	// 				// Read in the hyperviscosity power
	// 				sys_vars->HYPER_VISC_POW = atof(optarg);
	// 				if (sys_vars->HYPER_VISC_POW <= 0.0) {
	// 					fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided hyperviscosity power: [%lf] must be strictly positive\n-->> Exiting!\n\n", sys_vars->HYPER_VISC_POW);		
	// 					exit(1);
	// 				}
	// 				// visc_flag = 3;
	// 				break;
	// 			}
	// 			break;
	// 		case 'd':
	// 			if (drag_flag == 0) {
	// 				// Read in the drag
	// 				sys_vars->EKMN_ALPHA = atof(optarg);
	// 				if (sys_vars->EKMN_ALPHA < 0) {
	// 					fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided Ekman Drag: [%lf] must be positive\n-->> Exiting!\n\n", sys_vars->EKMN_ALPHA);		
	// 					exit(1);
	// 				}
	// 				drag_flag = 1;
	// 				break;
	// 			}
	// 			else if (drag_flag == 1) {
	// 				// Read in Ekman drag flag
	// 				sys_vars->EKMN_DRAG_FLAG = atoi(optarg);
	// 				if ((sys_vars->EKMN_DRAG_FLAG == 0) || (sys_vars->EKMN_DRAG_FLAG == 1)) {
	// 				}
	// 				else {
	// 					fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: Incorrect Ekman Drag flag: [%d] Must be either 0 or 1 -- Set to 0 (no drag) by default\n-->> Exiting!\n\n", sys_vars->EKMN_DRAG_FLAG);		
	// 					exit(1);
	// 				}
	// 				drag_flag = 2;
	// 				break;	
	// 			}
	// 			else if (drag_flag == 2) {
	// 				// Read in the hypodiffusivity power
	// 				sys_vars->EKMN_DRAG_POW = atof(optarg);
	// 				if (sys_vars->EKMN_DRAG_POW >= 0.0) {
	// 					fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The provided Hypodiffusibity power: [%lf] must be strictly negative\n-->> Exiting!\n\n", sys_vars->EKMN_DRAG_POW);		
	// 					exit(1);
	// 				}
	// 				// drag_flag = 3;
	// 				break;
	// 			}
	// 			break;
	// 		case 'i':
	// 			// Read in the initial conditions
	// 			if (!(strcmp(optarg,"TG_VEL"))) {
	// 				// The Taylor Green vortex - starting with the velocity
	// 				strncpy(sys_vars->u0, "TG_VEL", 64);
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"TG_VORT"))) {
	// 				// The Taylor Green vortex - starting with the vorticity
	// 				strncpy(sys_vars->u0, "TG_VORT", 64);
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"DOUBLE_SHEAR_LAYER"))) {
	// 				// Double Shear Layer - for testing the Euler solver
	// 				strncpy(sys_vars->u0, "DOUBLE_SHEAR_LAYER", 64);
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"TESTING"))) {
	// 				// Specific ICs for testing - powerlaw amps and phases = pi/4
	// 				strncpy(sys_vars->u0, "TESTING", 64);
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"DECAY_TURB_BB"))) {
	// 				// Decay Turbulence -> Broadband
	// 				strncpy(sys_vars->u0, "DECAY_TURB_BB", 64);
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"DECAY_TURB_NB"))) {
	// 				// Decay Turbulence -> Narrow Band
	// 				strncpy(sys_vars->u0, "DECAY_TURB_NB", 64);
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"DECAY_TURB_ALT"))) {
	// 				// Decay Turbulence -> McWilliams 1984
	// 				strncpy(sys_vars->u0, "DECAY_TURB_ALT", 64);
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"DECAY_TURB_EXP"))) {
	// 				// Decay Turbulence -> Exponential Spectrum
	// 				strncpy(sys_vars->u0, "DECAY_TURB_EXP", 64);
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"GAUSS_DECAY_TURB"))) {
	// 				// Decay Turbulence -> M. Brachet et al. 1988 
	// 				strncpy(sys_vars->u0, "GAUSS_DECAY_TURB", 64);
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"EXTRM_ENS"))) {
	// 				// Extreme Enstrophy Spectrum Initial condition
	// 				strncpy(sys_vars->u0, "EXTRM_ENS", 64);
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"MAX_PALIN"))) {
	// 				// Maximal Palinstorphy growth initial condition -> Ayala & Protas 2014
	// 				strncpy(sys_vars->u0, "MAX_PALIN", 64);
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"RING"))) {
	// 				// Ring initial condition
	// 				strncpy(sys_vars->u0, "RING", 64);
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"RANDOM"))) {
	// 				// Random initial conditions
	// 				strncpy(sys_vars->u0, "RANDOM", 64);
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"GAUSS_BLOB"))) {
	// 				// Gaussian Blob conditions
	// 				strncpy(sys_vars->u0, "GAUSS_BLOB", 64);
	// 				break;
	// 			}
	// 			else {
	// 				// No initial conditions specified -> this will default to random initial conditions
	// 				strncpy(sys_vars->u0, "NONE", 64);
	// 				break;
	// 			}
	// 			break;	
	// 		case 't':
	// 			// Read in output directory tag
	// 			strncpy(file_info->output_tag, optarg, 64);	
	// 			break;
	// 		case 'T':
	// 			if (trans_iter_flag == 0) {
	// 				// Read in transient iterations flag
	// 				sys_vars->TRANS_ITERS_FLAG = atoi(optarg);
	// 				if ((sys_vars->TRANS_ITERS_FLAG == 0) || (sys_vars->TRANS_ITERS_FLAG == 1)) {
	// 				}
	// 				else {
	// 					fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: Incorrect transient iterations flag: [%d] Must be either 0 or 1 -- Set to 0 (no transient iterations) by default\n-->> Exiting!\n\n", sys_vars->TRANS_ITERS_FLAG);		
	// 					exit(1);
	// 				}
	// 				trans_iter_flag = 1;
	// 				break;	
	// 			}
	// 			else if (trans_iter_flag == 1) {
	// 				// Read in transient iterations fraction of total iterations
	// 				sys_vars->TRANS_ITERS_FRAC = atof(optarg);
	// 				if ((sys_vars->TRANS_ITERS_FRAC >= 1.0) || (sys_vars->TRANS_ITERS_FRAC <= 0.0)) {
	// 					fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: Incorrect transient iterations fraction: [%lf] Must be between 0 and 1 -- Set to 0.2 by default\n-->> Exiting!\n\n", sys_vars->TRANS_ITERS_FRAC);		
	// 					exit(1);
	// 				}
	// 				// trans_iter_flag = 2;
	// 				break;
	// 			}
	// 			break;				
	// 		case 'p':
	// 			// Read in how often to print to file
	// 			sys_vars->SAVE_EVERY = atoi(optarg);
	// 			break;
	// 		case 'f':
	// 			// Read in the forcing type
	// 			if (!(strcmp(optarg,"ZERO")) && (force_flag == 0)) {
	// 				// Killing certain modes
	// 				strncpy(sys_vars->forcing, "ZERO", 64);
	// 				force_flag = 1;
	// 				break;
	// 			}
	// 			if (!(strcmp(optarg,"CONST_GAUSS")) && (force_flag == 0)) {
	// 				// Deterministic, constant in time, forcing of low wavenumbers
	// 				strncpy(sys_vars->forcing, "CONST_GAUSS", 64);
	// 				force_flag = 1;
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"KOLM"))  && (force_flag == 0)) {
	// 				// Kolmogorov forcing
	// 				strncpy(sys_vars->forcing, "KOLM", 64);
	// 				force_flag = 1;
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"NONE"))  && (force_flag == 0)) {
	// 				// No forcing
	// 				strncpy(sys_vars->forcing, "NONE", 64);
	// 				force_flag = 1;
	// 				break;
	// 			}
	// 			else if (!(strcmp(optarg,"STOC"))  && (force_flag == 0)) {
	// 				// Stochastic forcing
	// 				strncpy(sys_vars->forcing, "STOC", 64);
	// 				force_flag = 1;
	// 				break;
	// 			}
	// 			else if ((force_flag == 1)) {
	// 				// Get the forcing wavenumber
	// 				sys_vars->force_k = atoi(optarg);
	// 				force_flag = 2;
	// 				break;
	// 			}
	// 			else if ((force_flag == 2)) {
	// 				// Get the force scaling variable
	// 				sys_vars->force_scale_var = atof(optarg);
	// 				break;
	// 			}
	// 			break;
	// 		case 'z':
	// 			strncpy((file_info->input_file_name), optarg, 512);	// copy the input file name given as a command line argument
	// 			if ( access((file_info->input_file_name), F_OK) != 0) {
	// 				fprintf(stderr, "\n["RED"ERROR"RESET"] Parsing of Command Line Arguements Failed: The input file [%s] cannot be found, please ensure correct path to file is specified.\n", (file_info->input_file_name));		
	// 				exit(1);					
	// 			}
	// 			break;
	// 		default:
	// 			fprintf(stderr, "\n["RED"ERROR"RESET"] Incorrect command line flag encountered\n");		
	// 			fprintf(stderr, "Use"YELLOW" -o"RESET" to specify the output directory\n");
	// 			fprintf(stderr, "Use"YELLOW" -n"RESET" to specify the size of each dimension in the system\n");
	// 			fprintf(stderr, "Use"YELLOW" -s"RESET" to specify the start time of the simulation\n");
	// 			fprintf(stderr, "Use"YELLOW" -e"RESET" to specify the end time of the simulation\n");
	// 			fprintf(stderr, "Use"YELLOW" -h"RESET" to specify the timestep\n");
	// 			fprintf(stderr, "Use"YELLOW" -c"RESET" to specify the CFL constant for the adaptive stepping\n");
	// 			fprintf(stderr, "Use"YELLOW" -v"RESET" to specify the system viscosity\n");
	// 			fprintf(stderr, "Use"YELLOW" -i"RESET" to specify the initial condition\n");
	// 			fprintf(stderr, "Use"YELLOW" -t"RESET" to specify the tag name to be used in the output file directory\n");
	// 			fprintf(stderr, "Use"YELLOW" -f"RESET" to specify the forcing type\n");
	// 			fprintf(stderr, "Use"YELLOW" -T"RESET" to specify if transient iterations are to be performed\n");
	// 			fprintf(stderr, "Use"YELLOW" -p"RESET" to specify how often to print to file\n");
	// 			fprintf(stderr, "Use"YELLOW" -z"RESET" to specify an input file to read parameters from\n");
	// 			fprintf(stderr, "\nExample usage:\n"CYAN"\tmpirun -n 4 ./bin/main -o \"../Data/Tmp\" -n 64 -n 64 -s 0.0 -e 1.0 -h 0.0001 -v 1.0 -i \"TG_VORT\" -t \"TEMP_RUN\" \n"RESET);
	// 			fprintf(stderr, "-->> Now Exiting!\n\n");
	// 			exit(1);
	// 	}
	// }

	return 0;
}
// /**
//  * Function that prints the summary details of the simulation to a .txt
//  * @param sim_time The execution time of the simulation 
//  * @param argc     The number of command line arguments
//  * @param argv     Array conatining the command line arguments
//  */
// void PrintSimulationDetails(int argc, char** argv, double sim_time) {

// 	// Initialize variables
// 	FILE *sim_file;
// 	char sys_type[64];
// 	char solv_type[64];
// 	char model_type[64];
// 	char dealias_type[64];
// 	char file_path[512];

// 	// -------------------------------
// 	// Open File
// 	// -------------------------------
// 	strcpy(file_path, file_info->output_dir);
// 	strcat(file_path, "SimulationDetails.txt");
// 	sim_file = fopen(file_path, "w");

// 	// -------------------------------
// 	// Print Executing Command
// 	// -------------------------------
// 	fprintf(sim_file, "Executing Command:"); 
// 	fprintf(sim_file, "\n\n\tmpirun -n %d ", sys_vars->num_procs);
// 	for (int i = 0; i < argc; ++i) {
// 		fprintf(sim_file, "%s ", argv[i]);
// 	}
// 	fprintf(sim_file, "\n\n");

// 	// -------------------------------
// 	// Print Simulation Details
// 	// -------------------------------
// 	// Simulation Mode
// 	#if defined(__NAVIER)
// 	sprintf(sys_type, "%s", "NAVIER");
// 	#elif defined(__EULER)
// 	sprintf(sys_type, "%s", "EULER");
// 	#else
// 	sprintf(sys_type, "%s", "SYS_UNKN");
// 	#endif
// 	#if defined(__RK4)
// 	sprintf(solv_type, "%s", "RK4");
// 	#elif defined(__RK5)
// 	sprintf(solv_type, "%s", "RK5");
// 	#elif defined(__DPRK5)
// 	sprintf(solv_type, "%s", "DP5");
// 	#else 
// 	sprintf(solv_type, "%s", "SOLV_UKN");
// 	#endif
// 	#if defined(__PHASE_ONLY)
// 	sprintf(model_type, "%s", "PHAEONLY");
// 	#else
// 	sprintf(model_type, "%s", "FULL");
// 	#endif
// 	fprintf(sim_file, "Systen Type: %s\nSolver Type: %s\nModel Type: %s\n", sys_type, solv_type, model_type);

// 	// System Params
// 	fprintf(sim_file, "Viscosity: %1.6lf\n", sys_vars->NU);
// 	fprintf(sim_file, "Re: %5.1lf\n", 1.0 / sys_vars->NU);
// 	#if defined(__EKMN_DRAG)
// 	fprintf(sim_file, "Ekman Drag: YES\n");
// 	fprintf(sim_file, "Ekman Alpha: %1.4lf\n", sys_vars->EKMN_ALPHA);
// 	fprintf(sim_file, "Ekman Power: %d\n", EKMN_POW);
// 	#else
// 	fprintf(sim_file, "Ekman Drag: NO\n");
// 	#endif
// 	#if defined(__HYPER)
// 	fprintf(sim_file, "Hyperviscosity: YES\n");
// 	fprintf(sim_file, "Hyperviscosity Power: %1.1lf\n", VIS_POW);	
// 	#else
// 	fprintf(sim_file, "Hyperviscosity: NO\n");
// 	#endif

// 	// Spatial details
// 	fprintf(sim_file, "\nCollocation Points: [%ld, %ld]\n", sys_vars->N[0], sys_vars->N[1]);
// 	fprintf(sim_file, "Spatial Increment: [%1.4lf, %1.4lf]\n", sys_vars->dx, sys_vars->dy);
// 	fprintf(sim_file, "Fourier Modes: [%ld, %ld]\n\n", sys_vars->N[0], sys_vars->N[1]/2 + 1);

// 	// Initial Conditions
// 	fprintf(sim_file, "Initial Conditions: %s\n", sys_vars->u0);
	
// 	// Dealising details
// 	#if defined(__DEALIAS_HOU_LI)
// 	sprintf(dealias_type, "%s", "HOU-LI");
// 	#elif defined(__DEALIAS_23)
// 	sprintf(dealias_type, "%s", "2/3 RDS");
// 	#else
// 	sprintf(dealias_type, "%s", "UNKNOWN");
// 	#endif
// 	fprintf(sim_file, "Dealiasing: %s\n", dealias_type);
	
// 	// Forcing
// 	fprintf(sim_file, "Forcing Type: %s\n\n", sys_vars->forcing);

// 	// Time details
// 	fprintf(sim_file, "Time Range: [%1.1lf - %1.1lf]\n", sys_vars->t0, sys_vars->T);
// 	fprintf(sim_file, "Finishing Timestep: %1.10lf\n", sys_vars->dt);
// 	fprintf(sim_file, "CFL No.: %1.5lf\n", sys_vars->CFL_CONST);
// 	if (sys_vars->ADAPT_STEP_FLAG == ADAPTIVE_STEP) {
// 		fprintf(sim_file, "Adaptive Stepping: YES\n");
// 			if (sys_vars->CFL_COND_FLAG == CFL_STEP) {
// 				fprintf(sim_file, "CFL Stepping Mode: YES\n");
// 			}
// 			else {
// 				fprintf(sim_file, "CFL Stepping Mode: NO\n");	
// 			}	
// 		fprintf(sim_file, "Total Timesteps: %ld\n", sys_vars->num_t_steps);
// 		fprintf(sim_file, "Total Timesteps Executed: %ld\n", sys_vars->tot_iters);
// 		fprintf(sim_file, "Timestep Range [min - max]: [%1.10lf - %1.10lf]\n", sys_vars->min_dt, sys_vars->max_dt);
// 	}
// 	else {
// 		fprintf(sim_file, "Adaptive Stepping: NO\n");
// 		fprintf(sim_file, "Total Timesteps: %ld\n", sys_vars->num_t_steps);
// 	}
	
// 	// Printing
// 	fprintf(sim_file, "Data Saved Every: %d\n", sys_vars->print_every);
// 	fprintf(sim_file, "Total Saving Steps: %ld\n", sys_vars->tot_save_steps);

// 	// Flux subset details
// 	#if defined(__SPECT)
// 	fprintf(sim_file, "\nFlux Subset Details: \n\tLower Limit: %d\n\tUpper Limit: %d\n", LWR_SBST_LIM, UPR_SBST_LIM);
// 	#endif

	
// 	// -------------------------------
// 	// Print Execution Time to File
// 	// -------------------------------
// 	int hh = (int) sim_time / 3600;
// 	int mm = ((int )sim_time - hh * 3600) / 60;
// 	int ss = sim_time - hh * 3600 - mm * 60;
// 	fprintf(sim_file, "\n\nTotal Execution Time: %5.10lf --> %d hrs : %d mins : %d secs\n\n", sim_time, hh, mm, ss);

// 	// -------------------------------
// 	// Close File
// 	// -------------------------------
// 	fclose(sim_file);
// }
// /**
//  * Utility function that gathers all the local real space vorticity arrays on the master process and prints to screen 
//  * @param N Array containing the dimensions of the system
//  */
// void PrintVorticityReal(const long int* N) {

// 	// Initialize variables
// 	const long int Nx = N[0];
// 	const long int Ny = N[1];

// 	// Allocate memory to recieve lcoal vorticity arrays
// 	double* w0 = (double* )fftw_malloc(sizeof(double) * Nx * (Ny + 2));
// 	for (int i = 0; i < Nx; ++i) {
// 		for (int j = 0; j < Ny + 2; ++j) {
// 			w0[i * (Ny + 2) + j] = 0.0;
// 		}
// 	}

// 	// Gather all the local arrays into w0
// 	MPI_Gather(run_data->w, sys_vars->local_Nx * (Ny + 2), MPI_DOUBLE, w0, sys_vars->local_Nx * (Ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);

// 	// On the master rank print the result
// 	if ( !(sys_vars->rank) ) {
// 		for (int i = 0; i < Nx; ++i) {
// 			for (int j = 0; j < Ny; ++j) {
// 				printf("w[%ld]: %+5.16lf\t", i * (Ny) + j, w0[i * (Ny + 2) + j]);
// 			}
// 			printf("\n");
// 		}
// 		printf("\n\n");
// 	}

// 	// Free temporary memory
// 	fftw_free(w0);
// }
// /**
//  * Utility function that gathers all the local Fourier space vorticity arrays on the master process and prints to screen 
//  * @param N Array containing the dimensions of the system
//  */
// void PrintVorticityFourier(const long int* N) {

// 	// Initialize variables
// 	const long int Nx 		  = N[0];
// 	const long int Ny_Fourier = N[1] / 2 + 1;

// 	// Allocate memory to recieve lcoal vorticity arrays
// 	fftw_complex* w_hat0 = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * Nx * (Ny_Fourier));
// 	for (int i = 0; i < Nx; ++i) {
// 		for (int j = 0; j < Ny_Fourier; ++j) {
// 			w_hat0[i * (Ny_Fourier) + j] = 0.0 + 0.0 * I;
// 		}
// 	}
// 	// Gather all the local arrays into w0
// 	MPI_Gather(run_data->w_hat, sys_vars->local_Nx * Ny_Fourier, MPI_C_DOUBLE_COMPLEX, w_hat0, sys_vars->local_Nx * Ny_Fourier, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

// 	// On the master rank print the result
// 	if ( !(sys_vars->rank) ) {
// 		for (int i = 0; i < Nx; ++i) {
// 			for (int j = 0; j < Ny_Fourier; ++j) {
// 				printf("wh[%ld]: %+5.16lf %+5.16lfI\t", i * (Ny_Fourier) + j, creal(w_hat0[i * Ny_Fourier + j]), cimag(w_hat0[i * Ny_Fourier + j]));
// 			}
// 			printf("\n");
// 		}
// 		printf("\n\n");
// 	}

// 	// Free temporary memory
// 	fftw_free(w_hat0);
// }
// /**
//  * Utility function that gathers all the local real space velocity arrays on the master process and prints to screen 
//  * @param N Array containing the dimensions of the system
//  */
// void PrintVelocityReal(const long int* N) {

// 	// Initialize variables
// 	const long int Nx = N[0];
// 	const long int Ny = N[1];

// 	// Allocate memory to recieve lcoal vorticity arrays
// 	double* u0 = (double* )fftw_malloc(sizeof(double) * Nx * (Ny + 2) * SYS_DIM);
// 	for (int i = 0; i < Nx; ++i) {
// 		for (int j = 0; j < Ny + 2; ++j) {
// 			u0[SYS_DIM * (i * (Ny + 2) + j) + 0] = 0.0;
// 			u0[SYS_DIM * (i * (Ny + 2) + j) + 1] = 0.0;
// 		}
// 	}

// 	// Gather all the local arrays into u0 and v0
// 	MPI_Gather(run_data->u, sys_vars->local_Nx * (Ny + 2) * SYS_DIM, MPI_DOUBLE, u0, sys_vars->local_Nx * (Ny + 2)* SYS_DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
// 	// On the master rank print the result
// 	if ( !(sys_vars->rank) ) {
// 		for (int i = 0; i < Nx; ++i) {
// 			for (int j = 0; j < Ny; ++j) {
// 				printf("u[%ld]: %+5.16lf\t", i * (Ny) + j, u0[SYS_DIM * (i * (Ny + 2) + j) + 0]);
// 			}
// 			printf("\n");
// 		}
// 		printf("\n\n");

// 		for (int i = 0; i < Nx; ++i) {
// 			for (int j = 0; j < Ny; ++j) {
// 				printf("v[%ld]: %+5.16lf\t", i * (Ny) + j, u0[SYS_DIM * (i * (Ny + 2) + j) + 1]);
// 			}
// 			printf("\n");
// 		}
// 		printf("\n\n");
// 	}

// 	// Free temporary memory
// 	fftw_free(u0);
// }
// /**
//  * Utility function that gathers all the local Fourier space velocity arrays on the master process and prints to screen 
//  * @param N Array containing the dimensions of the system
//  */
// void PrintVelocityFourier(const long int* N) {

// 	// Initialize variables
// 	const long int Nx 		  = N[0];
// 	const long int Ny_Fourier = N[1] / 2 + 1;

// 	// Allocate memory to recieve lcoal vorticity arrays
// 	fftw_complex* u_hat0 = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * Nx * Ny_Fourier * SYS_DIM);
// 	for (int i = 0; i < Nx; ++i) {
// 		for (int j = 0; j < Ny_Fourier; ++j) {
// 			u_hat0[SYS_DIM * (i * (Ny_Fourier) + j) + 0] = 0.0 + 0.0 * I;
// 			u_hat0[SYS_DIM * (i * (Ny_Fourier) + j) + 1] = 0.0 + 0.0 * I;
// 		}
// 	}

// 	// Gather all the local arrays into u_hat0 and v_hat0
// 	MPI_Gather(run_data->u_hat, sys_vars->local_Nx * Ny_Fourier * SYS_DIM, MPI_C_DOUBLE_COMPLEX, u_hat0, sys_vars->local_Nx * Ny_Fourier* SYS_DIM, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	
// 	// On the master rank print the result
// 	if ( !(sys_vars->rank) ) {
// 		for (int i = 0; i < Nx; ++i) {
// 			for (int j = 0; j < Ny_Fourier; ++j) {
// 				printf("uh[%ld]: %+5.10lf %+5.10lfI\t", i * (Ny_Fourier) + j, creal(u_hat0[SYS_DIM * (i * Ny_Fourier + j) + 0]), cimag(u_hat0[SYS_DIM * (i * Ny_Fourier + j) + 0]));
// 			}
// 			printf("\n");
// 		}
// 		printf("\n\n");

// 		for (int i = 0; i < Nx; ++i) {
// 			for (int j = 0; j < Ny_Fourier; ++j) {
// 				printf("vh[%ld]: %+5.10lf %+5.10lf I\t", i * (Ny_Fourier) + j, creal(u_hat0[SYS_DIM * (i * Ny_Fourier + j) + 1]), cimag(u_hat0[SYS_DIM * (i * Ny_Fourier + j) + 1]));
// 			}
// 			printf("\n");
// 		}
// 		printf("\n\n");
// 	}

// 	// Free temporary memory
// 	fftw_free(u_hat0);
// }
// /**
//  * Function to print scalar array in Fourier space
//  * @param data The data to be printed
//  * @param N    The size of the dimensions
//  */
// void PrintScalarFourier(fftw_complex* data, const long int* N, char* arr_name) {

// 	// Initialize variables
// 	const long int Nx 		  = N[0];
// 	const long int Ny_Fourier = N[1] / 2 + 1;

// 	// Allocate memory to recieve lcoal vorticity arrays
// 	fftw_complex* w_hat0 = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * Nx * (Ny_Fourier));
// 	for (int i = 0; i < Nx; ++i) {
// 		for (int j = 0; j < Ny_Fourier; ++j) {
// 			w_hat0[(i * (Ny_Fourier) + j)] = 0.0 + 0.0 * I;
// 		}
// 	}

// 	// Gather all the local arrays into w0
// 	MPI_Gather(data, sys_vars->local_Nx * Ny_Fourier, MPI_C_DOUBLE_COMPLEX, w_hat0, sys_vars->local_Nx * Ny_Fourier, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

// 	// On the master rank print the result
// 	if ( !(sys_vars->rank) ) {
// 		for (int i = 0; i < Nx; ++i) {
// 			for (int j = 0; j < Ny_Fourier; ++j) {
// 				printf("%s[%ld]: %1.16lf %1.16lf I\n", arr_name, i * (Ny_Fourier) + j, creal(w_hat0[i * Ny_Fourier + j]), cimag(w_hat0[i * Ny_Fourier + j]));
// 			}
// 			printf("\n");
// 		}
// 		printf("\n\n");
// 	}

// 	// Free temporary memory
// 	fftw_free(w_hat0);
// }
// /**
//  * Function to print scalar array in Real space
//  * @param data The data to be printed
//  * @param N    The size of the dimensions
//  */
// void PrintScalarReal(double* data, const long int* N, char* arr_name) {

// 	// Initialize variables
// 	const long int Nx = N[0];
// 	const long int Ny = N[1];

// 	// Allocate memory to recieve lcoal vorticity arrays
// 	double* w_hat0 = (double* )fftw_malloc(sizeof(double) * Nx * (Ny + 2));
// 	for (int i = 0; i < Nx; ++i) {
// 		for (int j = 0; j < Ny + 2; ++j) {
// 			w_hat0[(i * (Ny + 2) + j)] = 0.0;
// 		}
// 	}

// 	// Gather all the local arrays into w0
// 	MPI_Gather(data, sys_vars->local_Nx * (Ny + 2), MPI_DOUBLE, w_hat0, sys_vars->local_Nx * (Ny + 2), MPI_DOUBLE, 0, MPI_COMM_WORLD);

// 	// On the master rank print the result
// 	if ( !(sys_vars->rank) ) {
// 		for (int i = 0; i < Nx; ++i) {
// 			for (int j = 0; j < Ny; ++j) {
// 				printf("%s[%ld]: %g ", arr_name, i * (Ny) + j, w_hat0[i * (Ny + 2) + j]);
// 			}
// 			printf("\n");
// 		}
// 		printf("\n\n");
// 	}

// 	// Free temporary memory
// 	fftw_free(w_hat0);
// }
// /**
//  * Function for printing a vector array in Real space
//  * @param data data to be printed
//  * @param N    The size of the dimensions
//  */
// void PrintVectorReal(double* data, const long int* N, char* arr_name1, char* arr_name2) {

// 	// Initialize variables
// 	const long int Nx = N[0];
// 	const long int Ny = N[1];

// 	// Allocate memory to recieve lcoal vorticity arrays
// 	double* u0 = (double* )fftw_malloc(sizeof(double) * Nx * (Ny + 2) * SYS_DIM);
// 	for (int i = 0; i < Nx; ++i) {
// 		for (int j = 0; j < Ny + 2; ++j) {
// 			u0[SYS_DIM * (i * (Ny + 2) + j) + 0] = 0.0;
// 			u0[SYS_DIM * (i * (Ny + 2) + j) + 1] = 0.0;
// 		}
// 	}

// 	// Gather all the local arrays into u0 and v0
// 	MPI_Gather(data, sys_vars->local_Nx * (Ny + 2) * SYS_DIM, MPI_DOUBLE, u0, sys_vars->local_Nx * (Ny + 2)* SYS_DIM, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
// 	// On the master rank print the result
// 	if ( !(sys_vars->rank) ) {
// 		for (int i = 0; i < Nx; ++i) {
// 			for (int j = 0; j < Ny; ++j) {
// 				printf("%s[%ld]: %+5.16lf", arr_name1, i * (Ny) + j, u0[SYS_DIM * (i * (Ny + 2) + j) + 0]);
// 			}
// 			printf("\n");
// 		}
// 		printf("\n\n");

// 		for (int i = 0; i < Nx; ++i) {
// 			for (int j = 0; j < Ny; ++j) {
// 				printf("%s[%ld]: %+5.16lf", arr_name2, i * (Ny) + j, u0[SYS_DIM * (i * (Ny + 2) + j) + 1]);
// 			}
// 			printf("\n");
// 		}
// 		printf("\n\n");
// 	}

// 	// Free temporary memory
// 	fftw_free(u0);
// }
// void PrintVectorFourier(fftw_complex* data, const long int* N, char* arr_name1, char* arr_name2) {

// 	// Initialize variables
// 	const long int Nx = N[0];
// 	const long int Ny_Fourier = N[1] / 2 + 1;

// 	// Allocate memory to recieve lcoal vorticity arrays
// 	fftw_complex* u_hat0 = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * Nx * (Ny_Fourier) * SYS_DIM);
// 	for (int i = 0; i < Nx; ++i) {
// 		for (int j = 0; j < Ny_Fourier; ++j) {
// 			u_hat0[SYS_DIM * (i * (Ny_Fourier) + j) + 0] = 0.0 + 0.0 * I;
// 			u_hat0[SYS_DIM * (i * (Ny_Fourier) + j) + 1] = 0.0 + 0.0 * I;
// 		}
// 	}

// 	// Gather all the local arrays into u_hat0 and v_hat0
// 	MPI_Gather(data, sys_vars->local_Nx * (Ny_Fourier) * SYS_DIM, MPI_C_DOUBLE_COMPLEX, u_hat0, sys_vars->local_Nx * (Ny_Fourier)* SYS_DIM, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

// 	// On the master rank print the result
// 	if ( !(sys_vars->rank) ) {
// 		for (int i = 0; i < Nx; ++i) {
// 			for (int j = 0; j < Ny_Fourier; ++j) {
// 				printf("%s[%ld]: %5.16lf %5.16lf I", arr_name1, i * (Ny_Fourier) + j, creal(u_hat0[SYS_DIM * (i * (Ny_Fourier) + j) + 0]), cimag(u_hat0[SYS_DIM * (i * (Ny_Fourier) + j) + 0]));
// 			}
// 			printf("\n");
// 		}
// 		printf("\n\n");

// 		for (int i = 0; i < Nx; ++i) {
// 			for (int j = 0; j < Ny_Fourier; ++j) {
// 				printf("%s[%ld]: %5.16lf %5.16lf I", arr_name2, i * (Ny_Fourier) + j, creal(u_hat0[SYS_DIM * (i * (Ny_Fourier) + j) + 1]), cimag(u_hat0[SYS_DIM * (i * (Ny_Fourier) + j) + 0]));
// 			}
// 			printf("\n");
// 		}
// 		printf("\n\n");
// 	}

// 	// Free temporary memory
// 	fftw_free(u_hat0);
// }
// /**
//  * Function to print the Real and Fourier space variables
//  */
// void PrintSpaceVariables(const long int* N) {

// 	// Initialize variables
// 	const long int Nx = N[0];
// 	const long int Ny_Fourier = N[1] / 2 + 1; 

// 	// Allocate global array memory
// 	double* x0 = (double* )fftw_malloc(sizeof(double) * Nx);
// 	int* k0    = (int* )fftw_malloc(sizeof(int) * Nx);

// 	// Gather the data from each process onto master rank for printing
// 	MPI_Gather(run_data->x[0], sys_vars->local_Nx, MPI_DOUBLE, x0, sys_vars->local_Nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);
// 	MPI_Gather(run_data->k[0], sys_vars->local_Nx, MPI_INT, k0, sys_vars->local_Nx, MPI_INT, 0, MPI_COMM_WORLD);

// 	// Print results on the master process
// 	if ( !(sys_vars->rank) ) {
// 		for (int i = 0; i < Nx; ++i) {
// 			if (i < Ny_Fourier) {
// 				printf("x[%d]: %5.16lf\ty[%d]: %5.16lf\tkx[%d]: %d\tky[%d]: %d\n", i, x0[i], i, run_data->x[1][i], i, k0[i], i, run_data->k[1][i]);
// 			}
// 			else {
// 				printf("x[%d]: %5.16lf\ty[%d]: %5.16lf\tkx[%d]: %d\n", i, x0[i], i, run_data->x[1][i], i, k0[i]);
// 			}
// 		}
// 		printf("\n\n");
// 	}

// 	// Free up memory
// 	fftw_free(x0);
// 	fftw_free(k0);
// }
// /**
//  * Function used in testing /debugging to write Real datasets to file
//  * @param data        The data to be written to file
//  * @param dset_name   The name of the dataset
//  * @param dset_rank   The rank of the dataset
//  * @param dset_dims   Array containing the dimensions of the dataset to be written
//  * @param local_dim_x The size of the first dimension of the lcoal dataset
//  */
// void WriteTestDataReal(double* data, char* dset_name, int dset_rank, int* dset_dims, int local_dim_x) {

// 	// Initialize variables
// 	hsize_t rank = dset_rank;
//     hsize_t Dims[rank];
//     herr_t status;
    
//     // Allocate array for gathering full dataset on root process
//     int mem_size = 1;
//     for (int i = 0; i < (int)rank; ++i) {
//     	mem_size *= dset_dims[i];
//     }
// 	double* full_data = (double* )fftw_malloc(sizeof(double) * mem_size);
// 	if (full_data == NULL) {
// 		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for writing ["CYAN"%s"RESET"] to test data file\n-->> Exiting!!!\n", dset_name);
// 		exit(1);
// 	}	

// 	// Gather full dataset on root process
// 	int local_size = local_dim_x;
// 	for (int i = 1; i < (int)rank; ++i) {
// 		local_size *= dset_dims[i];
// 	}
// 	MPI_Gather(data, local_size, MPI_DOUBLE, full_data, local_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
// 	// Write data to file on root process
// 	if (!sys_vars->rank) {   
// 	    for (int i = 0; i < (int)rank; ++i) {
// 	    	Dims[i] = dset_dims[i];
// 	    }
// 	    status = H5LTmake_dataset(file_info->test_file_handle, dset_name, rank, Dims, H5T_NATIVE_DOUBLE, full_data);
// 	    if (status < 0) {
// 	    	fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to write ["CYAN"%s"RESET"] to test data file\n-->> Exiting!!!\n", dset_name);
// 	    	exit(1);
// 	    }
// 	}

// 	// Free memory
// 	fftw_free(full_data);
// }
// /**
//  * Function used in testing / debugging to write a Fourier dataset to file
//  * @param data        The data to be written to file
//  * @param dset_name   The name of the dataset to be written	
//  * @param dim_x       The size of the first dimension of the global dataset
//  * @param dim_y       The size of the second dimension of the global dataset
//  * @param local_dim_x The size of first dimension of the local (to each process) dataset
//  */
// void WriteTestDataFourier(fftw_complex* data, char* dset_name, int dim_x, int dim_y, int local_dim_x) {
	
// 	// Initialzie variables
// 	hsize_t rank = 2;
//     hsize_t Dims[rank];
//     herr_t status;

//     // Allocate array for gathering full dataset on root process
// 	fftw_complex* full_data = (fftw_complex* )fftw_malloc(sizeof(fftw_complex) * dim_x * dim_y);
// 	if (full_data == NULL) {
// 		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to allocate memory for writing ["CYAN"%s"RESET"] to test data file\n-->> Exiting!!!\n", dset_name);
// 		exit(1);
// 	}	

// 	// Gather full dataset on root process
// 	MPI_Gather(data, local_dim_x * dim_y, MPI_C_DOUBLE_COMPLEX, full_data, local_dim_x * dim_y, MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	
// 	// Write dataset to file
// 	if (!sys_vars->rank) {   
// 	    Dims[0] = dim_x;
// 	    Dims[1] = dim_y;
// 	    status = H5LTmake_dataset(file_info->test_file_handle, dset_name, rank, Dims, file_info->COMPLEX_DTYPE, full_data);
// 	    if (status < 0) {
// 	    	fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to write ["CYAN"%s"RESET"] to test data file\n-->> Exiting!!!\n", dset_name);
// 	    	exit(1);
// 	    }	
// 	}

// 	// Free memory
// 	fftw_free(full_data);
// }
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------