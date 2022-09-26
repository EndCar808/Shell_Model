/**
* @file hdf5_funcs.c  
* @author Enda Carroll
* @date Sept 2022
* @brief File containing HDF5 function wrappers for creating, opening, wrtining to and closing output file
*/
// ---------------------------------------------------------------------
//  Standard Libraries and Headers
// ---------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
// ---------------------------------------------------------------------
//  User Libraries and Headers
// ---------------------------------------------------------------------
#include "data_types.h"
#include "hdf5_funcs.h"
#include "utils.h"

// ---------------------------------------------------------------------
//  Function Definitions
// ---------------------------------------------------------------------
/**
 * Wrapper function that creates the ouput directory creates and opens the main output file using parallel access and the spectra file using normal serial access
 */
void CreateOutputFilesWriteICs(const long int N) {

	// Initialize variables
	herr_t status;

	// Create compound datatype for the complex datasets
	file_info->COMPLEX_DTYPE = CreateComplexDatatype();


	///////////////////////////
	/// Create & Open Files
	///////////////////////////
	// -----------------------------------
	// Create Output Directory and Path
	// -----------------------------------
	GetOutputDirPath();

	// ---------------------------------
	// Create the output files
	// ---------------------------------
	// Create the main output file
	file_info->output_file_handle = H5Fcreate(file_info->output_file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (file_info->output_file_handle < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"]  --- Could not create main HDF5 output file at: "CYAN"%s"RESET" \n-->>Exiting....\n", file_info->output_file_name);
		exit(1);
	}

	////////////////////////////////
	/// Create Slabbed Datasets
	////////////////////////////////
	// create hdf5 dimension arrays for creating the hyperslabs
	const int Dim = 2;
	hsize_t dims[Dim];      // array to hold dims of full evolution data
	hsize_t maxdims[Dim];   // array to hold max dims of full evolution data
	hsize_t chunkdims[Dim]; // array to hold dims of the hyperslab chunks

	// Initialize the hyperslab arrays
	dims[0]      = sys_vars->num_print_steps;    // Number of timesteps
	dims[1]      = N;              				 // Number of shells
	maxdims[0]   = H5S_UNLIMITED;                // Setting max time index to unlimited means we must chunk our data
	maxdims[1]   = N;                            // Same as before = number of shells
	chunkdims[0] = 1;                            // 1D chunk to be saved 
	chunkdims[1] = N;                            // 1D chunk of size number of shells

	///--------------------------------------- Velocity Modes
	#if defined(__VEL) && !(defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
	CreateSlabbedDSet(0.0, 0, "VelModes", &(file_info->file_space[DSET_VEL]), &(file_info->data_set[DSET_VEL]), &(file_info->mem_space[DSET_VEL]), file_info->COMPLEX_DTYPE, dims, maxdims, chunkdims, Dim);
	#endif

	///--------------------------------------- Velocity Amplitudes
	#if defined(__VEL_AMP) && (defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
	CreateSlabbedDSet(0.0, 0, "VelAmps", &(file_info->file_space[DSET_VEL_AMP]), &(file_info->data_set[DSET_VEL_AMP]), &(file_info->mem_space[DSET_VEL_AMP]), H5T_NATIVE_DOUBLE, dims, maxdims, chunkdims, Dim);
	#endif
	///--------------------------------------- Velocity Phases
	#if defined(__VEL_PHI) && (defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
	CreateSlabbedDSet(0.0, 0, "VelPhases", &(file_info->file_space[DSET_VEL_PHI]), &(file_info->data_set[DSET_VEL_PHI]), &(file_info->mem_space[DSET_VEL_PHI]), H5T_NATIVE_DOUBLE, dims, maxdims, chunkdims, Dim);
	#endif

	#if defined(__MAGNETO)
	///--------------------------------------- Magnetic Modes
	#if defined(__MAG) && !(defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
	CreateSlabbedDSet(0.0, 0, "MagModes", &(file_info->file_space[DSET_MAG]), &(file_info->data_set[DSET_MAG]), &(file_info->mem_space[DSET_MAG]), file_info->COMPLEX_DTYPE, dims, maxdims, chunkdims, Dim);
	#endif

	///--------------------------------------- Magnetic Amplitudes
	#if defined(__MAG_AMP) && (defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
	CreateSlabbedDSet(0.0, 0, "MagAmps", &(file_info->file_space[DSET_MAG_AMP]), &(file_info->data_set[DSET_MAG_AMP]), &(file_info->mem_space[DSET_MAG_AMP]), H5T_NATIVE_DOUBLE, dims, maxdims, chunkdims, Dim);
	#endif
	///--------------------------------------- Magnetic Phases
	#if defined(__MAG_PSI) && (defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
	CreateSlabbedDSet(0.0, 0, "MagPhases", &(file_info->file_space[DSET_MAG_PSI]), &(file_info->data_set[DSET_MAG_PSI]), &(file_info->mem_space[DSET_MAG_PSI]), H5T_NATIVE_DOUBLE, dims, maxdims, chunkdims, Dim);
	#endif
	#endif

	///--------------------------------------- Energy Flux
	#if defined(__ENRG_FLUX)
	// Create slabbed dataset for the flux
	CreateSlabbedDSet(0.0, 0, "EnergyFlux", &(file_info->file_space[DSET_ENRG_FLUX]), &(file_info->data_set[DSET_ENRG_FLUX]), &(file_info->mem_space[DSET_ENRG_FLUX]), H5T_NATIVE_DOUBLE, dims, maxdims, chunkdims, Dim);
	// Create slabbed dataset for the diss
	CreateSlabbedDSet(0.0, 0, "EnergyDiss", &(file_info->file_space[DSET_ENRG_DISS]), &(file_info->data_set[DSET_ENRG_DISS]), &(file_info->mem_space[DSET_ENRG_DISS]), H5T_NATIVE_DOUBLE, dims, maxdims, chunkdims, Dim);
	#endif

	////////////////////////////////
	/// Write Initial Condtions
	////////////////////////////////
	if (sys_vars->TRANS_ITERS_FLAG != TRANSIENT_ITERS) {
		///--------------------------------------- Velocity Modes
		#if defined(__VEL) && !(defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
		WriteSlabbedDataFourier(0.0, 0, file_info->file_space[DSET_VEL], file_info->data_set[DSET_VEL], file_info->mem_space[DSET_VEL], file_info->COMPLEX_DTYPE, &(run_data->u[2]), "VelModes", N, 0);
		#endif

		///--------------------------------------- Velocity Amplitudes
		#if defined(__VEL_AMP) && (defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_VEL_AMP], file_info->data_set[DSET_VEL_AMP], file_info->mem_space[DSET_VEL_AMP], H5T_NATIVE_DOUBLE, &(run_data->a_n[2]), "VelAmps", N, 0);
		#endif
		///--------------------------------------- Velocity Phases
		#if defined(__VEL_PHI) && (defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT)) 
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_VEL_PHI], file_info->data_set[DSET_VEL_PHI], file_info->mem_space[DSET_VEL_PHI], H5T_NATIVE_DOUBLE, &(run_data->phi_n[2]), "VelPhases", N, 0);
		#endif

		#if defined(__MAGNETO)
		///--------------------------------------- Magnetic Modes
		#if defined(__MAG) && !(defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
		WriteSlabbedDataFourier(0.0, 0, file_info->file_space[DSET_MAG], file_info->data_set[DSET_MAG], file_info->mem_space[DSET_MAG], file_info->COMPLEX_DTYPE, &(run_data->b[2]), "MagModes", N, 0);
		#endif

		///--------------------------------------- Magnetic Amplitudes
		#if defined(__MAG_AMP) && (defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_MAG_AMP], file_info->data_set[DSET_MAG_AMP], file_info->mem_space[DSET_MAG_AMP], H5T_NATIVE_DOUBLE, &(run_data->b_n[2]), "MagAmps", N, 0);
		#endif
		///--------------------------------------- Magnetic Phases
		#if defined(__MAG_PSI) && (defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_MAG_PSI], file_info->data_set[DSET_MAG_PSI], file_info->mem_space[DSET_MAG_PSI], H5T_NATIVE_DOUBLE, &(run_data->psi_n[2]), "MagPhases", N, 0);
		#endif
		#endif

		///--------------------------------------- Energy Flux
		#if defined(__ENRG_FLUX)
		// Create slabbed dataset for the flux
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_ENRG_FLUX], file_info->data_set[DSET_ENRG_FLUX], file_info->mem_space[DSET_ENRG_FLUX], H5T_NATIVE_DOUBLE, run_data->energy_flux, "EnergyFlux", N, 0);
		// Create slabbed dataset for the diss
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_ENRG_DISS], file_info->data_set[DSET_ENRG_DISS], file_info->mem_space[DSET_ENRG_DISS], H5T_NATIVE_DOUBLE, run_data->energy_diss, "EnergyDiss", N, 0);
		#endif
	}

	// ------------------------------------
	// Close Identifiers - also close file
	// ------------------------------------
	status = H5Fclose(file_info->output_file_handle);
	if (status < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to close output file ["CYAN"%s"RESET"] at: Iter = ["CYAN"%d"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", file_info->output_file_name, 0, 0.0);
		exit(1);		
	}
}
/**
 * Function that creates the output file paths and directories
 */
void GetOutputDirPath(void) {

	// Initialize variables
	char sys_type[64];
	char solv_type[64];
	char model_type[64];
	char tmp_path[512];
	char file_data[512];  
	struct stat st = {0};	// this is used to check whether the output directories exist or not.

	// ----------------------------------
	// Check if Provided Directory Exists
	// ----------------------------------
	// Check if output directory exists
	if (stat(file_info->output_dir, &st) == -1) {
		printf("\n["YELLOW"NOTE"RESET"] --- Provided Output directory doesn't exist, now creating it...\n");
		// If not then create it
		if ((mkdir(file_info->output_dir, 0700)) == -1) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to create provided output directory ["CYAN"%s"RESET"]\n-->> Exiting...\n", file_info->output_dir);
			exit(1);
		}
	}

	////////////////////////////////////////////
	// Check if Output File Only is Requested
	////////////////////////////////////////////
	if (file_info->file_only) {
		// Update to screen that file only output option is selected
		printf("\n["YELLOW"NOTE"RESET"] --- File only output option selected...\n");
		
		// ----------------------------------
		// Get Simulation Details
		// ----------------------------------
		#if defined(__MAGNETO)
		sprintf(sys_type, "%s", "MAG_HYDRO");
		#elif !defined(__MAGNETO)
		sprintf(sys_type, "%s", "HYDRO");
		#else
		sprintf(sys_type, "%s", "UKN");
		#endif
		#if defined(INT_FAC_RK4)
		sprintf(solv_type, "%s", "INTFACRK4");
		#elif defined(RK4)
		sprintf(solv_type, "%s", "RK4");
		#elif defined(__RK5)
		sprintf(solv_type, "%s", "RK5");
		#elif defined(__DPRK5)
		sprintf(solv_type, "%s", "DP5");
		#else 
		sprintf(solv_type, "%s", "UKN");
		#endif
		#if defined(PHASE_ONLY)
		sprintf(model_type, "%s", "PO");
		#elif defined(PHASE_ONLY_DIRECT)
		sprintf(model_type, "%s", "PO_D");
		#else
		sprintf(model_type, "%s", "FULL");
		#endif

		// -------------------------------------
		// Get File Label from Simulation Data
		// -------------------------------------
		#if defined(__MAGNETO)
		// Construct file label from simulation data
		sprintf(file_data, "_SIM[%s-%s-%s]_N[%ld]_T[%1.1lf-%g-%1.3lf]_NU[%1.8lf]_ETA[%1.8lf]_ALPHA[%1.3lf]_BETA[%1.3lf]_u0[%s].h5", sys_type, solv_type, model_type, sys_vars->N, sys_vars->t0, sys_vars->dt, sys_vars->T, sys_vars->NU, sys_vars->ETA, sys_vars->ALPHA, sys_vars->BETA, sys_vars->u0);
		#else
		sprintf(file_data, "_SIM[%s-%s-%s]_N[%ld]_T[%1.1lf-%g-%1.3lf]_NU[%1.8lf]_ALPHA[%1.3lf]_BETA[%1.3lf]_u0[%s].h5", sys_type, solv_type, model_type, sys_vars->N, sys_vars->t0, sys_vars->dt, sys_vars->T, sys_vars->NU, sys_vars->ALPHA, sys_vars->BETA, sys_vars->u0);
		#endif

		// ----------------------------------
		// Construct File Paths
		// ---------------------------------- 
		// Construct main file path
		strcpy(tmp_path, file_info->output_dir);
		strcat(tmp_path, "Main_HDF_Data"); 
		strcpy(file_info->output_file_name, tmp_path); 
		strcat(file_info->output_file_name, file_data);
		printf("\nMain Output File: "CYAN"%s"RESET"\n\n", file_info->output_file_name);
	}
	else {
		// ----------------------------------
		// Get Simulation Details
		// ----------------------------------
		#if defined(__MAGNETO)
		sprintf(sys_type, "%s", "MAG_HYDRO");
		#elif !defined(__MAGNETO)
		sprintf(sys_type, "%s", "HYDRO");
		#else
		sprintf(sys_type, "%s", "UKN");
		#endif
		#if defined(INT_FAC_RK4)
		sprintf(solv_type, "%s", "INTFACRK4");
		#elif defined(__RK5)
		sprintf(solv_type, "%s", "RK5");
		#elif defined(RK4)
		sprintf(solv_type, "%s", "RK4");
		#elif defined(__DPRK5)
		sprintf(solv_type, "%s", "DP5");
		#else 
		sprintf(solv_type, "%s", "UKN");
		#endif
		#if defined(PHASE_ONLY)
		sprintf(model_type, "%s", "PO");
		#elif defined(PHASE_ONLY_DIRECT)
		sprintf(model_type, "%s", "PO_D");
		#else
		sprintf(model_type, "%s", "FULL");
		#endif

		// ----------------------------------
		// Construct Output folder
		// ----------------------------------
		#if defined(__MAGNETO)
		// Construct file label from simulation data
		sprintf(file_data, "SIM[%s-%s-%s]_N[%ld]_T[%1.1lf-%g-%1.3lf]_NU[%1.8lf]_ETA[%1.8lf]_ALPHA[%1.3lf]_BETA[%1.3lf]_u0[%s]_TAG[%s]/", sys_type, solv_type, model_type, sys_vars->N, sys_vars->t0, sys_vars->dt, sys_vars->T, sys_vars->NU, sys_vars->ETA, sys_vars->ALPHA, sys_vars->BETA, sys_vars->u0, file_info->output_tag);
		#else
		sprintf(file_data, "SIM[%s-%s-%s]_N[%ld]_T[%1.1lf-%g-%1.3lf]_NU[%1.8lf]_ETA[%1.8lf]_ALPHA[%1.3lf]_BETA[%1.3lf]_u0[%s]_TAG[%s]/", sys_type, solv_type, model_type, sys_vars->N, sys_vars->t0, sys_vars->dt, sys_vars->T, sys_vars->NU, sys_vars->ALPHA, sys_vars->BETA, sys_vars->u0, file_info->output_tag);
		#endif
		
		// ----------------------------------
		// Check Existence of Output Folder
		// ----------------------------------
		strcat(file_info->output_dir, file_data);
		// Check if folder exists
		if (stat(file_info->output_dir, &st) == -1) {
			// If not create it
			if ((mkdir(file_info->output_dir, 0700)) == -1) {
				fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to create folder for output files ["CYAN"%s"RESET"]\n-->> Exiting...\n", file_info->output_dir);
				exit(1);
			}
		}

		// ----------------------------------
		// Construct File Paths
		// ---------------------------------- 
		// Construct main file path
		strcpy(file_info->output_file_name, file_info->output_dir); 
		strcat(file_info->output_file_name, "Main_HDF_Data.h5");
		printf("\nMain Output File: "CYAN"%s"RESET"\n\n", file_info->output_file_name);
	}

}
/**
 * Wrapper function that writes the data to file by openining it, creating a group for the current iteration and writing the data under this group. The file is then closed again 
 * @param t     The current time of the simulation
 * @param dt    The current timestep being used
 * @param iters The current iteration
 */
void WriteDataToFile(double t, double dt, long int iters) {

	// Initialize variables
	herr_t status;
	static const int d_set_rank2D = 2;
	hsize_t dset_dims2D[d_set_rank2D];        // array to hold dims of the dataset to be created
	hsize_t slab_dims2D[d_set_rank2D];	      // Array to hold the dimensions of the hyperslab
	hsize_t mem_space_dims2D[d_set_rank2D];   // Array to hold the dimensions of the memoray space - for real data this will be different to slab_dims due to 0 padding


	// --------------------------------------
	// Check if files exist and Open/Create
	// --------------------------------------
	// /// Check if main file exists - open it if it does if not create it
	if (access(file_info->output_file_name, F_OK) != 0) {
		file_info->output_file_handle = H5Fcreate(file_info->output_file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		if (file_info->output_file_handle < 0) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to create output file ["CYAN"%s"RESET"] at: Iter = ["CYAN"%ld"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", file_info->output_file_name, iters, t);
			exit(1);
		}
	}
	else {
		// Open file with parallel I/O access properties
		file_info->output_file_handle = H5Fopen(file_info->output_file_name, H5F_ACC_RDWR, H5P_DEFAULT);
		if (file_info->output_file_handle < 0) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to open output file ["CYAN"%s"RESET"] at: Iter = ["CYAN"%ld"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", file_info->output_file_name, iters, t);
			exit(1);
		}
	}

	// -------------------------------
	// Write Data to File
	// -------------------------------
	///--------------------------------------- Velocity Modes
	#if defined(__VEL) && !(defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT)) 
	WriteSlabbedDataFourier(t, iters, file_info->file_space[DSET_VEL], file_info->data_set[DSET_VEL], file_info->mem_space[DSET_VEL], file_info->COMPLEX_DTYPE, &(run_data->u[2]), "VelModes", sys_vars->N, iters);
	#endif

	///--------------------------------------- Velocity Amplitudes
	#if defined(__VEL_AMP) && (defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_VEL_AMP], file_info->data_set[DSET_VEL_AMP], file_info->mem_space[DSET_VEL_AMP], H5T_NATIVE_DOUBLE, &(run_data->a_n[2]), "VelAmps", sys_vars->N, iters);
	#endif
	///--------------------------------------- Velocity Phases
	#if defined(__VEL_PHI) && (defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_VEL_PHI], file_info->data_set[DSET_VEL_PHI], file_info->mem_space[DSET_VEL_PHI], H5T_NATIVE_DOUBLE, &(run_data->phi_n[2]), "VelPhases", sys_vars->N, iters);
	#endif

	#if defined(__MAGNETO)
	///--------------------------------------- Magnetic Modes
	#if defined(__MAG) && !(defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
	WriteSlabbedDataFourier(t, iters, file_info->file_space[DSET_MAG], file_info->data_set[DSET_MAG], file_info->mem_space[DSET_MAG], file_info->COMPLEX_DTYPE, &(run_data->b[2]), "MagModes", sys_vars->N, iters);
	#endif

	///--------------------------------------- Magnetic Amplitudes
	#if defined(__MAG_AMP) && (defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_MAG_AMP], file_info->data_set[DSET_MAG_AMP], file_info->mem_space[DSET_MAG_AMP], H5T_NATIVE_DOUBLE, &(run_data->b_n[2]), "MagAmps", sys_vars->N, iters);
	#endif
	///--------------------------------------- Magnetic Phases
	#if defined(__MAG_PSI) && (defined(PHASE_ONLY) || defined(PHASE_ONLY_DIRECT))
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_MAG_PSI], file_info->data_set[DSET_MAG_PSI], file_info->mem_space[DSET_MAG_PSI], H5T_NATIVE_DOUBLE, &(run_data->psi_n[2]), "MagPhases", sys_vars->N, iters);
	#endif
	#endif

	///--------------------------------------- Energy Flux
	#if defined(__ENRG_FLUX)
	// Create slabbed dataset for the flux
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_ENRG_FLUX], file_info->data_set[DSET_ENRG_FLUX], file_info->mem_space[DSET_ENRG_FLUX], H5T_NATIVE_DOUBLE, run_data->energy_flux, "EnergyFlux", sys_vars->N, iters);
	// Create slabbed dataset for the diss
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_ENRG_DISS], file_info->data_set[DSET_ENRG_DISS], file_info->mem_space[DSET_ENRG_DISS], H5T_NATIVE_DOUBLE, run_data->energy_diss, "EnergyDiss", sys_vars->N, iters);
	#endif

	// -------------------------------
	// Close identifiers and File
	// -------------------------------
	status = H5Fclose(file_info->output_file_handle);
	if (status < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to close output file ["CYAN"%s"RESET"] at: Iter = ["CYAN"%ld"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", file_info->output_file_name, iters, t);
		exit(1);
	}
}
/**
 * Wrapper function for creating a slabbed dataset in the ouput file
 * @param t     		  Current time of the simulation
 * @param iters       	  Current iteration
 * @param dset_name       Name of the dataset
 * @param file_space      Identifier for the filespace
 * @param data_set        Identifieer for the dataset in the file
 * @param mem_space       Identifier for the memory space
 * @param dtype           The type of data for the dataset
 * @param dset_dims       The size of the dimensions of the dataset
 * @param dset_max_dims   The maximum allowed dimensions of the dataset
 * @param dset_chunk_dims The size of the chuncks to be used
 * @param num_dims        The dimension of the dataset
 */
void CreateSlabbedDSet(double t, int iters, char* dset_name, hid_t* file_space, hid_t* data_set, hid_t* mem_space, hid_t dtype, hsize_t* dset_dims, hsize_t* dset_max_dims, hsize_t* dset_chunk_dims, const int num_dims) {

	// Initialize variables
	herr_t status;

	// ------------------------------------
	// Create File Space
	// ------------------------------------
	// Create the 2D dataspace - setting the no. of dimensions, expected and max size of the dimensions
	*file_space = H5Screate_simple(num_dims, dset_dims, dset_max_dims);
	if (*file_space < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to set file space for dataset ["CYAN"%s"RESET"] at: Iter = ["CYAN"%d"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", dset_name, iters, t);
		exit(1);
	}

	// Must create a propertly list to enable data chunking due to max dimension being unlimited
	// Create property list 
	hid_t plist = H5Pcreate(H5P_DATASET_CREATE);

	// Using this property list set the chuncking - stores the chunking info in plist
	status = H5Pset_chunk(plist, num_dims, dset_chunk_dims);
	if (status < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to create property list for dataset ["CYAN"%s"RESET"] at: Iter = ["CYAN"%d"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", dset_name, iters, t);
		exit(1);
	}

	// ------------------------------------
	// Create Dataset
	// ------------------------------------
	// Create the dataset in the previouosly created datafile - using the chunk enabled property list and new compound datatype
	*data_set = H5Dcreate(file_info->output_file_handle, dset_name, dtype, *file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
	if (*data_set < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to create dataset space for dataset ["CYAN"%s"RESET"] at: Iter = ["CYAN"%d"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", dset_name, iters, t);
		exit(1);
	}
	
	// Create the memory space for the slab
	// setting the max dims to NULL defaults to same size as dims
	*mem_space = H5Screate_simple(num_dims, dset_chunk_dims, NULL);
	if (*mem_space < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to create memory space for dataset ["CYAN"%s"RESET"] at: Iter = ["CYAN"%d"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", dset_name, iters, t);
		exit(1);
	}

	// Close the property list object
	status = H5Pclose(plist);
	if (status < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to close property for dataset ["CYAN"%s"RESET"] at: Iter = ["CYAN"%d"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", dset_name, iters, t);
		exit(1);
	}
}
/**
 * Wrapper function to select appropriate hyperslab and write this to file
 * @param t          Current time in the simulation
 * @param iters      Current iteration
 * @param file_space Identifier for the file space
 * @param data_set   Identifier for the Dataset
 * @param mem_space  Identifier for the memory
 * @param dtype      The data type to be written
 * @param data       The data
 * @param data_name  The name of the dataset
 * @param n          The size of the chunk to write
 * @param index      The index in the dataset to write chunk to
 */
void WriteSlabbedDataReal(double t, int iters, hid_t file_space, hid_t data_set, hid_t mem_space, hid_t dtype, double* data, char* dset_name, int n, int index) {

	// Initialize variables
	hsize_t start_index[2]; // stores the index in the hyperslabbed dataset to start writing to
	hsize_t count[2];       // stores the size of hyperslab to write to the dataset

	// ------------------------------------
	// Get Chunk Dimensions
	// ------------------------------------
	count[0]       = 1;		// 1D slab so first dim is 1
	count[1]       = n;		// 1D slab of size of data array
	start_index[0] = index;	// set the starting row index to index in the global dataset to write slab to
	start_index[1] = 0;		// set column index to 0 to start writing from the first column

	// ------------------------------------
	// Select Slab
	// ------------------------------------
	// Select appropriate hyperslab 
	if ((H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start_index, NULL, count, NULL)) < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to select hyperslab in file for datset ["CYAN"%s"RESET"] at: Iter = ["CYAN"%d"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", dset_name, iters, t);
		exit(1);
	}

	// ------------------------------------
	// Write Chunk to File
	// ------------------------------------
	// Then write the current chunk to this hyperslab
	if ((H5Dwrite(data_set, dtype, mem_space, file_space, H5P_DEFAULT, data)) < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to write data to datset ["CYAN"%s"RESET"] at: Iter = ["CYAN"%d"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", dset_name, iters, t);
		exit(1);
	}
}
/**
 * Wrapper function to select appropriate hyperslab and write this to file
 * @param t          Current time in the simulation
 * @param iters      Current iteration
 * @param file_space Identifier for the file space
 * @param data_set   Identifier for the Dataset
 * @param mem_space  Identifier for the memory
 * @param dtype      The data type to be written
 * @param data       The data
 * @param data_name  The name of the dataset
 * @param n          The size of the chunk to write
 * @param index      The index in the dataset to write chunk to
 */
void WriteSlabbedDataFourier(double t, int iters, hid_t file_space, hid_t data_set, hid_t mem_space, hid_t dtype, fftw_complex* data, char* dset_name, int n, int index) {

	// Initialize variables
	hsize_t start_index[2]; // stores the index in the hyperslabbed dataset to start writing to
	hsize_t count[2];       // stores the size of hyperslab to write to the dataset

	// ------------------------------------
	// Get Chunk Dimensions
	// ------------------------------------
	count[0]       = 1;		// 1D slab so first dim is 1
	count[1]       = n;		// 1D slab of size of data array
	start_index[0] = index;	// set the starting row index to index in the global dataset to write slab to
	start_index[1] = 0;		// set column index to 0 to start writing from the first column

	// ------------------------------------
	// Select Slab
	// ------------------------------------
	// Select appropriate hyperslab 
	if ((H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start_index, NULL, count, NULL)) < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to select hyperslab in file for datset ["CYAN"%s"RESET"] at: Iter = ["CYAN"%d"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", dset_name, iters, t);
		exit(1);
	}

	// ------------------------------------
	// Write Chunk to File
	// ------------------------------------
	// Then write the current chunk to this hyperslab
	if ((H5Dwrite(data_set, dtype, mem_space, file_space, H5P_DEFAULT, data)) < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to write data to datset ["CYAN"%s"RESET"] at: Iter = ["CYAN"%d"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", dset_name, iters, t);
		exit(1);
	}
}
/**
 * Wrapper function that writes all the non-slabbed/chunk datasets to file after integeration has finished 
 * @param N 			 Array containing the dimensions of the system
 * @param iters 	     The number of iterations performed by the simulation 
 * @param save_data_indx The number of saving steps performed by the simulation
 */
void FinalWriteAndCloseOutputFile(const long int N, int iters, int save_data_indx) {

	// Initialize variables
	herr_t status;
	static const hsize_t D1 = 1;
	hsize_t dims1D[D1];

	// Record total iterations
	sys_vars->tot_iters      = (long int)iters - 1;
	sys_vars->tot_save_steps = (long int)save_data_indx - 1;

	////////////////////////////////
	/// Repon and Write Datasets
	////////////////////////////////
	// Repon Output file with read/write permissions
	file_info->output_file_handle = H5Fopen(file_info->output_file_name, H5F_ACC_RDWR , H5P_DEFAULT);
	if (file_info->output_file_handle < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to reopen output file for writing non chunked/slabbed datasets! \n-->>Exiting....\n");
		exit(1);
	}

	// -------------------------------
	// Write Wavenumbers
	// -------------------------------
	#if defined(__WAVELIST)
	// Allocate array to gather the wavenumbers from each of the local arrays - in the x direction
	dims1D[0] = N;
	if ( (H5LTmake_dataset(file_info->output_file_handle, "k", D1, dims1D, H5T_NATIVE_INT, run_data->k)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "k");
	}
	#endif

	// -------------------------------
	// Write System Measures
	// -------------------------------
	// Time
	#if defined(__TIME)
	// Time array only on rank 0
	dims1D[0] = sys_vars->num_print_steps;
	if ( (H5LTmake_dataset(file_info->output_file_handle, "Time", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->time)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "Time");
	}
	#endif

	#if defined(__SYS_MEASURES)
	// Energy
	if ( (H5LTmake_dataset(file_info->output_file_handle, "TotalEnergy", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->tot_energy)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TotalEnergy");
	}
	#if defined(__MAGNETO)
	// Helicity
	if ( (H5LTmake_dataset(file_info->output_file_handle, "TotalHelicity", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->tot_hel)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TotalEnstrophy");
	}
	// Cross Helicity
	if ( (H5LTmake_dataset(file_info->output_file_handle, "TotalCrossHelicity", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->tot_cross_hel)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TotalPalinstrophy");
	}
	#endif
	#endif

	// -----------------------------------
	// Close Files for the final time
	// -----------------------------------
	status = H5Fclose(file_info->output_file_handle);
	if (status < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to close main output file: "CYAN"%s"RESET" \n-->> Exiting....\n", file_info->output_file_name);
		exit(1);
	}

	#if defined(__VORT_FOUR) || defined(__MODES)
	// Close the complex datatype identifier
	H5Tclose(file_info->COMPLEX_DTYPE);
	#endif
}
/**
 * Function to create a HDF5 datatype for complex data
 */
hid_t CreateComplexDatatype(void) {

	// Declare HDF5 datatype variable
	hid_t dtype;

	// error handling var
	herr_t status;
	
	// Create complex struct
	struct complex_type_tmp cmplex;
	cmplex.re = 0.0;
	cmplex.im = 0.0;

	// create complex compound datatype
	dtype  = H5Tcreate(H5T_COMPOUND, sizeof(cmplex));

	// Insert the real part of the datatype
  	status = H5Tinsert(dtype, "r", offsetof(complex_type_tmp,re), H5T_NATIVE_DOUBLE);
  	if (status < 0) {
  		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Could not insert real part for the Complex Compound Datatype!!\n-->> Exiting...\n");
  		exit(1);
  	}

  	// Insert the imaginary part of the datatype
  	status = H5Tinsert(dtype, "i", offsetof(complex_type_tmp,im), H5T_NATIVE_DOUBLE);
  	if (status < 0) {
  		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Could not insert imaginary part for the Complex Compound Datatype!!\n-->> Exiting...\n");
  		exit(1);
  	}

  	return dtype;
}
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------