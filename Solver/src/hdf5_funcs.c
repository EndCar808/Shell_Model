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
#include "stats.h"
#include "phase_sync.h"
#include "sys_msr.h"

// ---------------------------------------------------------------------
//  Function Definitions
// ---------------------------------------------------------------------
/**
 * Wrapper function that creates the ouput directory creates and opens the main output file using parallel access and the spectra file using normal serial access
 */
void CreateOutputFilesWriteICs(const long int N) {

	// Initialize variables
	herr_t status;

	#if defined(__VEL) || defined(__MAG) || defined(__Z_PLUS) || defined(__Z_MINUS) || defined(__FORCING) || defined(__STATS) || defined(__PHASE_SYNC)
	if (sys_vars->INPUT_FILE_FLAG == NO_INPUT_FILE) {
		// Create compound datatype for the complex datasets if it has not been created yet
		file_info->COMPLEX_DTYPE = CreateComplexDatatype();
	}
	#endif


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
	#if defined(__CONSERVED_PHASES) || defined(__PHASE_SYNC) || defined(__PHASE_SYNC)
	file_info->phase_sync_file_handle = H5Fcreate(file_info->phase_sync_file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (file_info->phase_sync_file_handle < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"]  --- Could not create main Phase Sync output file at: "CYAN"%s"RESET" \n-->>Exiting....\n", file_info->phase_sync_file_name);
		exit(1);
	}
	#endif

	////////////////////////////////
	/// Create Slabbed Datasets
	////////////////////////////////
	// create hdf5 dimension arrays for creating the hyperslabs
	const int Dim2D = 2;
	hsize_t dims2D[Dim2D];      	// array to hold dims of full evolution data
	hsize_t maxdims2D[Dim2D];   	// array to hold max dims of full evolution data
	hsize_t chunkdims2D[Dim2D]; 	// array to hold dims of the hyperslab chunks
	hsize_t index2D[Dim2D]; 		// stores the index in the hyperslabbed dataset to start writing to
	hsize_t count2D[Dim2D];       	// stores the size of hyperslab to write to the dataset
	#if (defined(__MAGNETO) || defined(__ELSASSAR_MHD)) && (defined(__CONSERVED_PHASES) || defined(__PHASE_SYNC))
	const int Dim3D = 3;
	hsize_t dims3D[Dim3D];      	// array to hold dims of full evolution data
	hsize_t maxdims3D[Dim3D];   	// array to hold max dims of full evolution data
	hsize_t chunkdims3D[Dim3D]; 	// array to hold dims of the hyperslab chunks
	hsize_t index3D[Dim3D]; 		// stores the index in the hyperslabbed dataset to start writing to
	hsize_t count3D[Dim3D];       	// stores the size of hyperslab to write to the dataset
	#endif

	// Initialize the hyperslab arrays
	dims2D[0]      = sys_vars->num_print_steps;  	// Number of timesteps
	dims2D[1]      = N;          				 	// Number of shells + 0 mode
	maxdims2D[0]   = H5S_UNLIMITED;              	// Setting max time index to unlimited means we must chunk our data
	maxdims2D[1]   = N;                        		// Same as before = number of shells + 0 mode
	chunkdims2D[0] = 1;                          	// 1D chunk to be saved 
	chunkdims2D[1] = N;                        		// 1D chunk of size number of shells + 0 mode
	count2D[0]     = 1;								// 1D slab so first dim is 1
	count2D[1]     = N;								// 1D slab of size of data array
	index2D[0]     = 0;								// set the starting row index to index in the global dataset to write slab to -> 0th row for 0 iter
	index2D[1]     = 0;								// set column index to 0 to start writing from the first column

	///--------------------------------------- Velocity Modes
	#if defined(__VEL) && !(defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "VelModes", &(file_info->file_space[DSET_VEL]), &(file_info->data_set[DSET_VEL]), &(file_info->mem_space[DSET_VEL]), file_info->COMPLEX_DTYPE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	#endif

	///--------------------------------------- Velocity Amplitudes
	#if defined(__VEL_AMP) && (defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "VelAmps", &(file_info->file_space[DSET_VEL_AMP]), &(file_info->data_set[DSET_VEL_AMP]), &(file_info->mem_space[DSET_VEL_AMP]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	#endif
	///--------------------------------------- Velocity Phases
	#if defined(__VEL_PHI) && (defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "VelPhases", &(file_info->file_space[DSET_VEL_PHI]), &(file_info->data_set[DSET_VEL_PHI]), &(file_info->mem_space[DSET_VEL_PHI]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	#endif

	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	///--------------------------------------- Magnetic Modes
	#if defined(__MAG) && !(defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "MagModes", &(file_info->file_space[DSET_MAG]), &(file_info->data_set[DSET_MAG]), &(file_info->mem_space[DSET_MAG]), file_info->COMPLEX_DTYPE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	#endif

	#if defined(__ELSASSAR_MHD)
	///--------------------------------------- Z_Plus Elsassar Variable
	#if defined(__Z_PLUS) && !(defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "ZPlus", &(file_info->file_space[DSET_Z_PLUS]), &(file_info->data_set[DSET_Z_PLUS]), &(file_info->mem_space[DSET_Z_PLUS]), file_info->COMPLEX_DTYPE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	#endif
	///--------------------------------------- Z_Minus Elsassar Variable
	#if defined(__Z_MINUS) && !(defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "ZMinus", &(file_info->file_space[DSET_Z_MINUS]), &(file_info->data_set[DSET_Z_MINUS]), &(file_info->mem_space[DSET_Z_MINUS]), file_info->COMPLEX_DTYPE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	#endif
	#endif

	///--------------------------------------- Magnetic Amplitudes
	#if defined(__MAG_AMP) && (defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "MagAmps", &(file_info->file_space[DSET_MAG_AMP]), &(file_info->data_set[DSET_MAG_AMP]), &(file_info->mem_space[DSET_MAG_AMP]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	#endif
	///--------------------------------------- Magnetic Phases
	#if defined(__MAG_PSI) && (defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "MagPhases", &(file_info->file_space[DSET_MAG_PSI]), &(file_info->data_set[DSET_MAG_PSI]), &(file_info->mem_space[DSET_MAG_PSI]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	#endif
	#endif

	///--------------------------------------- Forcing
	#if defined(__FORCING)
	if(!(strcmp(sys_vars->forcing, "STOC"))) {
		// Velocity Forcing
		CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "VelocityForcingInTime", &(file_info->file_space[DSET_FORCING_U]), &(file_info->data_set[DSET_FORCING_U]), &(file_info->mem_space[DSET_FORCING_U]), file_info->COMPLEX_DTYPE, dims2D, maxdims2D, chunkdims2D, Dim2D);	
		// Magnetic Forcing
		CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "MagneticForcingInTime", &(file_info->file_space[DSET_FORCING_B]), &(file_info->data_set[DSET_FORCING_B]), &(file_info->mem_space[DSET_FORCING_B]), file_info->COMPLEX_DTYPE, dims2D, maxdims2D, chunkdims2D, Dim2D);	
	}
	#endif

	///--------------------------------------- Spectra
	#if defined(__ENRG_SPECT)
	// Total Energy Spectrum
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "EnergySpectrum", &(file_info->file_space[DSET_ENRG_SPECT]), &(file_info->data_set[DSET_ENRG_SPECT]), &(file_info->mem_space[DSET_ENRG_SPECT]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);	
	#endif
	#if defined(__DISS_SPECT)
	// Dissipation Spectrum
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "DissipationSpectrum", &(file_info->file_space[DSET_DISS_SPECT]), &(file_info->data_set[DSET_DISS_SPECT]), &(file_info->mem_space[DSET_DISS_SPECT]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	#endif
	#if defined(__KIN_ENRG_SPECT)
	// Kinetic Energy Spectrum
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "KineticEnergySpectrum", &(file_info->file_space[DSET_KIN_ENRG_SPECT]), &(file_info->data_set[DSET_KIN_ENRG_SPECT]), &(file_info->mem_space[DSET_KIN_ENRG_SPECT]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);	
	#endif
	#if defined(__MAG_ENRG_SPECT) && (defined(__MAGNETO) || defined(__ELSASSAR_MHD))
	// Magnetic Energy Spectrum
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "MagneticEnergySpectrum", &(file_info->file_space[DSET_MAG_ENRG_SPECT]), &(file_info->data_set[DSET_MAG_ENRG_SPECT]), &(file_info->mem_space[DSET_MAG_ENRG_SPECT]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);	
	#endif

	///--------------------------------------- Energy Flux
	#if defined(__ENRG_FLUX)
	// Create slabbed dataset for the flux
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "EnergyFlux", &(file_info->file_space[DSET_ENRG_FLUX]), &(file_info->data_set[DSET_ENRG_FLUX]), &(file_info->mem_space[DSET_ENRG_FLUX]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	// Create slabbed dataset for the diss for the velocity field
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "VelEnergyDiss", &(file_info->file_space[DSET_ENRG_DISS_VEL]), &(file_info->data_set[DSET_ENRG_DISS_VEL]), &(file_info->mem_space[DSET_ENRG_DISS_VEL]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	// Create slabbed dataset for the input for the velocity field
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "VelEnergyInput", &(file_info->file_space[DSET_ENRG_INPT_VEL]), &(file_info->data_set[DSET_ENRG_INPT_VEL]), &(file_info->mem_space[DSET_ENRG_INPT_VEL]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	#if defined(__MAGNETO)
	// Create slabbed dataset for the diss for the magnetic field
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "MagEnergyDiss", &(file_info->file_space[DSET_ENRG_DISS_MAG]), &(file_info->data_set[DSET_ENRG_DISS_MAG]), &(file_info->mem_space[DSET_ENRG_DISS_MAG]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	// Create slabbed dataset for the input for the magnetic field
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "MagEnergyInput", &(file_info->file_space[DSET_ENRG_INPT_MAG]), &(file_info->data_set[DSET_ENRG_INPT_MAG]), &(file_info->mem_space[DSET_ENRG_INPT_MAG]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	#endif
	#endif

	///--------------------------------------- Triads and Phase Differences
	#if defined(__CONSERVED_PHASES)
	// Create slabbed dataset for the velocity triads
	dims2D[1]      = phase_sync->num_triads;
	maxdims2D[1]   = phase_sync->num_triads;
	chunkdims2D[1] = phase_sync->num_triads;
	CreateSlabbedDSet(0.0, 0, file_info->phase_sync_file_handle, "VelTriads", &(file_info->file_space[DSET_VEL_TRIADS]), &(file_info->data_set[DSET_VEL_TRIADS]), &(file_info->mem_space[DSET_VEL_TRIADS]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	
	// Create slabbed dataset for the velocity phase differences
	dims2D[1]      = phase_sync->num_phase_diff;
	maxdims2D[1]   = phase_sync->num_phase_diff;
	chunkdims2D[1] = phase_sync->num_phase_diff;
	CreateSlabbedDSet(0.0, 0, file_info->phase_sync_file_handle, "VelPhaseDifferences", &(file_info->file_space[DSET_VEL_PHASE_DIFF]), &(file_info->data_set[DSET_VEL_PHASE_DIFF]), &(file_info->mem_space[DSET_VEL_PHASE_DIFF]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	// Create slabbed dataset for the magnetic triads
	dims3D[0]      = sys_vars->num_print_steps;
	dims3D[1]      = NUM_MAG_TRIAD_TYPES;
	dims3D[2]      = phase_sync->num_triads;
	maxdims3D[0]   = H5S_UNLIMITED;
	maxdims3D[1]   = NUM_MAG_TRIAD_TYPES;
	maxdims3D[2]   = phase_sync->num_triads;
	chunkdims3D[0] = 1;                          
	chunkdims3D[1] = NUM_MAG_TRIAD_TYPES;
	chunkdims3D[2] = phase_sync->num_triads;
	CreateSlabbedDSet(0.0, 0, file_info->phase_sync_file_handle, "MagTriads", &(file_info->file_space[DSET_MAG_TRIADS]), &(file_info->data_set[DSET_MAG_TRIADS]), &(file_info->mem_space[DSET_MAG_TRIADS]), H5T_NATIVE_DOUBLE, dims3D, maxdims3D, chunkdims3D, Dim3D);

	// Create slabbed dataset for the magnetic phase differences
	dims2D[1]      = phase_sync->num_phase_diff;
	maxdims2D[1]   = phase_sync->num_phase_diff;
	chunkdims2D[1] = phase_sync->num_phase_diff;
	CreateSlabbedDSet(0.0, 0, file_info->phase_sync_file_handle, "MagPhaseDifferences", &(file_info->file_space[DSET_MAG_PHASE_DIFF]), &(file_info->data_set[DSET_MAG_PHASE_DIFF]), &(file_info->mem_space[DSET_MAG_PHASE_DIFF]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	#endif
	#endif

	///--------------------------------------- Phase Order Parameters
	#if defined(__PHASE_SYNC)
	// Create slabbed dataset for the velocity triad order parameter
	dims2D[1]      = phase_sync->num_triads;
	maxdims2D[1]   = phase_sync->num_triads;
	chunkdims2D[1] = phase_sync->num_triads;
	CreateSlabbedDSet(0.0, 0, file_info->phase_sync_file_handle, "VelTriadOrderParameter", &(file_info->file_space[DSET_VEL_TRIAD_ORDER]), &(file_info->data_set[DSET_VEL_TRIAD_ORDER]), &(file_info->mem_space[DSET_VEL_TRIAD_ORDER]), file_info->COMPLEX_DTYPE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	// Create slabbed dataset for the velocity phase difference order parameter
	dims2D[1]      = phase_sync->num_phase_diff;
	maxdims2D[1]   = phase_sync->num_phase_diff;
	chunkdims2D[1] = phase_sync->num_phase_diff;
	CreateSlabbedDSet(0.0, 0, file_info->phase_sync_file_handle, "VelPhaseDifferenceOrderParameter", &(file_info->file_space[DSET_VEL_PHASE_DIFF_ORDER]), &(file_info->data_set[DSET_VEL_PHASE_DIFF_ORDER]), &(file_info->mem_space[DSET_VEL_PHASE_DIFF_ORDER]), file_info->COMPLEX_DTYPE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	// Initialize the hyperslab arrays for the magnetic triad order parameter
	dims3D[0]      = sys_vars->num_print_steps;
	dims3D[1]      = NUM_MAG_TRIAD_TYPES;
	dims3D[2]      = phase_sync->num_triads;
	maxdims3D[0]   = H5S_UNLIMITED;
	maxdims3D[1]   = NUM_MAG_TRIAD_TYPES;
	maxdims3D[2]   = phase_sync->num_triads;
	chunkdims3D[0] = 1;                          
	chunkdims3D[1] = NUM_MAG_TRIAD_TYPES;
	chunkdims3D[2] = phase_sync->num_triads;
	CreateSlabbedDSet(0.0, 0, file_info->phase_sync_file_handle, "MagTriadOrderParameter", &(file_info->file_space[DSET_MAG_TRIAD_ORDER]), &(file_info->data_set[DSET_MAG_TRIAD_ORDER]), &(file_info->mem_space[DSET_MAG_TRIAD_ORDER]), file_info->COMPLEX_DTYPE, dims3D, maxdims3D, chunkdims3D, Dim3D);
	dims2D[1]      = phase_sync->num_phase_diff;
	maxdims2D[1]   = phase_sync->num_phase_diff;
	chunkdims2D[1] = phase_sync->num_phase_diff;
	// Create slabbed dataset for the magnetic phase difference order parameter
	CreateSlabbedDSet(0.0, 0, file_info->phase_sync_file_handle, "MagPhaseDifferenceOrderParameter", &(file_info->file_space[DSET_MAG_PHASE_DIFF_ORDER]), &(file_info->data_set[DSET_MAG_PHASE_DIFF_ORDER]), &(file_info->mem_space[DSET_MAG_PHASE_DIFF_ORDER]), file_info->COMPLEX_DTYPE, dims2D, maxdims2D, chunkdims2D, Dim2D);
	#endif
	#endif



	////////////////////////////////
	/// Write Initial Condtions
	////////////////////////////////
	if (sys_vars->TRANS_ITERS_FLAG != TRANSIENT_ITERS) {
		///--------------------------------------- Velocity Modes
		#if defined(__VEL) && !(defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
		WriteSlabbedDataFourier(0.0, 0, file_info->file_space[DSET_VEL], file_info->data_set[DSET_VEL], file_info->mem_space[DSET_VEL], file_info->COMPLEX_DTYPE, &(run_data->u[2]), "VelModes", count2D, index2D);
		#endif

		///--------------------------------------- Velocity Amplitudes
		#if defined(__VEL_AMP) && (defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_VEL_AMP], file_info->data_set[DSET_VEL_AMP], file_info->mem_space[DSET_VEL_AMP], H5T_NATIVE_DOUBLE, &(run_data->a_n[2]), "VelAmps", count2D, index2D);
		#endif
		///--------------------------------------- Velocity Phases
		#if defined(__VEL_PHI) && (defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY)) 
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_VEL_PHI], file_info->data_set[DSET_VEL_PHI], file_info->mem_space[DSET_VEL_PHI], H5T_NATIVE_DOUBLE, &(run_data->phi_n[2]), "VelPhases", count2D, index2D);
		#endif

		#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
		///--------------------------------------- Magnetic Modes
		#if defined(__MAG) && !(defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
		WriteSlabbedDataFourier(0.0, 0, file_info->file_space[DSET_MAG], file_info->data_set[DSET_MAG], file_info->mem_space[DSET_MAG], file_info->COMPLEX_DTYPE, &(run_data->b[2]), "MagModes", count2D, index2D);
		#endif

		#if defined(__ELSASSAR_MHD)
		///--------------------------------------- Z_Plus Elsassar Variable
		#if defined(__Z_PLUS) && !(defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
		WriteSlabbedDataFourier(0.0, 0, file_info->file_space[DSET_Z_PLUS], file_info->data_set[DSET_Z_PLUS], file_info->mem_space[DSET_Z_PLUS], file_info->COMPLEX_DTYPE, &(run_data->z_plus[2]), "ZPlus", count2D, index2D);
		#endif
		///--------------------------------------- Z_Minus Elsassar Variables
		#if defined(__Z_MINUS) && !(defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
		WriteSlabbedDataFourier(0.0, 0, file_info->file_space[DSET_Z_MINUS], file_info->data_set[DSET_Z_MINUS], file_info->mem_space[DSET_Z_MINUS], file_info->COMPLEX_DTYPE, &(run_data->z_minus[2]), "ZMinus", count2D, index2D);
		#endif
		#endif

		///--------------------------------------- Magnetic Amplitudes
		#if defined(__MAG_AMP) && (defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_MAG_AMP], file_info->data_set[DSET_MAG_AMP], file_info->mem_space[DSET_MAG_AMP], H5T_NATIVE_DOUBLE, &(run_data->b_n[2]), "MagAmps", count2D, index2D);
		#endif
		///--------------------------------------- Magnetic Phases
		#if defined(__MAG_PSI) && (defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_MAG_PSI], file_info->data_set[DSET_MAG_PSI], file_info->mem_space[DSET_MAG_PSI], H5T_NATIVE_DOUBLE, &(run_data->psi_n[2]), "MagPhases", count2D, index2D);
		#endif
		#endif

		///--------------------------------------- Forcing
		#if defined(__FORCING)
		if(!(strcmp(sys_vars->forcing, "STOC"))) {
			// Velocity Forcing
			WriteSlabbedDataFourier(0.0, 0, file_info->file_space[DSET_FORCING_U], file_info->data_set[DSET_FORCING_U], file_info->mem_space[DSET_FORCING_U], file_info->COMPLEX_DTYPE, &(run_data->forcing_u[2]), "VelocityForcingInTime", count2D, index2D);
			// Magnetic Forcing
			WriteSlabbedDataFourier(0.0, 0, file_info->file_space[DSET_FORCING_B], file_info->data_set[DSET_FORCING_B], file_info->mem_space[DSET_FORCING_B], file_info->COMPLEX_DTYPE, &(run_data->forcing_b[2]), "MagneticForcingInTime", count2D, index2D);
		}
		#endif

		///--------------------------------------- Energy Spectrum
		#if defined(__ENRG_SPECT)
		// Energy Spectrum
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_ENRG_SPECT], file_info->data_set[DSET_ENRG_SPECT], file_info->mem_space[DSET_ENRG_SPECT], H5T_NATIVE_DOUBLE, run_data->energy_spect, "EnergySpectrum", count2D, index2D);
		#endif
		#if defined(__KIN_ENRG_SPECT)
		// Kinetic Energy Spectrum
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_KIN_ENRG_SPECT], file_info->data_set[DSET_KIN_ENRG_SPECT], file_info->mem_space[DSET_KIN_ENRG_SPECT], H5T_NATIVE_DOUBLE, run_data->kin_enrg_spect, "KineticEnergySpectrum", count2D, index2D);
		#endif
		#if defined(__MAG_ENRG_SPECT) && (defined(__ELSASSAR_MHD) || defined(__MAGNETO))
		// Magnetic Energy Spectrum
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_MAG_ENRG_SPECT], file_info->data_set[DSET_MAG_ENRG_SPECT], file_info->mem_space[DSET_MAG_ENRG_SPECT], H5T_NATIVE_DOUBLE, run_data->mag_enrg_spect, "MagneticEnergySpectrum", count2D, index2D);
		#endif
		#if defined(__DISS_SPECT)
		// Dissipation Spectrum
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_DISS_SPECT], file_info->data_set[DSET_DISS_SPECT], file_info->mem_space[DSET_DISS_SPECT], H5T_NATIVE_DOUBLE, run_data->diss_spect, "DissipationSpectrum", count2D, index2D);
		#endif

		///--------------------------------------- Energy Flux
		#if defined(__ENRG_FLUX)
		// Write slabbed dataset for the flux
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_ENRG_FLUX], file_info->data_set[DSET_ENRG_FLUX], file_info->mem_space[DSET_ENRG_FLUX], H5T_NATIVE_DOUBLE, run_data->energy_flux, "EnergyFlux", count2D, index2D);
		// Write slabbed dataset for the diss for the velocity field
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_ENRG_DISS_VEL], file_info->data_set[DSET_ENRG_DISS_VEL], file_info->mem_space[DSET_ENRG_DISS_VEL], H5T_NATIVE_DOUBLE, run_data->energy_diss_u, "VelEnergyDiss", count2D, index2D);
		// Write slabbed dataset for the input for the velocity field
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_ENRG_INPT_VEL], file_info->data_set[DSET_ENRG_INPT_VEL], file_info->mem_space[DSET_ENRG_INPT_VEL], H5T_NATIVE_DOUBLE, run_data->energy_input_u, "VelEnergyInput", count2D, index2D);
		#if defined(__MAGNETO)
		// Write slabbed dataset for the diss for the magnetic field
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_ENRG_DISS_MAG], file_info->data_set[DSET_ENRG_DISS_MAG], file_info->mem_space[DSET_ENRG_DISS_MAG], H5T_NATIVE_DOUBLE, run_data->energy_diss_b, "MagEnergyDiss", count2D, index2D);
		// Write slabbed dataset for the input for the magnetic field
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_ENRG_INPT_MAG], file_info->data_set[DSET_ENRG_INPT_MAG], file_info->mem_space[DSET_ENRG_INPT_MAG], H5T_NATIVE_DOUBLE, run_data->energy_input_b, "MagEnergyInput", count2D, index2D);
		#endif
		#endif

		//--------------------------------------- Triads & Phase Differences
		#if defined(__CONSERVED_PHASES)
		// Write slabbed dataset for the velocity triads
		count2D[1] = phase_sync->num_triads;
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_VEL_TRIADS], file_info->data_set[DSET_VEL_TRIADS], file_info->mem_space[DSET_VEL_TRIADS], H5T_NATIVE_DOUBLE, phase_sync->triads_u, "VelTriads", count2D, index2D);
		// Write slabbed dataset for the velocity phase differrences
		count2D[1] = phase_sync->num_phase_diff;
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_VEL_PHASE_DIFF], file_info->data_set[DSET_VEL_PHASE_DIFF], file_info->mem_space[DSET_VEL_PHASE_DIFF], H5T_NATIVE_DOUBLE, phase_sync->phase_diff_u, "VelPhaseDifferences", count2D, index2D);
		#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
		// Write slabbed dataset for the magnetic triads
		count3D[0] = 1;
		count3D[1] = NUM_MAG_TRIAD_TYPES;
		count3D[2] = phase_sync->num_triads;
		index3D[0] = 0;
		index3D[1] = 0;
		index3D[2] = 0;
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_MAG_TRIADS], file_info->data_set[DSET_MAG_TRIADS], file_info->mem_space[DSET_MAG_TRIADS], H5T_NATIVE_DOUBLE, phase_sync->triads_b, "MagTriads", count3D, index3D);
		// Write slabbed dataset for the magnetic phase differences
		count2D[1] = phase_sync->num_phase_diff;
		WriteSlabbedDataReal(0.0, 0, file_info->file_space[DSET_MAG_PHASE_DIFF], file_info->data_set[DSET_MAG_PHASE_DIFF], file_info->mem_space[DSET_MAG_PHASE_DIFF], H5T_NATIVE_DOUBLE, phase_sync->phase_diff_b, "MagPhaseDifferences", count2D, index2D);
		#endif
		#endif

		//--------------------------------------- Phase Sync
		#if defined(__PHASE_SYNC)
		// Record the time averaged order parameters here
		double complex* tmp_vel_triad_order = (double complex* )malloc(sizeof(double complex) * phase_sync->num_triads);
		double complex* tmp_vel_phase_diff_order = (double complex* )malloc(sizeof(double complex) * phase_sync->num_phase_diff);
		for (int i = 0; i < phase_sync->num_triads; ++i) {
			tmp_vel_triad_order[i] = phase_sync->triad_u_order[i] / phase_sync->num_phase_sync_steps;
			if (i < phase_sync->num_phase_diff) {
				tmp_vel_phase_diff_order[i] = phase_sync->phase_diff_u_order[i] / phase_sync->num_phase_sync_steps;
			}
		}

		// Write slabbed dataset for the velocity triad order parameter
		count2D[1] = phase_sync->num_triads;
		WriteSlabbedDataFourier(0.0, 0, file_info->file_space[DSET_VEL_TRIAD_ORDER], file_info->data_set[DSET_VEL_TRIAD_ORDER], file_info->mem_space[DSET_VEL_TRIAD_ORDER], file_info->COMPLEX_DTYPE, tmp_vel_triad_order, "VelTriadOrderParameter", count2D, index2D);

		// Write slabbed dataset for the velocity phase differrence order parameter
		count2D[1] = phase_sync->num_phase_diff;
		WriteSlabbedDataFourier(0.0, 0, file_info->file_space[DSET_VEL_PHASE_DIFF_ORDER], file_info->data_set[DSET_VEL_PHASE_DIFF_ORDER], file_info->mem_space[DSET_VEL_PHASE_DIFF_ORDER], file_info->COMPLEX_DTYPE, tmp_vel_phase_diff_order, "VelPhaseDifferenceOrderParameter", count2D, index2D);
		
		// Free tmp memory
		free(tmp_vel_triad_order);
		free(tmp_vel_phase_diff_order);
		#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
		// Record the time averaged order parameters here
		double complex* tmp_mag_triad_order = (double complex* )malloc(sizeof(double complex) * NUM_MAG_TRIAD_TYPES * phase_sync->num_triads);
		double complex* tmp_mag_phase_diff_order = (double complex* )malloc(sizeof(double complex) * phase_sync->num_phase_diff);
		for (int i = 0; i < phase_sync->num_triads; ++i) {
			if (i < phase_sync->num_phase_diff) {
				tmp_mag_phase_diff_order[i] = phase_sync->phase_diff_b_order[i] / phase_sync->num_phase_sync_steps;
			}
			for (int type = 0; type < NUM_MAG_TRIAD_TYPES; ++type) {
				tmp_mag_triad_order[type * phase_sync->num_triads + i] = phase_sync->triad_b_order[type * phase_sync->num_triads + i] / phase_sync->num_phase_sync_steps;
			}
		}

		// Write slabbed dataset for the magnetic triad order parameter
		count3D[0] = 1;
		count3D[1] = NUM_MAG_TRIAD_TYPES;
		count3D[2] = phase_sync->num_triads;
		index3D[0] = 0;
		index3D[1] = 0;
		index3D[2] = 0;
		WriteSlabbedDataFourier(0.0, 0, file_info->file_space[DSET_MAG_TRIAD_ORDER], file_info->data_set[DSET_MAG_TRIAD_ORDER], file_info->mem_space[DSET_MAG_TRIAD_ORDER], file_info->COMPLEX_DTYPE, tmp_mag_triad_order, "MagTriadOrderParameter", count3D, index3D);

		// Write slabbed dataset for the magnetic phase difference order parameter
		count2D[1] = phase_sync->num_phase_diff;
		WriteSlabbedDataFourier(0.0, 0, file_info->file_space[DSET_MAG_PHASE_DIFF_ORDER], file_info->data_set[DSET_MAG_PHASE_DIFF_ORDER], file_info->mem_space[DSET_MAG_PHASE_DIFF_ORDER], file_info->COMPLEX_DTYPE, tmp_mag_phase_diff_order, "MagPhaseDifferenceOrderParameter", count2D, index2D);
		
		// Free temporary memory
		free(tmp_mag_triad_order);
		free(tmp_mag_phase_diff_order);
		#endif
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
	#if defined(__CONSERVED_PHASES) || defined(__PHASE_SYNC) || defined(__PHASE_SYNC)
	status = H5Fclose(file_info->phase_sync_file_handle);
	if (status < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to close phase sync file ["CYAN"%s"RESET"] at: Iter = ["CYAN"%d"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", file_info->phase_sync_file_name, 0, 0.0);
		exit(1);		
	}
	#endif
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
		#if defined(__MAGNETO) && !defined(__ELSASSAR_MHD)
		sprintf(sys_type, "%s", "MAGHYDRO");
		#elif defined(__ELSASSAR_MHD) && !defined(__MAGNETO)
		sprintf(sys_type, "%s", "ELSASSARMHD");
		#elif !defined(__MAGNETO) || !defined(__ELSASSAR_MHD)
		sprintf(sys_type, "%s", "HYDRO");
		#else
		sprintf(sys_type, "%s", "UKN");
		#endif
		#if defined(INT_FAC_RK4)
		sprintf(solv_type, "%s", "INTFACRK4");
		#elif defined(RK4)
		sprintf(solv_type, "%s", "RK4");
		#elif defined(AB4CN)
		sprintf(solv_type, "%s", "AB4CN");
		#else 
		sprintf(solv_type, "%s", "UKN");
		#endif
		#if defined(PHASE_ONLY_FXD_AMP)
		sprintf(model_type, "%s", "POFXDAMP");
		#elif defined(PHASE_ONLY)
		sprintf(model_type, "%s", "PO");
		#else
		sprintf(model_type, "%s", "FULL");
		#endif

		// -------------------------------------
		// Get File Label from Simulation Data
		// -------------------------------------
		#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
		// Construct file label from simulation data
		sprintf(file_data, "_SIM_DATA_[%s-%s-%s]_N[%ld]_T[%1.1lf,%g,%1.3lf]_NU[%g]_ETA[%g]_ALPHA[%1.3lf]_BETA[%1.3lf]_K[%1.3lf,%1.3lf]_EPS[%1.2lf,%1.2lf]_FORC[%s,%d,%1.3lf]_u0[%s].h5", 
							sys_type, solv_type, model_type, 
							sys_vars->N, 
							sys_vars->t0, sys_vars->dt, sys_vars->T, 
							sys_vars->NU, 
							sys_vars->ETA, 
							sys_vars->ALPHA, 
							sys_vars->BETA, 
							sys_vars->k_0, sys_vars->Lambda, 
							sys_vars->EPS, sys_vars->EPS_M,
							sys_vars->forcing, sys_vars->force_k, sys_vars->force_scale_var, 
							sys_vars->u0);
		#else
		sprintf(file_data, "_SIM_DATA_[%s-%s-%s]_N[%ld]_T[%1.1lf,%g,%1.3lf]_NU[%g]_ALPHA[%1.3lf]_K[%1.3lf,%1.3lf]_EPS[%1.2lf]_FORC[%s,%d,%1.3lf]_u0[%s].h5", 
							sys_type, solv_type, model_type, 
							sys_vars->N, 
							sys_vars->t0, sys_vars->dt, sys_vars->T, 
							sys_vars->NU, 
							sys_vars->ALPHA, 
							sys_vars->k_0, sys_vars->Lambda, 
							sys_vars->EPS,
							sys_vars->forcing, sys_vars->force_k, sys_vars->force_scale_var, 
							sys_vars->u0);
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
		#if defined(__MAGNETO) && !defined(__ELSASSAR_MHD)
		sprintf(sys_type, "%s", "MAGHYDRO");
		#elif defined(__ELSASSAR_MHD) && !defined(__MAGNETO)
		sprintf(sys_type, "%s", "ELSASSARMHD");
		#elif !defined(__MAGNETO) || !defined(__ELSASSAR_MHD)
		sprintf(sys_type, "%s", "HYDRO");
		#else
		sprintf(sys_type, "%s", "UKN");
		#endif
		#if defined(INT_FAC_RK4)
		sprintf(solv_type, "%s", "INTFACRK4");
		#elif defined(AB4CN)
		sprintf(solv_type, "%s", "AB4CN");
		#elif defined(RK4)
		sprintf(solv_type, "%s", "RK4");
		#else 
		sprintf(solv_type, "%s", "UKN");
		#endif
		#if defined(PHASE_ONLY_FXD_AMP)
		sprintf(model_type, "%s", "POFXDAMP");
		#elif defined(PHASE_ONLY)
		sprintf(model_type, "%s", "PO");
		#else
		sprintf(model_type, "%s", "FULL");
		#endif

		// ----------------------------------
		// Construct Output folder
		// ----------------------------------
		#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
		// Construct file label from simulation data
		sprintf(file_data, "SIM_DATA_[%s-%s-%s]_N[%ld]_T[%1.1lf,%g,%1.3lf]_NU[%g]_ETA[%g]_ALPHA[%1.3lf]_BETA[%1.3lf]_K[%1.3lf,%1.3lf]_EPS[%1.2lf,%1.2lf]_FORC[%s,%d,%1.3lf]_u0[%s]_TAG[%s]/", 
							sys_type, solv_type, model_type, 
							sys_vars->N, 
							sys_vars->t0, sys_vars->dt, sys_vars->T, 
							sys_vars->NU, 
							sys_vars->ETA, 
							sys_vars->ALPHA,
							sys_vars->BETA, 
							sys_vars->k_0, sys_vars->Lambda, 
							sys_vars->EPS, sys_vars->EPS_M,
							sys_vars->forcing, sys_vars->force_k, sys_vars->force_scale_var, 
							sys_vars->u0, 
							file_info->output_tag);
		#else
		sprintf(file_data, "SIM_DATA_[%s-%s-%s]_N[%ld]_T[%1.1lf,%g,%1.3lf]_NU[%g]_ALPHA[%1.3lf]_K[%1.3lf,%1.3lf]_EPS[%1.2lf]_FORC[%s,%d,%1.3lf]_u0[%s]_TAG[%s]/", 
							sys_type, solv_type, model_type, 
							sys_vars->N, 
							sys_vars->t0, sys_vars->dt, sys_vars->T, 
							sys_vars->NU, sys_vars->ALPHA, 
							sys_vars->k_0, sys_vars->Lambda, 
							sys_vars->EPS,
							sys_vars->forcing, sys_vars->force_k, sys_vars->force_scale_var, 
							sys_vars->u0, 
							file_info->output_tag);
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
		#if defined(__STATS)
		strcpy(file_info->stats_file_name, file_info->output_dir); 
		strcat(file_info->stats_file_name, "Stats_HDF_Data.h5");
		printf("\nStats File: "CYAN"%s"RESET"\n\n", file_info->stats_file_name);
		#endif
		#if defined(__SYS_MEASURES)
		strcpy(file_info->system_msr_file_name, file_info->output_dir); 
		strcat(file_info->system_msr_file_name, "System_Measure_HDF_Data.h5");
		printf("\nSytem Measures File: "CYAN"%s"RESET"\n\n", file_info->system_msr_file_name);
		#endif
		#if defined(__CONSERVED_PHASES) || defined(__PHASE_SYNC) || defined(__PHASE_SYNC)
		strcpy(file_info->phase_sync_file_name, file_info->output_dir); 
		strcat(file_info->phase_sync_file_name, "Phase_Sync_HDF_Data.h5");
		printf("\nPhase Sync File: "CYAN"%s"RESET"\n\n", file_info->phase_sync_file_name);
		#endif
	}

}
/**
 * Wrapper function that writes the data to file by openining it, creating a group for the current iteration and writing the data under this group. The file is then closed again 
 * @param t     	The current time of the simulation
 * @param iters 	The current iteration
 * @param save_indx The current iteration
 */
void WriteDataToFile(double t, const long int iters, const long int save_indx) {

	// Initialize variables
	herr_t status;

	// --------------------------------------
	// Check if files exist and Open/Create
	// --------------------------------------
	// Check if main file exists - open it if it does if not create it
	if (access(file_info->output_file_name, F_OK) != 0) {
		file_info->output_file_handle = H5Fcreate(file_info->output_file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		if (file_info->output_file_handle < 0) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to create output file ["CYAN"%s"RESET"] at: Iter = ["CYAN"%ld"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", file_info->output_file_name, iters, t);
			exit(1);
		}
	}
	else {
		// Open file with serial I/O access properties
		file_info->output_file_handle = H5Fopen(file_info->output_file_name, H5F_ACC_RDWR, H5P_DEFAULT);
		if (file_info->output_file_handle < 0) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to open output file ["CYAN"%s"RESET"] at: Iter = ["CYAN"%ld"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", file_info->output_file_name, iters, t);
			exit(1);
		}
	}
	#if defined(__CONSERVED_PHASES) || defined(__PHASE_SYNC) || defined(__PHASE_SYNC)
	if (access(file_info->phase_sync_file_name, F_OK) != 0) {
		file_info->phase_sync_file_handle = H5Fcreate(file_info->phase_sync_file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		if (file_info->phase_sync_file_handle < 0) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to create phase_sync file ["CYAN"%s"RESET"] at: Iter = ["CYAN"%ld"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", file_info->phase_sync_file_name, iters, t);
			exit(1);
		}
	}
	else {
		// Open file with serial I/O access properties
		file_info->phase_sync_file_handle = H5Fopen(file_info->phase_sync_file_name, H5F_ACC_RDWR, H5P_DEFAULT);
		if (file_info->phase_sync_file_handle < 0) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to open phase_sync file ["CYAN"%s"RESET"] at: Iter = ["CYAN"%ld"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", file_info->phase_sync_file_name, iters, t);
			exit(1);
		}
	}
	#endif

	// -------------------------------
	// Write Data to File
	// -------------------------------
	// Initialize slab dimension arrays
	const int Dim2D = 2;
	hsize_t index2D[Dim2D]; 		// stores the index in the hyperslabbed dataset to start writing to
	hsize_t count2D[Dim2D];       	// stores the size of hyperslab to write to the dataset
	#if (defined(__MAGNETO) || defined(__ELSASSAR_MHD)) && (defined(__CONSERVED_PHASES) || defined(__PHASE_SYNC))
	const int Dim3D = 3;
	hsize_t index3D[Dim3D]; 		// stores the index in the hyperslabbed dataset to start writing to
	hsize_t count3D[Dim3D];       	// stores the size of hyperslab to write to the dataset
	#endif

	// Initialize the hyperslab arrays
	count2D[0] = 1;					// 1D slab so first dim is 1
	count2D[1] = sys_vars->N;		// 1D slab of size of data array
	index2D[0] = save_indx;			// set the starting row index to index in the global dataset to write slab to -> 0th row for 0 iter
	index2D[1] = 0;					// set column index to 0 to start writing from the first column

	///--------------------------------------- Velocity Modes
	#if defined(__VEL) && !(defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY)) 
	WriteSlabbedDataFourier(t, iters, file_info->file_space[DSET_VEL], file_info->data_set[DSET_VEL], file_info->mem_space[DSET_VEL], file_info->COMPLEX_DTYPE, &(run_data->u[2]), "VelModes", count2D, index2D);
	#endif

	///--------------------------------------- Velocity Amplitudes
	#if defined(__VEL_AMP) && (defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_VEL_AMP], file_info->data_set[DSET_VEL_AMP], file_info->mem_space[DSET_VEL_AMP], H5T_NATIVE_DOUBLE, &(run_data->a_n[2]), "VelAmps", count2D, index2D);
	#endif
	///--------------------------------------- Velocity Phases
	#if defined(__VEL_PHI) && (defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_VEL_PHI], file_info->data_set[DSET_VEL_PHI], file_info->mem_space[DSET_VEL_PHI], H5T_NATIVE_DOUBLE, &(run_data->phi_n[2]), "VelPhases", count2D, index2D);
	#endif

	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	///--------------------------------------- Magnetic Modes
	#if defined(__MAG) && !(defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
	WriteSlabbedDataFourier(t, iters, file_info->file_space[DSET_MAG], file_info->data_set[DSET_MAG], file_info->mem_space[DSET_MAG], file_info->COMPLEX_DTYPE, &(run_data->b[2]), "MagModes", count2D, index2D);
	#endif

	#if defined(__ELSASSAR_MHD)
	///--------------------------------------- Z_Plus Elsassar Variable
	#if defined(__Z_PLUS) && !(defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
	WriteSlabbedDataFourier(t, iters, file_info->file_space[DSET_Z_PLUS], file_info->data_set[DSET_Z_PLUS], file_info->mem_space[DSET_Z_PLUS], file_info->COMPLEX_DTYPE, &(run_data->z_plus[2]), "ZPlus", count2D, index2D);
	#endif
	///--------------------------------------- Z_Minus Elsassar Variables
	#if defined(__Z_MINUS) && !(defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
	WriteSlabbedDataFourier(t, iters, file_info->file_space[DSET_Z_MINUS], file_info->data_set[DSET_Z_MINUS], file_info->mem_space[DSET_Z_MINUS], file_info->COMPLEX_DTYPE, &(run_data->z_minus[2]), "ZMinus", count2D, index2D);
	#endif
	#endif

	///--------------------------------------- Magnetic Amplitudes
	#if defined(__MAG_AMP) && (defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_MAG_AMP], file_info->data_set[DSET_MAG_AMP], file_info->mem_space[DSET_MAG_AMP], H5T_NATIVE_DOUBLE, &(run_data->b_n[2]), "MagAmps", count2D, index2D);
	#endif
	///--------------------------------------- Magnetic Phases
	#if defined(__MAG_PSI) && (defined(PHASE_ONLY_FXD_AMP) || defined(PHASE_ONLY))
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_MAG_PSI], file_info->data_set[DSET_MAG_PSI], file_info->mem_space[DSET_MAG_PSI], H5T_NATIVE_DOUBLE, &(run_data->psi_n[2]), "MagPhases", count2D, index2D);
	#endif
	#endif

	///--------------------------------------- Forcing
	#if defined(__FORCING)
	if(!(strcmp(sys_vars->forcing, "STOC"))) {
		// Velocity Forcing
		WriteSlabbedDataFourier(t, iters, file_info->file_space[DSET_FORCING_U], file_info->data_set[DSET_FORCING_U], file_info->mem_space[DSET_FORCING_U], file_info->COMPLEX_DTYPE, &(run_data->forcing_u[2]), "VelocityForcingInTime", count2D, index2D);
		// Magnetic Forcing
		WriteSlabbedDataFourier(t, iters, file_info->file_space[DSET_FORCING_B], file_info->data_set[DSET_FORCING_B], file_info->mem_space[DSET_FORCING_B], file_info->COMPLEX_DTYPE, &(run_data->forcing_b[2]), "MagneticForcingInTime", count2D, index2D);
	}
	#endif

	///--------------------------------------- Energy Spectrum
	#if defined(__ENRG_SPECT)
	// Energy Spectrum
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_ENRG_SPECT], file_info->data_set[DSET_ENRG_SPECT], file_info->mem_space[DSET_ENRG_SPECT], H5T_NATIVE_DOUBLE, run_data->energy_spect, "EnergySpectrum", count2D, index2D);
	#endif
	#if defined(__KIN_ENRG_SPECT)
	// Kinetic Energy Spectrum
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_KIN_ENRG_SPECT], file_info->data_set[DSET_KIN_ENRG_SPECT], file_info->mem_space[DSET_KIN_ENRG_SPECT], H5T_NATIVE_DOUBLE, run_data->kin_enrg_spect, "KineticEnergySpectrum", count2D, index2D);
	#endif
	#if defined(__MAG_ENRG_SPECT) && (defined(__ELSASSAR_MHD) || defined(__MAGNETO))
	// Magnetic Energy Spectrum
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_MAG_ENRG_SPECT], file_info->data_set[DSET_MAG_ENRG_SPECT], file_info->mem_space[DSET_MAG_ENRG_SPECT], H5T_NATIVE_DOUBLE, run_data->mag_enrg_spect, "MagneticEnergySpectrum", count2D, index2D);
	#endif
	#if defined(__ENRG_SPECT)
	// Dissipation Spectrum
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_DISS_SPECT], file_info->data_set[DSET_DISS_SPECT], file_info->mem_space[DSET_DISS_SPECT], H5T_NATIVE_DOUBLE, run_data->diss_spect, "DissipationSpectrum", count2D, index2D);
	#endif

	///--------------------------------------- Energy Flux
	#if defined(__ENRG_FLUX)
	// Write slabbed dataset for the flux
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_ENRG_FLUX], file_info->data_set[DSET_ENRG_FLUX], file_info->mem_space[DSET_ENRG_FLUX], H5T_NATIVE_DOUBLE, run_data->energy_flux, "EnergyFlux", count2D, index2D);
	// Write slabbed dataset for the diss for the velocity field
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_ENRG_DISS_VEL], file_info->data_set[DSET_ENRG_DISS_VEL], file_info->mem_space[DSET_ENRG_DISS_VEL], H5T_NATIVE_DOUBLE, run_data->energy_diss_u, "VelEnergyDiss", count2D, index2D);
	// Write slabbed dataset for the input for the velocity field
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_ENRG_INPT_VEL], file_info->data_set[DSET_ENRG_INPT_VEL], file_info->mem_space[DSET_ENRG_INPT_VEL], H5T_NATIVE_DOUBLE, run_data->energy_input_u, "VelEnergyInput", count2D, index2D);
	#if defined(__MAGNETO)
	// Write slabbed dataset for the diss for the magnetic field
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_ENRG_DISS_MAG], file_info->data_set[DSET_ENRG_DISS_MAG], file_info->mem_space[DSET_ENRG_DISS_MAG], H5T_NATIVE_DOUBLE, run_data->energy_diss_b, "MagEnergyDiss", count2D, index2D);
	// Write slabbed dataset for the input for the magnetic field
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_ENRG_INPT_MAG], file_info->data_set[DSET_ENRG_INPT_MAG], file_info->mem_space[DSET_ENRG_INPT_MAG], H5T_NATIVE_DOUBLE, run_data->energy_input_b, "MagEnergyInput", count2D, index2D);
	#endif
	#endif


	//--------------------------------------- Triads & Phase Differences
	#if defined(__CONSERVED_PHASES)
	// Write slabbed dataset for the velocity triads
	count2D[0] = 1;
	count2D[1] = phase_sync->num_triads;
	index2D[0] = save_indx;
	index2D[1] = 0;
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_VEL_TRIADS], file_info->data_set[DSET_VEL_TRIADS], file_info->mem_space[DSET_VEL_TRIADS], H5T_NATIVE_DOUBLE, phase_sync->triads_u, "VelTriads", count2D, index2D);
	// Write slabbed dataset for the velocity phase differrences
	count2D[1] = phase_sync->num_phase_diff;
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_VEL_PHASE_DIFF], file_info->data_set[DSET_VEL_PHASE_DIFF], file_info->mem_space[DSET_VEL_PHASE_DIFF], H5T_NATIVE_DOUBLE, phase_sync->phase_diff_u, "VelPhaseDifferences", count2D, index2D);
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	// Write slabbed dataset for the magnetic triads
	count3D[0] = 1;
	count3D[1] = NUM_MAG_TRIAD_TYPES;
	count3D[2] = phase_sync->num_triads;
	index3D[0] = save_indx;
	index3D[1] = 0;
	index3D[2] = 0;
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_MAG_TRIADS], file_info->data_set[DSET_MAG_TRIADS], file_info->mem_space[DSET_MAG_TRIADS], H5T_NATIVE_DOUBLE, phase_sync->triads_b, "MagTriads", count3D, index3D);
	// Write slabbed dataset for the magnetic phase differences
	count3D[1] = phase_sync->num_phase_diff;
	WriteSlabbedDataReal(t, iters, file_info->file_space[DSET_MAG_PHASE_DIFF], file_info->data_set[DSET_MAG_PHASE_DIFF], file_info->mem_space[DSET_MAG_PHASE_DIFF], H5T_NATIVE_DOUBLE, phase_sync->phase_diff_b, "MagPhaseDifferences", count3D, index3D);
	#endif
	#endif
	
	//--------------------------------------- Phase Sync
	#if defined(__PHASE_SYNC)
	// Record the time averaged order parameters here
	double complex* tmp_vel_triad_order_n = (double complex* )malloc(sizeof(double complex) * phase_sync->num_triads);
	double complex* tmp_vel_phase_diff_order_n = (double complex* )malloc(sizeof(double complex) * phase_sync->num_phase_diff);
	for (int i = 0; i < phase_sync->num_triads; ++i) {
		tmp_vel_triad_order_n[i] = phase_sync->triad_u_order[i] / phase_sync->num_phase_sync_steps;
		if (i < phase_sync->num_phase_diff) {
			tmp_vel_phase_diff_order_n[i] = phase_sync->phase_diff_u_order[i] / phase_sync->num_phase_sync_steps;
		}
	}

	// Write slabbed dataset for the velocity triad order parameter
	count2D[1]       = phase_sync->num_triads;
	WriteSlabbedDataFourier(t, iters, file_info->file_space[DSET_VEL_TRIAD_ORDER], file_info->data_set[DSET_VEL_TRIAD_ORDER], file_info->mem_space[DSET_VEL_TRIAD_ORDER], file_info->COMPLEX_DTYPE, tmp_vel_triad_order_n, "VelTriadOrderParameter", count2D, index2D);

	// Write slabbed dataset for the velocity phase differrence order parameter
	count2D[1]       = phase_sync->num_phase_diff;
	WriteSlabbedDataFourier(t, iters, file_info->file_space[DSET_VEL_PHASE_DIFF_ORDER], file_info->data_set[DSET_VEL_PHASE_DIFF_ORDER], file_info->mem_space[DSET_VEL_PHASE_DIFF_ORDER], file_info->COMPLEX_DTYPE, tmp_vel_phase_diff_order_n, "VelPhaseDifferenceOrderParameter", count2D, index2D);

	// Free tmp memory
	free(tmp_vel_triad_order_n);
	free(tmp_vel_phase_diff_order_n);
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	// Record the time averaged order parameters here
	double complex* tmp_mag_triad_order_n = (double complex* )malloc(sizeof(double complex) * NUM_MAG_TRIAD_TYPES * phase_sync->num_triads);
	double complex* tmp_mag_phase_diff_order_n = (double complex* )malloc(sizeof(double complex) * phase_sync->num_phase_diff);
	for (int i = 0; i < phase_sync->num_triads; ++i) {
		if (i < phase_sync->num_phase_diff) {
			tmp_mag_phase_diff_order_n[i] = phase_sync->phase_diff_b_order[i] / phase_sync->num_phase_sync_steps;
		}
		for (int type = 0; type < NUM_MAG_TRIAD_TYPES; ++type) {
			tmp_mag_triad_order_n[type * phase_sync->num_triads + i] = phase_sync->triad_b_order[type * phase_sync->num_triads + i] / phase_sync->num_phase_sync_steps;
		}
	}

	// Write slabbed dataset for the magnetic triad order parameter
	count3D[0] = 1;
	count3D[1] = NUM_MAG_TRIAD_TYPES;
	count3D[2] = phase_sync->num_triads;
	index3D[0] = save_indx;
	index3D[1] = 0;
	index3D[2] = 0;
	WriteSlabbedDataFourier(t, iters, file_info->file_space[DSET_MAG_TRIAD_ORDER], file_info->data_set[DSET_MAG_TRIAD_ORDER], file_info->mem_space[DSET_MAG_TRIAD_ORDER], file_info->COMPLEX_DTYPE, tmp_mag_triad_order_n, "MagTriadOrderParameter", count3D, index3D);
	
	// Write slabbed dataset for the magnetic phase difference order parameter
	count2D[1]       = phase_sync->num_phase_diff;
	WriteSlabbedDataFourier(t, iters, file_info->file_space[DSET_MAG_PHASE_DIFF_ORDER], file_info->data_set[DSET_MAG_PHASE_DIFF_ORDER], file_info->mem_space[DSET_MAG_PHASE_DIFF_ORDER], file_info->COMPLEX_DTYPE, tmp_mag_phase_diff_order_n, "MagPhaseDifferenceOrderParameter", count2D, index2D);

	// Free tmp memory
	free(tmp_mag_triad_order_n);
	free(tmp_mag_phase_diff_order_n);
	#endif
	#endif


	// -------------------------------
	// Close identifiers and File
	// -------------------------------
	status = H5Fclose(file_info->output_file_handle);
	if (status < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to close output file ["CYAN"%s"RESET"] at: Iter = ["CYAN"%ld"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", file_info->output_file_name, iters, t);
		exit(1);
	}
	#if defined(__CONSERVED_PHASES) || defined(__PHASE_SYNC) || defined(__PHASE_SYNC)
	status = H5Fclose(file_info->phase_sync_file_handle);
	if (status < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to close phase sync file ["CYAN"%s"RESET"] at: Iter = ["CYAN"%ld"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", file_info->phase_sync_file_name, iters, t);
		exit(1);		
	}
	#endif
}
/**
 * Function to open and read data from an input file - In phase only model mode the function will attempt to read the amplitudes and phases, if phases don't exist, 
 * uniformly randomly generated phases will be used. In full mode the function will try the modes first, if not it will then try the amplitudes and phases. 
 * If no amplitudes exist in both modes the programmes exits!
 * @param N The number of shells/size of the input arrays
 */
void ReadInputFile(const long int N) {

	// Initialize variables
	int n;
	herr_t status;

	// -------------------------------
	// Open Input File
	// -------------------------------
	// Open file with serial I/O access properties
	file_info->input_file_handle = H5Fopen(file_info->input_file_name, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_info->input_file_handle < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to open input file ["CYAN"%s"RESET"]\n-->> Exiting...\n", file_info->input_file_name);
		exit(1);
	}

	// -------------------------------
	// Read In Initial Condition
	// -------------------------------
	#if defined(PHASE_ONLY)
	///----------------------------- Read in initial velocity amps and phases
	// Create tmp array to read in data
	double* tmp_u_amp = (double* )malloc(sizeof(double) * N);
	if ( (H5LTread_dataset(file_info->input_file_handle, "VelAmps", H5T_NATIVE_DOUBLE, tmp_u_amp)) < 0) {
		printf("\n["RED"ERROR"RESET"] --- Input Dataset ["CYAN"%s"RESET"] does not exist\n---> Exiting!!!\n", "VelAmps");
		exit(1);
	}
	double* tmp_u_phase = (double* )malloc(sizeof(double) * N);
	if ( (H5LTread_dataset(file_info->input_file_handle, "VelPhases", H5T_NATIVE_DOUBLE, tmp_u_phase)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to read input dataset ["CYAN"%s"RESET"] ---> Using uniformly random generated phases instead\n", "VelPhases");
		for (int i = 0; i < N; ++i) {
			tmp_u_phase[i] = (double)rand() / (double)RAND_MAX * 2.0 * M_PI;
		}
	}
	
	// Write tmp arrays to velocity modes / mode amplitudes and phases
	for (int i = 0; i < N + 4; ++i) {
		if (i >= 2 && i < N + 2) {
			run_data->a_n[i]   = tmp_u_amp[i - 2];
			run_data->phi_n[i] = tmp_u_phase[i - 2];
			run_data->u[i]     = tmp_u_amp[i - 2] * cexp(I * tmp_u_phase[i - 2]);
		}
		else {
			run_data->a_n[i]   = 0.0;
			run_data->phi_n[i] = 0.0;
			run_data->u[i]     = 0.0  + 0.0 * I;	
		}
	}

	// Free temp memory
	free(tmp_u_amp);
	free(tmp_u_phase);

	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	///----------------------------- Read in initial magnetic amps and phases
	// Create tmp array to read in data
	double* tmp_b_amp = (double* )malloc(sizeof(double) * N);
	if ( (H5LTread_dataset(file_info->input_file_handle, "MagAmps", H5T_NATIVE_DOUBLE, tmp_b_amp)) < 0) {
		printf("\n["RED"ERROR"RESET"] --- Input Dataset ["CYAN"%s"RESET"] does not exist\n---> Exiting!!!\n", "MagAmps");
		exit(1);
	}
	double* tmp_b_phase = (double* )malloc(sizeof(double) * N);
	if ( (H5LTread_dataset(file_info->input_file_handle, "MagPhases", H5T_NATIVE_DOUBLE, tmp_b_phase)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to read input dataset ["CYAN"%s"RESET"] ---> Using uniformly random generated phases instead\n", "MagPhases");
		for (int i = 0; i < N; ++i) {
			tmp_b_phase[i] = (double)rand() / (double)RAND_MAX * 2.0 * M_PI;
		}
	}
	
	// Write tmp array to magnetic modes / mode amplitudes and phases
	for (int i = 0; i < N + 4; ++i) {
		if (i >= 2 && i < N + 2) {
			run_data->b_n[i]   = tmp_b_amp[i - 2];
			run_data->psi_n[i] = tmp_b_phase[i - 2];
			run_data->b[i]     = tmp_b_amp[i - 2] * cexp(I * tmp_b_phase[i - 2]);
		}
		else {
			run_data->b_n[i]   = 0.0;
			run_data->psi_n[i] = 0.0;
			run_data->b[i]     = 0.0  + 0.0 * I;	
		}
	}

	// Free temp memory
	free(tmp_b_amp);
	free(tmp_b_phase);
	#endif
	#else
	// Create compound datatype for the complex datasets
	file_info->COMPLEX_DTYPE = CreateComplexDatatype();

	///----------------------------- Read in initial velocity modes
	// Create tmp array to read in data
	double complex* tmp_u = (double complex* )malloc(sizeof(double complex) * N);
	// if ( (H5LTread_dataset(file_info->input_file_handle, "VelModes", file_info->COMPLEX_DTYPE, tmp_u)) < 0) {
	// 	// Print warning to screen
	// 	printf("\n["MAGENTA"WARNING"RESET"] --- Failed to read input dataset ["CYAN"%s"RESET"] ---> Trying Amplitude/Phase datasets\n", "VelModes");

	// Attempt to read in the initial amplitudes
	double* tmp_u_amp = (double* )malloc(sizeof(double) * N);
	if ( (H5LTread_dataset(file_info->input_file_handle, "VelAmps", H5T_NATIVE_DOUBLE, tmp_u_amp)) < 0) {
		printf("\n["RED"ERROR"RESET"] --- Input dataset ["CYAN"%s"RESET"]\n---> Exiting!!!", "VelAmps");
		exit(1);
	}
	// Attemp to read in the intial phases if not use randomly generated ones
	double* tmp_u_phase = (double* )malloc(sizeof(double) * N);
	// if ( (H5LTread_dataset(file_info->input_file_handle, "VelPhases", H5T_NATIVE_DOUBLE, tmp_u_phase)) < 0) {
	// 	// Print warning to screen
	// 	printf("\n["MAGENTA"WARNING"RESET"] --- Failed to read input dataset ["CYAN"%s"RESET"] ---> Using uniformly random generated phases instead!\n", "VelPhases");

	// Uniformly randomly generated phases
	for (int i = 0; i < N; ++i) {
		tmp_u_phase[i] = (double)rand()/(double)RAND_MAX * 2.0 * M_PI;
	}			
	// }
		
		// Create the modes from the phase/amplitudes
		for (int i = 0; i < N; ++i) {
			tmp_u[i] = tmp_u_amp[i] * cexp(I * tmp_u_phase[i]);
		}

		// Free tmp memory
		free(tmp_u_amp);
		free(tmp_u_phase);
	// }
	
	// Write tmp_u array to velocity modes / mode amplitudes and phases
	for (int i = 0; i < N + 4; ++i) {
		if (i >= 2 && i < N + 2) {
			run_data->u[i] = tmp_u[i - 2];
			#if defined(PHASE_ONLY_FXD_AMP)
			run_data->a_n[i]   = cabs(tmp_u[i - 2]);
			run_data->phi_n[i] = carg(tmp_u[i - 2]);
			#endif
		}
		else {
			run_data->u[i] = 0.0 + 0.0 * I;
			#if defined(PHASE_ONLY_FXD_AMP)
			run_data->a_n[i]   = 0.0 + 0.0 * I;
			run_data->phi_n[i] = 0.0 + 0.0 * I;
			#endif
		}
	}

	// Free temp memory
	free(tmp_u);

	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	///----------------------------- Read in initial magnetic modes
	// Create tmp array to read in data
	double complex* tmp_b = (double complex* )malloc(sizeof(double complex) * N);
	if ( (H5LTread_dataset(file_info->input_file_handle, "MagModes", file_info->COMPLEX_DTYPE, tmp_b)) < 0) {
		// Print warning to screen
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to read input dataset ["CYAN"%s"RESET"] ---> Trying Amplitude/Phase datasets\n", "MagModes");
		
		// Attempt to read in the initial amplitudes
		double* tmp_b_amp = (double* )malloc(sizeof(double) * N);
		if ( (H5LTread_dataset(file_info->input_file_handle, "MagAmps", H5T_NATIVE_DOUBLE, tmp_b_amp)) < 0) {
			printf("\n["RED"ERROR"RESET"] --- Input dataset ["CYAN"%s"RESET"]\n---> Exiting!!!", "MagAmps");
			exit(1);
		}
		// Attemp to read in the intial phases if not use randomly generated ones
		double* tmp_b_phase = (double* )malloc(sizeof(double) * N);
		if ( (H5LTread_dataset(file_info->input_file_handle, "MagPhases", H5T_NATIVE_DOUBLE, tmp_b_phase)) < 0) {
			// Print warning to screen
			printf("\n["MAGENTA"WARNING"RESET"] --- Failed to read input dataset ["CYAN"%s"RESET"] ---> Using uniformly random generated phases instead!\n", "MagPhases");

			// Uniformly randomly generated phases
			for (int i = 0; i < N; ++i) {
				tmp_b_phase[i] = (double)rand()/(double)RAND_MAX * 2.0 * M_PI;
			}
		}		
		
		// Create the modes from the phase/amplitudes
		for (int i = 0; i < N; ++i) {
			tmp_b[i] = tmp_b_amp[i] * cexp(I * tmp_b_phase[i]);
		}

		// Free tmp memory
		free(tmp_b_amp);
		free(tmp_b_phase);
	}
	
	// Write tmp_b array to magnetic modes / mode amplitudes and phases
	for (int i = 0; i < N + 4; ++i) {
		if (i >= 2 && i < N + 2) {
			run_data->b[i] = tmp_b[i - 2];
			#if defined(PHASE_ONLY_FXD_AMP)
			run_data->b_n[i]   = cabs(tmp_b[i - 2]);
			run_data->psi_n[i] = carg(tmp_b[i - 2]);
			#endif
		}
		else {
			run_data->b[i] = 0.0 + 0.0 * I;
			#if defined(PHASE_ONLY_FXD_AMP)
			run_data->b_n[i]   = 0.0 + 0.0 * I;
			run_data->psi_n[i] = 0.0 + 0.0 * I;
			#endif
		}
	}

	// Free temp memory
	free(tmp_b);
	#endif	
	#endif

	// -------------------------------
	// Close identifiers and File
	// -------------------------------
	status = H5Fclose(file_info->input_file_handle);
	if (status < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to close input file ["CYAN"%s"RESET"]\n-->> Exiting...\n", file_info->input_file_name);
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
void CreateSlabbedDSet(double t, int iters, hid_t file_handle, char* dset_name, hid_t* file_space, hid_t* data_set, hid_t* mem_space, hid_t dtype, hsize_t* dset_dims, hsize_t* dset_max_dims, hsize_t* dset_chunk_dims, const int num_dims) {

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
	*data_set = H5Dcreate(file_handle, dset_name, dtype, *file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
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
 * @param n          Array containing The size of the chunk to write
 * @param index      Array containing The index in the dataset to write chunk to
 */
void WriteSlabbedDataReal(double t, int iters, hid_t file_space, hid_t data_set, hid_t mem_space, hid_t dtype, double* data, char* dset_name, hsize_t* count, hsize_t* index) {

	// ------------------------------------
	// Select Slab
	// ------------------------------------
	// Select appropriate hyperslab 
	if ((H5Sselect_hyperslab(file_space, H5S_SELECT_SET, index, NULL, count, NULL)) < 0) {
		// Print error message to error stream 
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to select hyperslab in file for datset ["CYAN"%s"RESET"] at: Iter = ["CYAN"%d"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", dset_name, iters, t);

		// Save current state of system to file before exiting
		FinalWriteAndCloseOutputFile(sys_vars->N, iters, index[0]);
		
		// Exit programme
		exit(1);
	}

	// ------------------------------------
	// Write Chunk to File
	// ------------------------------------
	// Then write the current chunk to this hyperslab
	if ((H5Dwrite(data_set, dtype, mem_space, file_space, H5P_DEFAULT, data)) < 0) {
		// Print error message to error stream 
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to write data to datset ["CYAN"%s"RESET"] at: Iter = ["CYAN"%d"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", dset_name, iters, t);

		// Save current state of system to file before exiting
		FinalWriteAndCloseOutputFile(sys_vars->N, iters, index[0]);

		// Exit programme
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
 * @param n          Array containing The size of the chunk to write
 * @param index      Array containing The index in the dataset to write chunk to
 */
void WriteSlabbedDataFourier(double t, int iters, hid_t file_space, hid_t data_set, hid_t mem_space, hid_t dtype, double complex* data, char* dset_name, hsize_t* count, hsize_t* index) {

	// ------------------------------------
	// Select Slab
	// ------------------------------------
	// Select appropriate hyperslab 
	if ((H5Sselect_hyperslab(file_space, H5S_SELECT_SET, index, NULL, count, NULL)) < 0) {
		// Print error message to error stream
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to select hyperslab in file for datset ["CYAN"%s"RESET"] at: Iter = ["CYAN"%d"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", dset_name, iters, t);
		
		// Save current state of system to file before exiting
		FinalWriteAndCloseOutputFile(sys_vars->N, iters, index[0]);
		
		// Exit programme
		exit(1);
	}

	// ------------------------------------
	// Write Chunk to File
	// ------------------------------------
	// Then write the current chunk to this hyperslab
	if ((H5Dwrite(data_set, dtype, mem_space, file_space, H5P_DEFAULT, data)) < 0) {
		// Print error message to error stream
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to write data to datset ["CYAN"%s"RESET"] at: Iter = ["CYAN"%d"RESET"] t = ["CYAN"%lf"RESET"]\n-->> Exiting...\n", dset_name, iters, t);
		
		// Save current state of system to file before exiting
		FinalWriteAndCloseOutputFile(sys_vars->N, iters, index[0]);
		
		// Exit programme
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

	herr_t status;

	// Record total iterations
	sys_vars->tot_iters      = (long int)iters - 1;
	sys_vars->tot_save_steps = (long int)save_data_indx - 1;

	////////////////////////////////
	/// Open and Write Datasets
	////////////////////////////////
	#if defined(__SYS_MEASURES)
	file_info->sys_msr_file_handle = H5Fcreate(file_info->system_msr_file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (file_info->sys_msr_file_handle < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"]  --- Could not create main HDF5 system measure file at: "CYAN"%s"RESET" \n-->>Exiting....\n", file_info->system_msr_file_name);
		exit(1);
	}
	#endif
	#if defined(__STATS)
	file_info->stats_file_handle = H5Fcreate(file_info->stats_file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (file_info->stats_file_handle < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"]  --- Could not create main HDF5 stats file at: "CYAN"%s"RESET" \n-->>Exiting....\n", file_info->stats_file_name);
		exit(1);
	}
	#endif
	#if defined(__CONSERVED_PHASES) || defined(__PHASE_SYNC) || defined(__PHASE_SYNC)
	if (access(file_info->phase_sync_file_name, F_OK) != 0) {
		file_info->phase_sync_file_handle = H5Fcreate(file_info->phase_sync_file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		if (file_info->phase_sync_file_handle < 0) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to create phase_sync file ["CYAN"%s"RESET"] at final write\n-->> Exiting...\n", file_info->phase_sync_file_name);
			exit(1);
		}
	}
	else {
		// Open file with serial I/O access properties
		file_info->phase_sync_file_handle = H5Fopen(file_info->phase_sync_file_name, H5F_ACC_RDWR, H5P_DEFAULT);
		if (file_info->phase_sync_file_handle < 0) {
			fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to open phase_sync file ["CYAN"%s"RESET"] at final write\n-->> Exiting...\n", file_info->phase_sync_file_name);
			exit(1);
		}
	}
	#endif

	// -------------------------------
	// Write System Measurables
	// -------------------------------
	WriteSystemMeasuresToFile();

	// -------------------------------
	// Write Phase Sync Stats
	// -------------------------------
	#if defined(__PHASE_SYNC_STATS)
	WritePhaseSyncStatsToFile();
	#endif

	printf("\n\nAfter Sync Before STATSS\n\n\n");
	// -------------------------------
	// Write Stats
	// -------------------------------
	#if defined(__STATS)
	WriteStatsToFile();
	#endif

	printf("\n\nAfter STATS\n\n\n");
	
	// -----------------------------------
	// Close Files for the final time
	// -----------------------------------
	status = H5Fclose(file_info->sys_msr_file_handle);
	if (status < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to close main output file: "CYAN"%s"RESET" \n-->> Exiting....\n", file_info->output_file_name);
		exit(1);
	}
	#if defined(__STATS)
	status = H5Fclose(file_info->stats_file_handle);
	if (status < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to close main stats file: "CYAN"%s"RESET" \n-->> Exiting....\n", file_info->stats_file_name);
		exit(1);
	}
	#endif
	#if defined(__CONSERVED_PHASES) || defined(__PHASE_SYNC) || defined(__PHASE_SYNC)
	status = H5Fclose(file_info->phase_sync_file_handle);
	if (status < 0) {
		fprintf(stderr, "\n["RED"ERROR"RESET"] --- Unable to close main phase_sync file: "CYAN"%s"RESET" \n-->> Exiting....\n", file_info->phase_sync_file_name);
		exit(1);
	}
	#endif

	if (sys_vars->INPUT_FILE_FLAG == INPUT_FILE) {
		#if !defined(PHASE_ONLY)
		// Close the complex datatype identifier
		H5Tclose(file_info->COMPLEX_DTYPE);
		#endif
	}
	else {
		#if defined(__VEL) || defined(__MAG) || defined(__Z_PLUS) || defined(__Z_MINUS) || defined(__FORCING) || defined(__STATS) || defined(__PHASE_SYNC)
		// Close the complex datatype identifier
		H5Tclose(file_info->COMPLEX_DTYPE);
		#endif
	}
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