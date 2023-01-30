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

	#if defined(__VEL) || defined(__MAG) || defined(__Z_PLUS) || defined(__Z_MINUS) || defined(__FORCING) || defined(__STATS) || defined(__PHASE_SYNC)
	// Create compound datatype for the complex datasets
	file_info->COMPLEX_DTYPE = CreateComplexDatatype();
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
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "ZMinus", &(file_info->file_space[DSET_MINUS]), &(file_info->data_set[DSET_MINUS]), &(file_info->mem_space[DSET_MINUS]), file_info->COMPLEX_DTYPE, dims2D, maxdims2D, chunkdims2D, Dim2D);
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
	// Energy Spectrum
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "EnergySpectrum", &(file_info->file_space[DSET_ENRG_SPECT]), &(file_info->data_set[DSET_ENRG_SPECT]), &(file_info->mem_space[DSET_ENRG_SPECT]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);	
	#endif
	#if defined(__DISS_SPECT)
	// Dissipation Spectrum
	CreateSlabbedDSet(0.0, 0, file_info->output_file_handle, "DissipationSpectrum", &(file_info->file_space[DSET_DISS_SPECT]), &(file_info->data_set[DSET_DISS_SPECT]), &(file_info->mem_space[DSET_DISS_SPECT]), H5T_NATIVE_DOUBLE, dims2D, maxdims2D, chunkdims2D, Dim2D);
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
		sprintf(sys_type, "%s", "ELSASSAR_MHD");
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
		sprintf(sys_type, "%s", "ELSASSAR_MHD");
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
	if ( (H5LTread_dataset(file_info->stats_file_handle, "VelAmps", H5T_NATIVE_DOUBLE, tmp_u_amp)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelAmps");
	}
	double* tmp_u_phase = (double* )malloc(sizeof(double) * N);
	if ( (H5LTread_dataset(file_info->stats_file_handle, "VelPhases", H5T_NATIVE_DOUBLE, tmp_u_phase)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelPhases");
	}
	
	// Write tmp array to velocity modes / mode amplitudes and phases
	for (int i = 0; i < N; ++i) {
		n = i + 2;
		run_data->a_n[n]   = tmp_u_amp[i];
		run_data->phi_n[n] = tmp_u_phase[i];
	}

	// Free temp memory
	free(tmp_u_amp);
	free(tmp_u_phase);

	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	///----------------------------- Read in initial magnetic amps and phases
	// Create tmp array to read in data
	double* tmp_b_amp = (double* )malloc(sizeof(double) * N);
	if ( (H5LTread_dataset(file_info->stats_file_handle, "MagAmps", H5T_NATIVE_DOUBLE, tmp_b_amp)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagAmps");
	}
	double* tmp_b_phase = (double* )malloc(sizeof(double) * N);
	if ( (H5LTread_dataset(file_info->stats_file_handle, "MagPhases", H5T_NATIVE_DOUBLE, tmp_b_phase)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagPhases");
	}
	
	// Write tmp array to magnetic modes / mode amplitudes and phases
	for (int i = 0; i < N; ++i) {
		n = i + 2;
		run_data->b_n[n]   = tmp_b_amp[i];
		run_data->psi_n[n] = tmp_b_phase[i];
	}

	// Free temp memory
	free(tmp_b_amp);
	free(tmp_b_phase);
	#endif
	#else
	///----------------------------- Read in initial velocity modes
	// Create tmp array to read in data
	double complex* tmp_u = (double complex* )malloc(sizeof(double complex) * N);
	if ( (H5LTread_dataset(file_info->stats_file_handle, "VelModes", H5T_NATIVE_DOUBLE, tmp_u)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelModes");
	}
	
	// Write tmp_u array to velocity modes / mode amplitudes and phases
	for (int i = 0; i < N; ++i) {
		n = i + 2;
		run_data->u[n] = tmp_u[i];
		#if defined(PHASE_ONLY_FXD_AMP)
		run_data->a_n[n]   = cabs(tmp_u[i]);
		run_data->phi_n[n] = carg(tmp_u[i]);
		#endif
	}

	// Free temp memory
	free(tmp_u);

	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	///----------------------------- Read in initial magnetic modes
	// Create tmp array to read in data
	double complex* tmp_b = (double complex* )malloc(sizeof(double complex) * N);
	if ( (H5LTread_dataset(file_info->stats_file_handle, "MagModes", H5T_NATIVE_DOUBLE, tmp_b)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagModes");
	}
	
	// Write tmp_b array to magnetic modes / mode amplitudes and phases
	for (int i = 0; i < N; ++i) {
		n = i + 2;
		run_data->b[n] = tmp_b[i];
		#if defined(PHASE_ONLY_FXD_AMP)
		run_data->b_n[n]   = cabs(tmp_b[i]);
		run_data->psi_n[n] = carg(tmp_b[i]);
		#endif
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

	// Initialize variables
	herr_t status;
	static const hsize_t D1 = 1;
	hsize_t dims1D[D1];
	#if defined(STATS) || defined(__PHASE_SYNC_STATS)
	const hsize_t D2 = 2;
	hsize_t dims2D[D2];
	const hsize_t D3 = 3;
	hsize_t dims3D[D3];
	#endif

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
	// Write Wavenumbers
	// -------------------------------
	#if defined(__WAVELIST)
	// Allocate array to gather the wavenumbers from each of the local arrays - in the x direction
	dims1D[0] = N;
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "k", D1, dims1D, H5T_NATIVE_DOUBLE, &(run_data->k[2]))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "k");
	}
	#endif

	// -------------------------------
	// Write System Measures
	// -------------------------------
	// Time
	#if defined(__TIME)
	dims1D[0] = sys_vars->num_print_steps;
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "Time", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->time)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "Time");
	}
	#endif

	///----------------- System Measures
	#if defined(__SYS_MEASURES)
	// Energy
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TotalEnergy", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->tot_energy)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TotalEnergy");
	}
	// Velocity Helicity
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TotalVelocityHelicity", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->tot_hel_u)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TotalVelocityHelicity");
	}
	// Total Dissipation
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TotalVelocityDissipation", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->tot_diss_u)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TotalVelocityDissipation");
	}
	// Characteristic Velocity
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "CharacteristicVel", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->u_charact)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "CharacteristicVel");
	}
	// Integral Length Scale
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "IntegralLengthScale", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->int_scale)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "IntegralLengthScale");
	}
	// Reynolds no.
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "ReynoldsNo", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->reynolds_no)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "ReynoldsNo");
	}
	// Taylor Miscroscale
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TaylorMicroScale", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->taylor_micro_scale)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TaylorMicroScale");
	}
	// Kolmogorov Length Scale
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "KolmogorovLengthScale", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->kolmogorov_scale)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "KolmogorovLengthScale");
	}
	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	// Magnetic Helicity
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TotalMagneticHelicity", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->tot_hel_b)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TotalMagneticHelicity");
	}
	// Cross Helicity
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TotalCrossHelicity", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->tot_cross_hel)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TotalCrossHelicity");
	}
	// Total Dissipation
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TotalMagneticDissipation", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->tot_diss_b)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TotalMagneticDissipation");
	}
	#endif
	#endif


    ///----------------- Time Averaged 
	#if defined(__VEL_AMP_AVG)
	// Time average the data
	for (int i = 0; i < sys_vars->N; ++i) {
		run_data->a_n_t_avg[i] /= run_data->num_sys_msr_steps;
	}

	dims1D[0] = sys_vars->N;
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TimeAveragedVelocityAmplitudes", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->a_n_t_avg)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TimeAveragedVelocityAmplitudes");
	}
	#endif
	#if defined(__MAG_AMP_AVG) && (defined(__MAGNETO) || defined(__ELSASSAR_MHD))
	// Time average the data
	for (int i = 0; i < sys_vars->N; ++i) {
		run_data->b_n_t_avg[i] /= run_data->num_sys_msr_steps;
	}

	dims1D[0] = sys_vars->N;
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TimeAveragedMagneticAmplitudes", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->b_n_t_avg)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TimeAveragedMagneticAmplitudes");
	}
	#endif

	///----------------- Time Averaged Energy Spectrum
	#if defined(__ENRG_SPECT_AVG)
	// Time average the data
	for (int i = 0; i < sys_vars->N; ++i) {
		run_data->energy_spect_t_avg[i] /= run_data->num_sys_msr_steps;
	}

	dims1D[0] = sys_vars->N;
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TimeAveragedEnergySpectrum", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->energy_spect_t_avg)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TimeAveragedEnergySpectrum");
	}
	#endif

	///----------------- Time Averaged Dissipation Spectrum
	#if defined(__DISS_SPECT_AVG)
	// Time average the data
	for (int i = 0; i < sys_vars->N; ++i) {
		run_data->diss_spect_t_avg[i] /= run_data->num_sys_msr_steps;
	}

	dims1D[0] = sys_vars->N;
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TimeAveragedDissipationSpectrum", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->diss_spect_t_avg)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TimeAveragedDissipationSpectrum");
	}
	#endif

	///----------------- Time Averaged Dissipation Spectrum
	#if defined(__ENRG_FLUX_AVG)
	// Time average the data
	for (int i = 0; i < sys_vars->N; ++i) {
		run_data->energy_flux_t_avg[i] /= run_data->num_sys_msr_steps;
	}

	dims1D[0] = sys_vars->N;
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TimeAveragedEnergyFlux", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->energy_flux_t_avg)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TimeAveragedEnergyFlux");
	}
	#endif

	///----------------- Total Energy Variation Variables
	#if defined(__TOT_ENRG_FLUX)
	// Energy Flux
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TotalEnergyFlux", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->tot_energy_flux)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TotalEnergyFlux");
	}
	// Velocity Energy Dissipation
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TotalVelEnergyDiss", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->tot_energy_diss_u)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TotalVelEnergyDiss");
	}
	// Velocity Energy Input
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TotalVelEnergyInput", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->tot_energy_input_u)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TotalVelEnergyInput");
	}
	#if defined(__MAGNETO)
	// Magnetic Energy Dissipation
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TotalMagEnergyDiss", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->tot_energy_diss_b)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TotalMagEnergyDiss");
	}
	// Magnetic Energy Input
	if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "TotalMagEnergyInput", D1, dims1D, H5T_NATIVE_DOUBLE, run_data->tot_energy_input_b)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "TotalMagEnergyInput");
	}
	#endif
	#endif

	// -------------------------------
	// Write Forcing
	// -------------------------------
	#if defined(__FORCING)
	if(strcmp(sys_vars->forcing, "STOC") != 0) {
		dims1D[0] = N;

		// Velocity Forcing
		if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "VelocityForcing", D1, dims1D, file_info->COMPLEX_DTYPE, &(run_data->forcing_u[1]))) < 0) {
			printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelocityForcing");
		}
		#if defined(__MAGNETO)
		// Magnetic Forcing
		if ( (H5LTmake_dataset(file_info->sys_msr_file_handle, "MagneticForcing", D1, dims1D, file_info->COMPLEX_DTYPE, &(run_data->forcing_b[1]))) < 0) {
			printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagneticForcing");
		}
		#endif
	}
	#endif

	// -------------------------------
	// Write Phase Sync Stats
	// -------------------------------
	#if defined(__PHASE_SYNC_STATS)
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
	dims1D[0] = VEL_NUM_BINS + 1;
	if ( (H5LTmake_dataset(file_info->phase_sync_file_handle, "VelTriads_Ranges", D1, dims1D, H5T_NATIVE_DOUBLE, vel_triads_ranges)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelTriads_Ranges");
	}

	// Free temporary memory
	free(vel_triads_counts);
	free(vel_triads_ranges);

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
	#endif
	
	// -------------------------------
	// Write Stats
	// -------------------------------
	#if defined(STATS)
	///-------------------------- Write velocity field 
	#if defined(PHASE_ONLY)
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "VelAmps", D1, dims1D, H5T_NATIVE_DOUBLE, &(run_data->a_n[2]))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelAmps");
	}
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "VelPhases", D1, dims1D, H5T_NATIVE_DOUBLE, &(run_data->phi_n[2]))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelPhases");
	}
	#else
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "VelModes", D1, dims1D, file_info->COMPLEX_DTYPE, &(run_data->u[2]))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelModes");
	}
	#endif

	///-------------------------- Velocity field stats
	double* tmp_vel_stats = (double* )malloc(sizeof(double) * (NUM_RUN_STATS) * N);
	for (int i = 0; i < N; ++i) {
		tmp_vel_stats[i * (NUM_RUN_STATS) + 0] = gsl_rstat_mean(stats_data->real_vel_moments[i]);
		tmp_vel_stats[i * (NUM_RUN_STATS) + 1] = gsl_rstat_variance(stats_data->real_vel_moments[i]);
		tmp_vel_stats[i * (NUM_RUN_STATS) + 2] = gsl_rstat_skew(stats_data->real_vel_moments[i]);
		tmp_vel_stats[i * (NUM_RUN_STATS) + 3] = gsl_rstat_kurtosis(stats_data->real_vel_moments[i]);		
		tmp_vel_stats[i * (NUM_RUN_STATS) + 4] = gsl_rstat_rms(stats_data->real_vel_moments[i]);
		tmp_vel_stats[i * (NUM_RUN_STATS) + 5] = gsl_rstat_min(stats_data->real_vel_moments[i]);
		tmp_vel_stats[i * (NUM_RUN_STATS) + 6] = gsl_rstat_max(stats_data->real_vel_moments[i]);
	}

	// Write data 
	dims2D[0] = N;
	dims2D[1] = NUM_RUN_STATS;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "VelStats", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_vel_stats)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "VelStats");
	}

	// Free temp memory
	free(tmp_vel_stats);


	///-------------------------- Real Velocity Histogram
	#if defined(__VEL_HIST)
	double* tmp_vel_hist_bin    = (double* )malloc(sizeof(double) * (VEL_NUM_BINS) * N);
	double* tmp_vel_hist_ranges = (double* )malloc(sizeof(double) * (VEL_NUM_BINS + 1) * N);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < VEL_NUM_BINS + 1; ++j) {
			if (j < VEL_NUM_BINS) {
				tmp_vel_hist_bin[i * VEL_NUM_BINS + j] = stats_data->real_vel_hist[i]->bin[j];
			}
			tmp_vel_hist_ranges[i * (VEL_NUM_BINS + 1) + j] = stats_data->real_vel_hist[i]->range[j];
		}
	}

	// Write Bin Count data 
	dims2D[0] = N;
	dims2D[1] = VEL_NUM_BINS;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "RealVelHist_Counts", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_vel_hist_bin)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "RealVelHist_Counts");
	}
	// Write Bin Range data 
	dims2D[0] = N;
	dims2D[1] = VEL_NUM_BINS + 1;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "RealVelHist_Ranges", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_vel_hist_ranges)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "RealVelHist_Ranges");
	}

	// Free temporary memory
	free(tmp_vel_hist_bin);
	free(tmp_vel_hist_ranges);
	#endif

	///-------------------------- Velocity structure function
	#if defined(__STR_FUNC_VEL)
	// Allocate temporary contiguous array
	double* tmp_vel_str = (double* )malloc(sizeof(double) * NUM_POW * N);
	for (int i = 0; i < N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_vel_str[i * NUM_POW + p] = stats_data->vel_str_func[p][i] / stats_data->num_stats_steps;
		}
	}

	// Write data 
	dims2D[0] = N;
	dims2D[1] = NUM_POW;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionVel", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_vel_str)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionVel");
	}

	// Free temp memory
	free(tmp_vel_str);
	#endif

	///-------------------------- Velocity Flux structure function
	#if defined(__STR_FUNC_VEL_FLUX)
	// Allocate temporary contiguous array
	double* tmp_vel_str_flux = (double* )malloc(sizeof(double) * 2 * NUM_POW * N);
	for (int i = 0; i < N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_vel_str_flux[2 * (i * NUM_POW + p) + 0] = stats_data->vel_flux_str_func[0][p][i] / stats_data->num_stats_steps;
			tmp_vel_str_flux[2 * (i * NUM_POW + p) + 1] = stats_data->vel_flux_str_func[1][p][i] / stats_data->num_stats_steps;
		}
	}

	// Write data 
	dims3D[0] = N;
	dims3D[1] = NUM_POW;
	dims3D[2] = 2;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionVelFlux", D3, dims3D, H5T_NATIVE_DOUBLE, tmp_vel_str_flux)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionVelFlux");
	}

	for (int i = 0; i < N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_vel_str_flux[2 * (i * NUM_POW + p) + 0] = stats_data->vel_flux_str_func_abs[0][p][i] / stats_data->num_stats_steps;
			tmp_vel_str_flux[2 * (i * NUM_POW + p) + 1] = stats_data->vel_flux_str_func_abs[1][p][i] / stats_data->num_stats_steps;
		}
	}

	// Write data 
	dims3D[0] = N;
	dims3D[1] = NUM_POW;
	dims3D[2] = 2;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionVelFluxAbs", D3, dims3D, H5T_NATIVE_DOUBLE, tmp_vel_str_flux)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionVelFluxAbs");
	}

	// Free temp memory
	free(tmp_vel_str_flux);
	#endif


	#if defined(__MAGNETO) || defined(__ELSASSAR_MHD)
	///-------------------------- Write Magnetic field 
	#if defined(PHASE_ONLY)
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "MagAmps", D1, dims1D, H5T_NATIVE_DOUBLE, &(run_data->b_n[2]))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagAmps");
	}
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "MagPhases", D1, dims1D, H5T_NATIVE_DOUBLE, &(run_data->psi_n[2]))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagPhases");
	}
	#else
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "MagModes", D1, dims1D, file_info->COMPLEX_DTYPE, &(run_data->b[2]))) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagModes");
	}
	#endif

	///-------------------------- Real Velocity Histogram
	#if defined(__MAG_HIST)
	double* tmp_mag_hist_bin    = (double* )malloc(sizeof(double) * (VEL_NUM_BINS) * N);
	double* tmp_mag_hist_ranges = (double* )malloc(sizeof(double) * (VEL_NUM_BINS + 1) * N);
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < VEL_NUM_BINS + 1; ++j) {
			if (j < VEL_NUM_BINS) {
				tmp_mag_hist_bin[i * VEL_NUM_BINS + j] = stats_data->real_mag_hist[i]->bin[j];
			}
			tmp_mag_hist_ranges[i * (VEL_NUM_BINS + 1) + j] = stats_data->real_mag_hist[i]->range[j];
		}
	}

	// Write Bin Count data 
	dims2D[0] = N;
	dims2D[1] = VEL_NUM_BINS;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "RealMagHist_Counts", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_mag_hist_bin)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "RealMagHist_Counts");
	}
	// Write Bin Range data 
	dims2D[0] = N;
	dims2D[1] = VEL_NUM_BINS + 1;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "RealMagHist_Ranges", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_mag_hist_ranges)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "RealMagHist_Ranges");
	}

	// Free temporary memory
	free(tmp_mag_hist_bin);
	free(tmp_mag_hist_ranges);
	#endif

	///-------------------------- Magnetic field stats
	double* tmp_mag_stat = (double* )malloc(sizeof(double) * (NUM_RUN_STATS) * N);
	for (int i = 0; i < N; ++i) {
		tmp_mag_stat[i * (NUM_RUN_STATS) + 0] = gsl_rstat_mean(stats_data->real_mag_moments[i]);
		tmp_mag_stat[i * (NUM_RUN_STATS) + 1] = gsl_rstat_variance(stats_data->real_mag_moments[i]);
		tmp_mag_stat[i * (NUM_RUN_STATS) + 2] = gsl_rstat_skew(stats_data->real_mag_moments[i]);
		tmp_mag_stat[i * (NUM_RUN_STATS) + 3] = gsl_rstat_kurtosis(stats_data->real_mag_moments[i]);		
		tmp_mag_stat[i * (NUM_RUN_STATS) + 4] = gsl_rstat_rms(stats_data->real_mag_moments[i]);
		tmp_mag_stat[i * (NUM_RUN_STATS) + 5] = gsl_rstat_min(stats_data->real_mag_moments[i]);
		tmp_mag_stat[i * (NUM_RUN_STATS) + 6] = gsl_rstat_max(stats_data->real_mag_moments[i]);
	}

	// Write data 
	dims2D[0] = N;
	dims2D[1] = NUM_RUN_STATS;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "MagStats", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_mag_stat)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "MagStats");
	}

	// Free temp memory
	free(tmp_mag_stat);


	///-------------------------- Magnetic structure function
	#if defined(__STR_FUNC_MAG) 
	// Allocate temporary contiguous array
	double* tmp_mag_str = (double* )malloc(sizeof(double) * NUM_POW * N);
	for (int i = 0; i < N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_mag_str[i * NUM_POW + p] = stats_data->mag_str_func[p][i] / stats_data->num_stats_steps;
		}
	}

	// Write data 
	dims2D[0] = N;
	dims2D[1] = NUM_POW;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionMag", D2, dims2D, H5T_NATIVE_DOUBLE, tmp_mag_str)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionMag");
	}

	// Free temp memory
	free(tmp_mag_str);
	#endif

	///-------------------------- Magnetic flux structure function
	#if defined(__STR_FUNC_MAG_FLUX)
	// Allocate temporary contiguous array
	double* tmp_mag_str_flux = (double* )malloc(sizeof(double) * 2 * NUM_POW * N);
	for (int i = 0; i < N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_mag_str_flux[2 * (i * NUM_POW + p) + 0] = stats_data->mag_flux_str_func[0][p][i] / stats_data->num_stats_steps;
			tmp_mag_str_flux[2 * (i * NUM_POW + p) + 1] = stats_data->mag_flux_str_func[1][p][i] / stats_data->num_stats_steps;
		}
	}

	// Write data 
	dims3D[0] = N;
	dims3D[1] = NUM_POW;
	dims3D[2] = 2;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionMagFlux", D3, dims3D, H5T_NATIVE_DOUBLE, tmp_mag_str_flux)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionMagFlux");
	}

	for (int i = 0; i < N; ++i) {
		for (int p = 0; p < NUM_POW; ++p) {
			tmp_mag_str_flux[2 * (i * NUM_POW + p) + 0] = stats_data->mag_flux_str_func_abs[0][p][i] / stats_data->num_stats_steps;
			tmp_mag_str_flux[2 * (i * NUM_POW + p) + 1] = stats_data->mag_flux_str_func_abs[1][p][i] / stats_data->num_stats_steps;
		}
	}

	// Write data 
	dims3D[0] = N;
	dims3D[1] = NUM_POW;
	dims3D[2] = 2;
	if ( (H5LTmake_dataset(file_info->stats_file_handle, "StructureFunctionMagFluxAbs", D3, dims3D, H5T_NATIVE_DOUBLE, tmp_mag_str_flux)) < 0) {
		printf("\n["MAGENTA"WARNING"RESET"] --- Failed to make dataset ["CYAN"%s"RESET"]\n", "StructureFunctionMagFluxAbs");
	}

	// Free temp memory
	free(tmp_mag_str_flux);
	#endif
	#endif
	#endif
	
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

	#if defined(__VEL) || defined(__MAG) defined(__Z_PLUS) || defined(__Z_MINUS) || defined(__FORCING) || defined(__STATS) || defined(__PHASE_SYNC)
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