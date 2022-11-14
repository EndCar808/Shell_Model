/**
* @file hdf5_funcs.h  
* @author Enda Carroll
* @date Sept 2022
* @brief File containing function prototpyes for hdf5_funcs file
*/
// ---------------------------------------------------------------------
//  Standard Libraries and Headers
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
//  User Libraries and Headers
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
//  Function Prototpyes
// ---------------------------------------------------------------------
void CreateOutputFilesWriteICs(const long int N);
void GetOutputDirPath(void);
hid_t CreateComplexDatatype(void);
void ReadInputFile(const long int N);
void CreateSlabbedDSet(double t, int iters, hid_t file_handle, char* dset_name, hid_t* file_space, hid_t* data_set, hid_t* mem_space, hid_t dtype, hsize_t* dset_dims, hsize_t* dset_max_dims, hsize_t* dset_chunk_dims, const int num_dims);
void WriteSlabbedDataReal(double t, int iters, hid_t file_space, hid_t data_set, hid_t mem_space, hid_t dtype, double* data, char* dset_name, hsize_t* count, hsize_t* index);
void WriteSlabbedDataFourier(double t, int iters, hid_t file_space, hid_t data_set, hid_t mem_space, hid_t dtype, double complex* data, char* dset_name, hsize_t* count, hsize_t* index);
void WriteDataToFile(double t, const long int iters, const long int save_indx);
void FinalWriteAndCloseOutputFile(const long int N, int iters, int save_data_indx);
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------