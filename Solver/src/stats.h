/**
* @file stats.h
* @author Enda Carroll
* @date Sept 2022
* @brief Header file containging funciton prototypes for the stats.c file
*/
// ---------------------------------------------------------------------
//  Standard Libraries and Headers
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
//  User Libraries and Headers
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
//  Function Definitions
// ---------------------------------------------------------------------
void InitializeStats(void);
void ComputeStats(const long int iters, const long int save_data_indx);
void WriteStatsToFile(void);
void FreeStatsObjects(void);
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------