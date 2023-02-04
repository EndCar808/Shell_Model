/**
* @file stats.h
* @author Enda Carroll
* @date Sept 2022
* @brief Header file containging funciton prototypes for the phase_sync.c file
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
void ComputePhaseSyncData(const long int iter);
void InitializePhaseSyncObjects(void);
void FreePhaseSyncObjects(void);
void WritePhaseSyncStatsToFile(void);
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------