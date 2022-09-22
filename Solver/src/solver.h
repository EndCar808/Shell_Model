/**
* @file solver.h
* @author Enda Carroll
* @date Sept 2022
* @brief Header file containing the function prototypes for the solver.c file
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
// Main function for the solver
void Solve(void);
// Integration functions
#if defined(__INT_FAC_RK4)
void IntFacRK4Step(const double dt, const long int N, RK_data_struct* RK_data);
#endif
void NonlinearTerm(fftw_complex* u, fftw_complex* b, fftw_complex* u_nonlin, fftw_complex* b_nonlin, const long int N);
// Initialize the system functions
void InitializeShellWavenumbers(long int* k, const long int N);
void InitialConditions(const long int N);
void InitializeIntegrationVariables(double* t0, double* t, double* dt, double* T, long int* trans_steps);
// Memory Functions
void AllocateMemory(const long int N, RK_data_struct* RK_data);
void FreeMemory(RK_data_struct* RK_data);
// Print Update
void PrintUpdateToTerminal(int iters, double t, double dt, double T, int save_data_indx);
// Check System
void SystemCheck(double dt, int iters);
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------