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
#if defined(INT_FAC_RK4)
void IntFacRK4Step(const double dt, const long int N, RK_data_struct* RK_data);
#endif
#if defined(RK4)
void RK4Step(const double dt, const long int N, RK_data_struct* RK_data);
#endif
#if defined(AB4CN)
void RK4Step(const double dt, const long int N, RK_data_struct* RK_data);
void AB4CNStep(const double dt, const long int iters, const long int N, RK_data_struct* RK_data);
#endif
#if defined(PHASE_ONLY_DIRECT) && !defined(PHASE_ONLY)
void NonlinearTerm(double* u, double* b, double* u_nonlin, double* b_nonlin, const long int N);
#else
void NonlinearTerm(double complex* u, double complex* b, double complex* u_nonlin, double complex* b_nonlin, const long int N);
#endif
// Initialize the system functions
void InitializeShellWavenumbers(double* k, const long int N);
void InitialConditions(const long int N);
void InitializeIntegrationVariables(double* t0, double* t, double* dt, double* T, long int* trans_steps);
void InitializeForicing(const long int N);
// Forcing Functions
void ComputeForicing(const long int N);
// Memory Functions
void AllocateMemory(const long int N, RK_data_struct* RK_data);
void FreeMemory(RK_data_struct* RK_data);
// Print Update
void PrintUpdateToTerminal(long int iters, double t, double dt, double T, long int save_data_indx);
// Check System
void SystemCheck(double dt, long int iters, long int save_data_indx);
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------