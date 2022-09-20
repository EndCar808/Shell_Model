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
#ifdef __RK4
void RK4Step(const double dt, const long int* N, const ptrdiff_t local_Nx, RK_data_struct* RK_data);
#elif defined(__RK5) || defined(__DPRK5)
void RK5DPStep(const double dt, const long int* N, const int iters, const ptrdiff_t local_Nx, RK_data_struct* RK_data);
#endif
#ifdef __DPRK5
double DPMax(double a, double b);
double DPMin(double a, double b);
#endif
// Initialize the system functions
void InitializeShellWavenumbers(long int* k, const long int N);
// Memory Functions
void AllocateMemory(const long int N, RK_data_struct* RK_data);
void FreeMemory(RK_data_struct* RK_data);

// void NonlinearRHSBatch(fftw_complex* w_hat, fftw_complex* dw_hat_dt, double* nonlinear, double* u, double* w);
// void ApplyDealiasing(fftw_complex* array, int array_dim, const long int* N);
// void ForceConjugacy(fftw_complex* w_hat, const long int* N);
// void InitializeSpaceVariables(double** x, int** k, const long int* N);
// void InitializeIntegrationVariables(double* t0, double* t, double* dt, double* T, long int* trans_steps);
// void InitialConditions(fftw_complex* w_hat, double* u, fftw_complex* u_hat, const long int* N);
// void InitializeFFTWPlans(const long int* N);
// double GetMaxData(char* dtype);
// // Timestep
// void GetTimestep(double* dt);
// // Check System
// void SystemCheck(double dt, int iters);
// // Print Update
// void PrintUpdateToTerminal(int iters, double t, double dt, double T, int save_data_indx);
// // Testing
// void TestTaylorGreenVortex(const double t, const long int* N, double* norms);
// void TaylorGreenSoln(const double t, const long int* N);
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------