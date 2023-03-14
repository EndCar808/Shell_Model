/**
* @file utils.h
* @author Enda Carroll
* @date Sept 2022
* @brief Header file for the utils.c file
*/
// ---------------------------------------------------------------------
//  Standard Libraries and Headers
// ---------------------------------------------------------------------


// ---------------------------------------------------------------------
//  User Libraries and Headers
// ---------------------------------------------------------------------



// ---------------------------------------------------------------------
//  Function Prototypes
// ---------------------------------------------------------------------
// Command  Line Arguments
int GetCMLArgs(int argc, char** argv);
// Simulation Details
void PrintSimulationDetails(int argc, char** argv, double sim_time);
// Misc
double sgn(double x);
double log_lambda(double x);
double my_delta(double i, double j);
double my_mod_2pi(double phase);
void rng_amps(void);
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------