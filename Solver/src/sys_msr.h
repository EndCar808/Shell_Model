/**
* @file utils.h
* @author Enda Carroll
* @date Sept 2022
* @brief Header file for the sys_msr.c file
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
void InitializeSystemMeasurables(RK_data_struct* RK_data);
void ComputeSystemMeasurables(double t, const long int iter, const long int save_iter, RK_data_struct* RK_data);
void WriteSystemMeasuresToFile(void);
void FreeSystemMeasuresObjects(void);
// ---------------------------------------------------------------------
//  End of File
// ---------------------------------------------------------------------