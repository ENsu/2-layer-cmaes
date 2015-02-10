#ifndef _F14_shifted_rotated_expanded_scaffer
#define _F14_shifted_rotated_expanded_scaffer
#include"test_func.h"
class F14_shifted_rotated_expanded_scaffer:public test_func {
private:
	// Shifted global optimum
	double *m_o;
	double **m_matrix;

	// In order to avoid excessive memory allocation,
	// a fixed memory buffer is allocated for each function object.
	double *m_z;
	double *m_zM;
public:
	// Constructors
	F14_shifted_rotated_expanded_scaffer (int dimension, double bias);
	F14_shifted_rotated_expanded_scaffer (int dimension, double bias, string file_data, string file_m);

	// Function body
	double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME14 = "Shifted Rotated Expanded Scaffer's F6 Function";
const string DEFAULT_FILE_DATA14 = "supportData/E_ScafferF6_func_data.txt";
const string DEFAULT_FILE_MX_PREFIX14 = "supportData/E_ScafferF6_M_D";
const string DEFAULT_FILE_MX_SUFFIX14 = ".txt";
#endif
