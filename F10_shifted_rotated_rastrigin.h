#ifndef _F10_shifted_rotated_rastrigin
#define _F10_shifted_rotated_rastrigin
#include"test_func.h"
class F10_shifted_rotated_rastrigin:public test_func {
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
	F10_shifted_rotated_rastrigin (int dimension, double bias);
	F10_shifted_rotated_rastrigin (int dimension, double bias, string file_data, string file_m);

	// Function body
	double f(double* x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME10 = "Shifted Rotated Rastrigin's Function";
const string DEFAULT_FILE_DATA10 = "supportData/rastrigin_func_data.txt";
const string DEFAULT_FILE_MX_PREFIX10 = "supportData/rastrigin_M_D";
const string DEFAULT_FILE_MX_SUFFIX10 = ".txt";
#endif
