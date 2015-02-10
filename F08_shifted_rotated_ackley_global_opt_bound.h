#ifndef _F08_shifted_rotated_ackley_global_opt_bound
#define _F08_shifted_rotated_ackley_global_opt_bound
#include"test_func.h"
class F08_shifted_rotated_ackley_global_opt_bound:public test_func {
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
	F08_shifted_rotated_ackley_global_opt_bound (int dimension, double bias);
	F08_shifted_rotated_ackley_global_opt_bound (int dimension, double bias, string file_data, string file_m);

	// Function body
	double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME8 = "Shifted Rotated Ackley's Function with Global Optimum on Bounds";
const string DEFAULT_FILE_DATA8 = "supportData/ackley_func_data.txt";
const string DEFAULT_FILE_MX_PREFIX8 = "supportData/ackley_M_D";
const string DEFAULT_FILE_MX_SUFFIX8 = ".txt";
#endif
