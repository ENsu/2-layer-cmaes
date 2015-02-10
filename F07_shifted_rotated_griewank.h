#ifndef _F07_shifted_rotated_griewank
#define _F07_shifted_rotated_griewank
#include"test_func.h"
class F07_shifted_rotated_griewank:public test_func {
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
	F07_shifted_rotated_griewank (int dimension, double bias);
	F07_shifted_rotated_griewank (int dimension, double bias, string file_data, string file_m);

	// Function body
	double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME7 = "Shifted Rotated Griewank's Function without Bounds";
const string DEFAULT_FILE_DATA7 = "supportData/griewank_func_data.txt";
const string DEFAULT_FILE_MX_PREFIX7 = "supportData/griewank_M_D";
const string DEFAULT_FILE_MX_SUFFIX7 = ".txt";
#endif
