#ifndef _F13_shifted_expanded_griewank_rosenbrock
#define _F13_shifted_expanded_griewank_rosenbrock
#include"test_func.h"
class F13_shifted_expanded_griewank_rosenbrock:public test_func {
private:
	// Shifted global optimum
	double *m_o;

	// In order to avoid excessive memory allocation,
	// a fixed memory buffer is allocated for each function object.
	double *m_z;
public:
	// Constructors
	F13_shifted_expanded_griewank_rosenbrock (int dimension, double bias);
	F13_shifted_expanded_griewank_rosenbrock (int dimension, double bias, string file_data);

	// Function body
	double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME13 = "Shifted Expanded Griewank's plus Rosenbrock's Function";
const string DEFAULT_FILE_DATA13 = "supportData/EF8F2_func_data.txt";
#endif
