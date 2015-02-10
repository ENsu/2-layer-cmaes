#ifndef _F06_shifted_rosenbrock
#define _F06_shifted_rosenbrock
#include"test_func.h"

class F06_shifted_rosenbrock:public test_func {
private:
	// Shifted global optimum
	double *m_o;

	// In order to avoid excessive memory allocation,
	// a fixed memory buffer is allocated for each function object.
	double *m_z;

public:
	// Constructors
	F06_shifted_rosenbrock (int dimension, double bias);
	F06_shifted_rosenbrock (int dimension, double bias, string file_data);

	// Function body
	double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME6 = "Shifted Rosenbrock's Function";
const string DEFAULT_FILE_DATA6 = "supportData/rosenbrock_func_data.txt";
#endif
