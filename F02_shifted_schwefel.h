#ifndef _F02_shifted_schwefel
#define _F02_shifted_schwefel
#include"test_func.h"
class F02_shifted_schwefel:public test_func{

private:
	// Shifted global optimum
	double *m_o;

	// In order to avoid excessive memory allocation,
	// a fixed memory buffer is allocated for each function object.
	double *m_z;

public:

	// Constructors
	F02_shifted_schwefel (int dimension, double bias);
	F02_shifted_schwefel (int dimension, double bias, string file_data);

	// Function body
	double f(double *x,int length);
};

// Fixed (class) parameters
const string FUNCTION_NAME2="Shifted Schwefel's Problem 1.2";
const string DEFAULT_FILE_DATA2="supportData/schwefel_102_data.txt";
#endif
