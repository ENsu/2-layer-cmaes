#ifndef _F05_schwefel_global_opt_bound
#define _F05_schwefel_global_opt_bound
#include"test_func.h"
class F05_schwefel_global_opt_bound:public test_func {
private:
	// Shifted global optimum
	double *m_o;
	double **m_A;

	// In order to avoid excessive memory allocation,
	// a fixed memory buffer is allocated for each function object.
	double *m_B;
	double *m_z;
public:
	// Constructors
	F05_schwefel_global_opt_bound (int dimension, double bias);
	F05_schwefel_global_opt_bound (int dimension, double bias, string file_data);

	// Function body
	double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME5 = "Schwefel's Problem 2.6 with Global Optimum on Bounds";
const string DEFAULT_FILE_DATA5 = "supportData/schwefel_206_data.txt";
#endif
