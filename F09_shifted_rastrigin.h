#ifndef _F09_shifted_rastrigin
#define _F09_shifted_rastrigin
#include"test_func.h"
class F09_shifted_rastrigin:public test_func {
private:
	// Shifted global optimum
	double *m_o;

	// In order to avoid excessive memory allocation,
	// a fixed memory buffer is allocated for each function object.
	double *m_z;
public:
	// Constructors
	F09_shifted_rastrigin (int dimension, double bias);
	F09_shifted_rastrigin (int dimension, double bias, string file_data);

	// Function body
	double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME9 = "Shifted Rastrigin's Function";
const string DEFAULT_FILE_DATA9 = "supportData/rastrigin_func_data.txt";
#endif
