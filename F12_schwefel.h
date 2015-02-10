#ifndef _F12_schwefel
#define _F12_schwefel
#include"test_func.h"
class F12_schwefel:public test_func {
private:
	// Shifted global optimum
	double *m_o;
	double **m_a;
	double **m_b;
	
	// In order to avoid excessive memory allocation,
	// a fixed memory buffer is allocated for each function object.
	double *m_A;
	double *m_B;
public:

	// Constructors
	F12_schwefel (int dimension, double bias);
	F12_schwefel (int dimension, double bias, string file_data);

	// Function body
	double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME12 = "Schwefel's Problem 2.13";
const string DEFAULT_FILE_DATA12 = "supportData/schwefel_213_data.txt";
#endif
