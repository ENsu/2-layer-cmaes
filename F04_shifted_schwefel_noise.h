#ifndef _F04_shifted_schwefel_noise
#define _F04_shifted_schwefel_noise
#include"test_func.h"
class F04_shifted_schwefel_noise:public test_func {
private:
	// Shifted global optimum
	double *m_o;

	// In order to avoid excessive memory allocation,
	// a fixed memory buffer is allocated for each function object.
	double *m_z;
public:
	// Constructors
	F04_shifted_schwefel_noise (int dimension, double bias);
	F04_shifted_schwefel_noise (int dimension, double bias, string file_data);

	// Function body
	double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME4 = "Shifted Schwefel's Problem 1.2 with Noise in Fitness";
const string DEFAULT_FILE_DATA4 = "supportData/schwefel_102_data.txt";
	
#endif
