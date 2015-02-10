#ifndef _F17_rotated_hybrid_composition_1_noise
#define _F17_rotated_hybrid_composition_1_noise
#include"test_func.h"
#include"HCJob.h"
#include"benchmark.h"
class MyHCJob17:public HCJob {
public:
	double basic_func(int func_no, double *x,int length);
};
class F17_rotated_hybrid_composition_1_noise:public test_func {
private:
	MyHCJob17 theJob;

	// Shifted global optimum
	double **m_o;
	double ***m_M;
	double *m_testPoint;
	double *m_testPointM;
	double *m_fmax;

	// In order to avoid excessive memory allocation,
	// a fixed memory buffer is allocated for each function object.
	double *m_w;
	double **m_z;
	double **m_zM;
public:
	// Number of functions
	static const int NUM_FUNC = 10;
	// Constructors
	F17_rotated_hybrid_composition_1_noise (int dimension, double bias);
	F17_rotated_hybrid_composition_1_noise (int dimension, double bias, string file_data, string file_m);

	// Function body
	double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME17 = "Rotated Hybrid Composition Function 1 with Noise in Fitness";
const string DEFAULT_FILE_DATA17 = "supportData/hybrid_func1_data.txt";
const string DEFAULT_FILE_MX_PREFIX17 = "supportData/hybrid_func1_M_D";
const string DEFAULT_FILE_MX_SUFFIX17 = ".txt";

static double m_sigma17[] = {
	1.0,    1.0,    1.0,    1.0,    1.0,
	1.0,    1.0,    1.0,    1.0,    1.0
};
static double m_lambda17[] = {
	1.0,            1.0,            10.0,           10.0,
	5.0/60.0,       5.0/60.0,       5.0/32.0,       5.0/32.0,
	5.0/100.0,      5.0/100.0
};
static double m_func_biases17[] = {
	0.0,    100.0,  200.0,  300.0,  400.0,
	500.0,  600.0,  700.0,  800.0,  900.0
};
#endif
