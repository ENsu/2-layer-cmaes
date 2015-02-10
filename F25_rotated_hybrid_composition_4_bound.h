#ifndef _F25_rotated_hybrid_composition_4_bound
#define _F25_rotated_hybrid_composition_4_bound
#include"test_func.h"
#include"HCJob.h"
#include"benchmark.h"
class MyHCJob25:public HCJob {
	public:
		double basic_func(int func_no, double *x,int length);
};
class F25_rotated_hybrid_composition_4_bound:public test_func {
private:
	MyHCJob25 theJob;

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
	F25_rotated_hybrid_composition_4_bound (int dimension, double bias);
	F25_rotated_hybrid_composition_4_bound (int dimension, double bias, string file_data, string file_m);

	// Function body
	double f(double* x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME25 = "Rotated Hybrid Composition Function 4 without bounds";
const string DEFAULT_FILE_DATA25 = "supportData/hybrid_func4_data.txt";
const string DEFAULT_FILE_MX_PREFIX25 = "supportData/hybrid_func4_M_D";
const string DEFAULT_FILE_MX_SUFFIX25 = ".txt";
static double m_sigma25[] = {
	2.0,    2.0,    2.0,    2.0,    2.0,
	2.0,    2.0,    2.0,    2.0,    2.0
};
static double m_lambda25[] = {
	10.0,           5.0/20.0,       1.0,    5.0/32.0,       1.0,
	5.0/100.0,      5.0/50.0,       1.0,    5.0/100.0,      5.0/100.0
};
static double m_func_biases25[] = {
	0.0,    100.0,  200.0,  300.0,  400.0,
	500.0,  600.0,  700.0,  800.0,  900.0
};
#endif
