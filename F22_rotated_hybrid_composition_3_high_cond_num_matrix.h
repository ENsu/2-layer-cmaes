#ifndef _F22_rotated_hybrid_composition_3_high_cond_num_matrix
#define _F22_rotated_hybrid_composition_3_high_cond_num_matrix
#include"test_func.h"
#include"HCJob.h"
#include"benchmark.h"
class MyHCJob22:public HCJob {
public: 
	double basic_func(int func_no, double *x,int length);
};
class F22_rotated_hybrid_composition_3_high_cond_num_matrix:public test_func {
private:
	MyHCJob22 theJob;

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
	F22_rotated_hybrid_composition_3_high_cond_num_matrix (int dimension, double bias);
	F22_rotated_hybrid_composition_3_high_cond_num_matrix (int dimension, double bias, string file_data, string file_m);

	// Function body
	double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME22 = "Rotated Hybrid Composition Function 3 with High Condition Number Matrix";
const string DEFAULT_FILE_DATA22 = "supportData/hybrid_func3_data.txt";
const string DEFAULT_FILE_MX_PREFIX22 = "supportData/hybrid_func3_HM_D";
const string DEFAULT_FILE_MX_SUFFIX22 = ".txt";
static double m_sigma22[] = {
	1.0,    1.0,    1.0,    1.0,    1.0,
	2.0,    2.0,    2.0,    2.0,    2.0
};
static double m_lambda22[] = {
	5.0*5.0/100.0,  5.0/100.0,      5.0*1.0,        1.0,                    5.0*1.0,
	1.0,                    5.0*10.0,       10.0,           5.0*5.0/200.0,  5.0/200.0
};
static double m_func_biases22[] = {
	0.0,    100.0,  200.0,  300.0,  400.0,
	500.0,  600.0,  700.0,  800.0,  900.0
};
#endif
