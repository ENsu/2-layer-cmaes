#ifndef _F21_rotated_hybrid_composition_3
#define _F21_rotated_hybrid_composition_3
#include"test_func.h"
#include"HCJob.h"
#include"benchmark.h"
class MyHCJob21:public HCJob {
	public:
		double basic_func(int func_no, double *x,int length);
};
class F21_rotated_hybrid_composition_3:public test_func {
private:
	MyHCJob21 theJob;

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
	F21_rotated_hybrid_composition_3 (int dimension, double bias);
	F21_rotated_hybrid_composition_3 (int dimension, double bias, string file_data, string file_m);

	// Function body
	double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME21 = "Rotated Hybrid Composition Function 3";
const string DEFAULT_FILE_DATA21 = "supportData/hybrid_func3_data.txt";
const string DEFAULT_FILE_MX_PREFIX21 = "supportData/hybrid_func3_M_D";
const string DEFAULT_FILE_MX_SUFFIX21 = ".txt";
static double m_sigma21[] = {
	1.0,    1.0,    1.0,    1.0,    1.0,
	2.0,    2.0,    2.0,    2.0,    2.0
};
static double m_lambda21[] = {
	5.0*5.0/100.0,  5.0/100.0,      5.0*1.0,        1.0,                    5.0*1.0,
	1.0,                    5.0*10.0,       10.0,           5.0*5.0/200.0,  5.0/200.0
};
static double m_func_biases21[] = {
	0.0,    100.0,  200.0,  300.0,  400.0,
	500.0,  600.0,  700.0,  800.0,  900.0
};
#endif
