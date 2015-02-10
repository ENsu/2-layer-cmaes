#ifndef _F24_rotated_hybrid_composition_4
#define _F24_rotated_hybrid_composition_4
#include"test_func.h"
#include"HCJob.h"
#include"benchmark.h"
class MyHCJob24:public HCJob {
public:
	double basic_func(int func_no, double *x,int length);
};
class F24_rotated_hybrid_composition_4:public test_func {
private:
	MyHCJob24 theJob;

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
	F24_rotated_hybrid_composition_4 (int dimension, double bias);
	F24_rotated_hybrid_composition_4 (int dimension, double bias, string file_data, string file_m);

	// Function body
	double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME24 = "Rotated Hybrid Composition Function 4";
const string DEFAULT_FILE_DATA24 = "supportData/hybrid_func4_data.txt";
const string DEFAULT_FILE_MX_PREFIX24 = "supportData/hybrid_func4_M_D";
const string DEFAULT_FILE_MX_SUFFIX24 = ".txt";
static double m_sigma24[] = {
	2.0,    2.0,    2.0,    2.0,    2.0,
	2.0,    2.0,    2.0,    2.0,    2.0
};
static double m_lambda24[] = {
	10.0,           5.0/20.0,       1.0,    5.0/32.0,       1.0,
	5.0/100.0,      5.0/50.0,       1.0,    5.0/100.0,      5.0/100.0
};
static double m_func_biases24[] = {
	0.0,    100.0,  200.0,  300.0,  400.0,
	500.0,  600.0,  700.0,  800.0,  900.0
};
#endif
