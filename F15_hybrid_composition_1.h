#ifndef _F15_hybrid_composition_1
#define _F15_hybrid_composition_1
#include"test_func.h"
#include"HCJob.h"
#include"benchmark.h"
class MyHCJob15:public HCJob {
public:
	double basic_func(int func_no, double *x,int length);
};
class F15_hybrid_composition_1:public test_func {
private:
	// Number of functions
	const static int NUM_FUNC = 10;
	MyHCJob15 theJob;

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
	// Constructors
	F15_hybrid_composition_1 (int dimension, double bias);
	F15_hybrid_composition_1 (int dimension, double bias, string file_data);

	// Function body
	double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME15 = "Hybrid Composition Function 1";
const string DEFAULT_FILE_DATA15 = "supportData/hybrid_func1_data.txt";
static double m_sigma15[] = {
	1.0,    1.0,    1.0,    1.0,    1.0,
	1.0,    1.0,    1.0,    1.0,    1.0
};
static double m_lambda15[] = {
	1.0,            1.0,            10.0,           10.0,
	5.0/60.0,       5.0/60.0,       5.0/32.0,       5.0/32.0,
	5.0/100.0,      5.0/100.0
};
static double m_func_biases15[] = {
	0.0,    100.0,  200.0,  300.0,  400.0,
	500.0,  600.0,  700.0,  800.0,  900.0
};
#endif
