#ifndef _F19_rotated_hybrid_composition_2_narrow_basin_global_opt
#define _F19_rotated_hybrid_composition_2_narrow_basin_global_opt
#include"test_func.h"
#include"HCJob.h"
#include"benchmark.h"
class MyHCJob19:public HCJob {
	public:
		double basic_func(int func_no, double *x,int length);
};
class F19_rotated_hybrid_composition_2_narrow_basin_global_opt:public test_func {
private:
	MyHCJob19 theJob;

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
	F19_rotated_hybrid_composition_2_narrow_basin_global_opt (int dimension, double bias);
	F19_rotated_hybrid_composition_2_narrow_basin_global_opt (int dimension, double bias, string file_data, string file_m);

	// Function body
	double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME19 = "Rotated Hybrid Composition Function 2 with narrow basin global optimum";
const string DEFAULT_FILE_DATA19 = "supportData/hybrid_func2_data.txt";
const string DEFAULT_FILE_MX_PREFIX19 = "supportData/hybrid_func2_M_D";
const string DEFAULT_FILE_MX_SUFFIX19 = ".txt";
static double m_sigma19[] = {
	0.1,    2.0,    1.5,    1.5,    1.0,    1.0,
	1.5,    1.5,    2.0,    2.0
};
static double m_lambda19[] = {
	0.1*5.0/32.0,   5.0/32.0,       2.0*1,  1.0,            2.0*5.0/100.0,
	5.0/100.0,      2.0*10.0,       10.0,   2.0*5.0/60.0,   5.0/60.0
};
static double m_func_biases19[] = {
	0.0,    100.0,  200.0,  300.0,  400.0,
	500.0,  600.0,  700.0,  800.0,  900.0
};
#endif
