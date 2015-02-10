#ifndef _F03_shifted_rotated_high_cond_elliptic
#define _F03_shifted_rotated_high_cond_elliptic
#include"test_func.h"
class F03_shifted_rotated_high_cond_elliptic:public test_func {

private:
	// Shifted global optimum
	double *m_o;
	double **m_matrix;

	// In order to avoid excessive memory allocation,
	// a fixed memory buffer is allocated for each function object.
	double *m_z;
	double *m_zM;

	double constant;

public:
	// Constructors
	F03_shifted_rotated_high_cond_elliptic (int dimension, double bias);
	F03_shifted_rotated_high_cond_elliptic (int dimension, double bias, string file_data, string file_m);

	// Function body
	double f(double *x,int length);
};

// Fixed (class) parameters
const string FUNCTION_NAME3 = "Shifted Rotated High Conditioned Elliptic Function";
const string DEFAULT_FILE_DATA3 = "supportData/high_cond_elliptic_rot_data.txt";
const string DEFAULT_FILE_MX_PREFIX3 = "supportData/elliptic_M_D";
const string DEFAULT_FILE_MX_SUFFIX3 = ".txt";
#endif
