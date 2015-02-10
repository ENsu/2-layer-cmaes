#ifndef _F11_shifted_rotated_weierstrass
#define _F11_shifted_rotated_weierstrass
#include"test_func.h"
#include<math.h>
#include<gsl/gsl_math.h>
class F11_shifted_rotated_weierstrass:public test_func {
	private:
		// Shifted global optimum
		double *m_o;
		double **m_matrix;

		// In order to avoid excessive memory allocation,
		// a fixed memory buffer is allocated for each function object.
		double *m_z;
		double *m_zM;

	public:
		//static const double PIx2 = M_PI * 2.0;
		static const int Kmax = 20;
		//static const double a = 0.5;
		//static const double b = 3.0;

		// Constructors
		F11_shifted_rotated_weierstrass (int dimension, double bias);
		F11_shifted_rotated_weierstrass (int dimension, double bias, string file_data, string file_m);

		// Function body
		double f(double *x,int length);
};
// Fixed (class) parameters
const string FUNCTION_NAME11 = "Shifted Rotated Weierstrass Function";
const string DEFAULT_FILE_DATA11 = "supportData/weierstrass_data.txt";
const string DEFAULT_FILE_MX_PREFIX11 = "supportData/weierstrass_M_D";
const string DEFAULT_FILE_MX_SUFFIX11 = ".txt";
//static const double PIx2 = M_PI * 2.0;
//static const int Kmax = 20;
static const double a = 0.5;
static const double b = 3.0;
#endif
