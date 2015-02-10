#ifndef _benchmark
#define _benchmark
#include<cmath>
#include<string>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<iostream>
#include<iomanip>
#include<fstream>
#include"random.hpp"
#include"test_func.h"
#include"HCJob.h"

#include<gsl/gsl_math.h>

//#define no_noise

	// Class variables
	static double *m_iSqrt;

	// Instance variables
	static double *m_biases;
	
	// Fixed (class) parameters
	static const int NUM_TEST_FUNC = 25;
//	static const string DEFAULT_FILE_BIAS;
//	static const string test_func_class_names[25];

	// For certain functions, some essential data can be calculated beforehand.
	// Hence, a maximum supported number of dimensions should be specified.
	// Specifiy the number of dimensions here if you need more.
	static const int MAX_SUPPORT_DIM = 100;
	static const double PIx2 = M_PI * 2.0;

	// Formatter for the number representation
//	static DecimalFormat scientificFormatter = new DecimalFormat("0.0000000000000000E00");
//	static DecimalFormat numberFormatter = scientificFormatter;
//	static DecimalFormat percentageFormatter = new DecimalFormat("0.0000000000");

	// Random number generator
	// If you want to plan a specific seed, do it here.

	// Class loader & reflection
//	static final ClassLoader loader = ClassLoader.getSystemClassLoader();
//	static final Class[] test_func_arg_types = { int.class, double.class };



	void initial();
	
	//
	// Use this function to manufacture new test function objects
	//
extern	test_func* testFunctionFactory(int func_num, int dimension);

	// Run tests on the test functions
	//	< 0:	Error!
	//	= 0:	On all the functions
	//	> 0:	On the specified test function
	void runTest(int func_num);
	void runTest();


	//
	// Basic functions
	//

	// Sphere function
	double sphere(double *x,int length);

	// Sphere function with noise
	double sphere_noise(double *x,int length);

	// Schwefel's problem 1.2
	double schwefel_102(double *x,int length);

	// Rosenbrock's function
	double rosenbrock(double *x,int length);

	// F2: Rosenbrock's Function -- 2D version
	double F2(double x, double y);

	// Griewank's function
	double griewank(double *x,int length);

	// F8: Griewank's Function -- 1D version
	double F8(double x);
	
	// Ackley's function
	double ackley(double *x,int length);
	
	// Round function
	// 0. Use the Matlab version for rounding numbers
	double myRound(double x);
	// 1. "o" is provided
	double myXRound(double x, double o);
	// 2. "o" is not provided
	double myXRound(double x);

	// Rastrigin's function
	double rastrigin(double *x,int length);

	// Non-Continuous Rastrigin's function
	double rastriginNonCont(double *x,int length);

	// Weierstrass function
	double weierstrass(double *x,int length);

	double weierstrass(double *x, double a, double b, int Kmax,int length);

	// F8F2
	double F8F2(double *x,int length);

	// Scaffer's F6 function
	double ScafferF6(double x, double y);

	// Expanded Scaffer's F6 function
	double EScafferF6(double *x,int length);

	// Non-Continuous Expanded Scaffer's F6 function
	double EScafferF6NonCont(double *x,int length);

	// Elliptic
	double elliptic(double *x,int length);

	// Hybrid composition
	double hybrid_composition(double *x, HCJob &job,int length);
	
	//
	// Elementary operations
	//

	// Shift
	void shift(double *results, double *x, double *o,int length);

	// Rotate
	void rotate(double *results, double *x, double **matrix,int length);

	//
	// Matrix & vector operations
	//

	// (1xD) row vector * (Dx1) column vector = (1) scalar
	double xy(double *x, double *y,int length);

	// (1xD) row vector * (DxD) matrix = (1xD) row vector
	void xA(double *result, double *x, double **A,int length);

	// (DxD) matrix * (Dx1) column vector = (Dx1) column vector
	void Ax(double *result, double **A, double *x,int length);

	//
	// Utility functions for loading data from the given text file
	//
	void loadTestDataFromFile(string file, int num_test_points, int test_dimension, double **x, double *f);

	void loadRowVectorFromFile(string file, int columns, double *row);

	void loadRowVector(ifstream &brSrc, int columns, double *row);

	void loadColumnVectorFromFile(string file, int rows, double *column);
	
	void loadColumnVector(ifstream &brSrc, int rows, double *column);

	void loadNMatrixFromFile(string file, int N, int rows, int columns, double ***matrix);

	void loadMatrixFromFile(string file, int rows, int columns, double **matrix);

	void loadMatrix(ifstream &brSrc, int rows, int columns, double **matrix);

	double interval_f(double x, double y, double z);
	

const string DEFAULT_FILE_BIAS = "supportData/fbias_data.txt";
const string test_func_class_names[25] = {
		"F01_shifted_sphere",
		"F02_shifted_schwefel",
		"F03_shifted_rotated_high_cond_elliptic",
		"F04_shifted_schwefel_noise",
		"F05_schwefel_global_opt_bound",
		"F06_shifted_rosenbrock",
		"F07_shifted_rotated_griewank",
		"F08_shifted_rotated_ackley_global_opt_bound",
		"F09_shifted_rastrigin",
		"F10_shifted_rotated_rastrigin",
		"F11_shifted_rotated_weierstrass",
		"F12_schwefel",
		"F13_shifted_expanded_griewank_rosenbrock",
		"F14_shifted_rotated_expanded_scaffer",
		"F15_hybrid_composition_1",
		"F16_rotated_hybrid_composition_1",
		"F17_rotated_hybrid_composition_1_noise",
		"F18_rotated_hybrid_composition_2",
		"F19_rotated_hybrid_composition_2_narrow_basin_global_opt",
		"F20_rotated_hybrid_composition_2_global_opt_bound",
		"F21_rotated_hybrid_composition_3",
		"F22_rotated_hybrid_composition_3_high_cond_num_matrix",
		"F23_noncontinuous_rotated_hybrid_composition_3",
		"F24_rotated_hybrid_composition_4",
		"F25_rotated_hybrid_composition_4_bound"
	};

#endif
