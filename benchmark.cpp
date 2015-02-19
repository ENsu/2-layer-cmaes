#include "benchmark.h"
#include "F01_shifted_sphere.h"
#include "F02_shifted_schwefel.h"
#include "F03_shifted_rotated_high_cond_elliptic.h"
#include "F04_shifted_schwefel_noise.h"
#include "F05_schwefel_global_opt_bound.h"
#include "F06_shifted_rosenbrock.h"
#include "F07_shifted_rotated_griewank.h"
#include "F08_shifted_rotated_ackley_global_opt_bound.h"
#include "F09_shifted_rastrigin.h"
#include "F10_shifted_rotated_rastrigin.h"
#include "F11_shifted_rotated_weierstrass.h"
#include "F12_schwefel.h"
#include "F13_shifted_expanded_griewank_rosenbrock.h"
#include "F14_shifted_rotated_expanded_scaffer.h"
#include "F15_hybrid_composition_1.h"
#include "F16_rotated_hybrid_composition_1.h"
#include "F17_rotated_hybrid_composition_1_noise.h"
#include "F18_rotated_hybrid_composition_2.h"
#include "F19_rotated_hybrid_composition_2_narrow_basin_global_opt.h"
#include "F20_rotated_hybrid_composition_2_global_opt_bound.h"
#include "F21_rotated_hybrid_composition_3.h"
#include "F22_rotated_hybrid_composition_3_high_cond_num_matrix.h"
#include "F23_noncontinuous_rotated_hybrid_composition_3.h"
#include "F24_rotated_hybrid_composition_4.h"
#include "F25_rotated_hybrid_composition_4_bound.h"
randomG random_;

// Constructors
//	Load the data common to all test functions.
//	Get ready for creating test function instances.


void initial() {
	m_biases = new double[NUM_TEST_FUNC];
	m_iSqrt = new double[MAX_SUPPORT_DIM];
	random_.randomize(0);
	loadRowVectorFromFile(DEFAULT_FILE_BIAS, NUM_TEST_FUNC, m_biases);
	for (int i = 0 ; i < MAX_SUPPORT_DIM ; i ++) {
		m_iSqrt[i] = sqrt(((double )i) + 1.0);
	}
}


//
// Use this function to manufacture new test function objects
//
test_func* testFunctionFactory(int func_num, int dimension) {
	test_func *returnFunc=NULL;
	switch(func_num){
		case 1:
			returnFunc=new F01_shifted_sphere(dimension,m_biases[func_num-1]);
			break;
		case 2:
			returnFunc=new F02_shifted_schwefel(dimension,m_biases[func_num-1]);
			break;
		case 3:
			returnFunc=new F03_shifted_rotated_high_cond_elliptic(dimension,m_biases[func_num-1]);
			break;
		case 4:
			returnFunc=new F04_shifted_schwefel_noise(dimension,m_biases[func_num-1]);
			break;
		case 5:
			returnFunc=new F05_schwefel_global_opt_bound(dimension,m_biases[func_num-1]);
			break;
		case 6:
			returnFunc=new F06_shifted_rosenbrock(dimension,m_biases[func_num-1]);
			break;
		case 7:
			returnFunc=new F07_shifted_rotated_griewank(dimension,m_biases[func_num-1]);
			break;
		case 8:
			returnFunc=new F08_shifted_rotated_ackley_global_opt_bound(dimension,m_biases[func_num-1]);
			break;
		case 9:
			returnFunc=new F09_shifted_rastrigin(dimension,m_biases[func_num-1]);
			break;
		case 10:
			returnFunc=new F10_shifted_rotated_rastrigin(dimension,m_biases[func_num-1]);
			break;
		case 11:
			returnFunc=new F11_shifted_rotated_weierstrass(dimension,m_biases[func_num-1]);
			break;
		case 12:
			returnFunc=new F12_schwefel(dimension,m_biases[func_num-1]);
			break;
		case 13:
			returnFunc=new F13_shifted_expanded_griewank_rosenbrock(dimension,m_biases[func_num-1]);
			break;
		case 14:
			returnFunc=new F14_shifted_rotated_expanded_scaffer(dimension,m_biases[func_num-1]);
			break;
		case 15:
			returnFunc=new F15_hybrid_composition_1(dimension,m_biases[func_num-1]);
			break;
		case 16:
			returnFunc=new F16_rotated_hybrid_composition_1(dimension,m_biases[func_num-1]);
			break;
		case 17:
			returnFunc=new F17_rotated_hybrid_composition_1_noise(dimension,m_biases[func_num-1]);
			break;
		case 18:
			returnFunc=new F18_rotated_hybrid_composition_2(dimension,m_biases[func_num-1]);
			break;
		case 19:
			returnFunc=new F19_rotated_hybrid_composition_2_narrow_basin_global_opt(dimension,m_biases[func_num-1]);
			break;
		case 20:
			returnFunc=new F20_rotated_hybrid_composition_2_global_opt_bound(dimension,m_biases[func_num-1]);
			break;
		case 21:
			returnFunc=new F21_rotated_hybrid_composition_3(dimension,m_biases[func_num-1]);
			break;
		case 22:
			returnFunc=new F22_rotated_hybrid_composition_3_high_cond_num_matrix(dimension,m_biases[func_num-1]);
			break;
		case 23:
			returnFunc=new F23_noncontinuous_rotated_hybrid_composition_3(dimension,m_biases[func_num-1]);
			break;
		case 24:
			returnFunc=new F24_rotated_hybrid_composition_4(dimension,m_biases[func_num-1]);
			break;
		case 25:
			returnFunc=new F25_rotated_hybrid_composition_4_bound(dimension,m_biases[func_num-1]);
			break;
		default:
			returnFunc=new F01_shifted_sphere(dimension,m_biases[func_num-1]);
			break;
	}
	return returnFunc;
}

// Run tests on the test functions
//	< 0:	Error!
//	= 0:	On all the functions
//	> 0:	On the specified test function

void runTest() {
	runTest(0);
}

void runTest(int func_num) {
	if (func_num == 0) {
		for (int i = 1 ; i <= NUM_TEST_FUNC ; i ++)
			runTest(i);
	}
	else if ((func_num < 0) || (func_num > NUM_TEST_FUNC)) {
		puts("The specified func_num is out of range.");
		exit(-1);
	}
	else {
		// Run the test function against the check points
		int num_test_points = 10;
		int test_dimension = 50;

		char tmp[128];

		double *test_f = new double[num_test_points];

		double **test_x = new double *[num_test_points];
		for(int i=0;i<num_test_points;++i)
			test_x[i]=new double[test_dimension];

		sprintf(tmp,"testData/test_data_func%d.txt",func_num);

		string file_test(tmp);

		// Create the test function object
		test_func *aFunc = testFunctionFactory(func_num, test_dimension);

		cout << "Run tests on function " << func_num << " (" << aFunc->name() << "):" << endl;
		cout << "  " << num_test_points << " " << aFunc->dimension() << "-dimension check points" << endl;

		loadTestDataFromFile(file_test, num_test_points, test_dimension, test_x, test_f);

		for (int i = 0 ; i < num_test_points ; i ++) {
			// Execute the test function
			// Collect and compare the results
			double result = aFunc->f(test_x[i],50);
			double diff = result - test_f[i];
			double ratio = fabs(diff / test_f[i]);
			printf("    %.16E - %.16E = %.16E\n",result,test_f[i],diff);
			printf("        Difference ratio = %.16E",ratio);
			if (ratio != 0.0) {
				if (ratio <= 1e-12) {
					cout << " (<= 1E-12)" << endl;
				}
				else {
					cout << " (> 1E-12) *****" << endl;
				}
			}
			else {
				cout << endl;
			}
		}

		delete [] test_f;
		for(int i=0;i<num_test_points;++i)
			delete [] test_x[i];
	}
}

//
// Basic functions
//


// Sphere function
double sphere(double *x,int length) {

	double sum = 0.0;

	for (int i = 0 ; i < length ; i ++) {
	    sum += x[i] * x[i];
	}
	return (sum);
}

// Sphere function with noise
double sphere_noise(double *x,int length) {

	double sum = 0.0;

	for (int i = 0 ; i < length ; i ++) {
		sum += x[i] * x[i];
	}

	// NOISE
	// Comment the next line to remove the noise
#ifndef no_noise
	sum *= (1.0 + 0.1 * fabs(random_.normal01()));
#endif

	return (sum);
}

// Schwefel's problem 1.2
double schwefel_102(double *x,int length) {

	double prev_sum, curr_sum, outer_sum;

	curr_sum = x[0];
	outer_sum = (curr_sum * curr_sum);

	for (int i = 1 ; i < length ; i ++) {
		prev_sum = curr_sum;
		curr_sum = prev_sum + x[i];
		outer_sum += (curr_sum * curr_sum);
	}

	return (outer_sum);
}

// Rosenbrock's function
double rosenbrock(double *x,int length) {

	double sum = 0.0;

	for (int i = 0 ; i < (length-1) ; i ++) {
		double temp1 = (x[i] * x[i]) - x[i+1];
		double temp2 = x[i] - 1.0;
		sum += (100.0 * temp1 * temp1) + (temp2 * temp2);
	}

	return (sum);
}

// F2: Rosenbrock's Function -- 2D version
double F2(double x, double y) {
	double temp1 = (x * x) - y;
	double temp2 = x - 1.0;
	return ((100.0 * temp1 * temp1) + (temp2 * temp2));
}

// Griewank's function
double griewank(double *x,int length) {

	double sum = 0.0;
	double product = 1.0;

	for (int i = 0 ; i < length ; i ++) {
		sum += ((x[i] * x[i]) / 4000.0);
		product *= cos( x[i] / m_iSqrt[i] );
	}

	return (sum - product + 1.0);
}

// F8: Griewank's Function -- 1D version
double F8(double x) {
	return (((x * x) / 4000.0) - cos(x) + 1.0);
}

// Ackley's function
double ackley(double *x,int length) {

	double sum1 = 0.0;
	double sum2 = 0.0;

	for (int i = 0 ; i < length ; i ++) {
		sum1 += (x[i] * x[i]);
		sum2 += (cos(PIx2 * x[i]));
	}

	return (-20.0 * exp(-0.2 * sqrt(sum1 / ((double )length))) - exp(sum2 / ((double )length)) + 20.0 + M_E);
}

bool signbit(double x){
	if(x < 0)
		return true;
	else
		return false;
}

double round(double x){		//only for positive value x 
	if((x - (int)x) < 0.5)
		return (int)x;
	else
		return (int)x+1;
}

// Round function
// 0. Use the Matlab version for rounding numbers
double myRound(double x) {
	return ( (signbit(x)==0?1.0:-1.0) * round(fabs(x)));
}
// 1. "o" is provided
double myXRound(double x, double o) {
	return ((fabs(x - o) < 0.5) ? x : (myRound(2.0 * x) / 2.0));
}
// 2. "o" is not provided
double myXRound(double x) {
	return ((fabs(x) < 0.5) ? x : (myRound(2.0 * x) / 2.0));
}

// Rastrigin's function
double rastrigin(double *x,int length) {

	double sum = 0.0;

	for (int i = 0 ; i < length ; i ++) {
		sum += (x[i] * x[i]) - (10.0 * cos(PIx2 * x[i])) + 10.0;
	}

	return (sum);
}

// Non-Continuous Rastrigin's function
double rastriginNonCont(double *x,int length) {

	double sum = 0.0;
	double currX;

	for (int i = 0 ; i < length ; i ++) {
		currX = myXRound(x[i]);
		sum += (currX * currX) - (10.0 * cos(PIx2 * currX)) + 10.0;
	}

	return (sum);
}


// Weierstrass function
double weierstrass(double *x,int length) {
	return (weierstrass(x, 0.5, 3.0, 20,length));
}

double weierstrass(double *x, double a, double b, int Kmax,int length) {

	double sum1 = 0.0;
	for (int i = 0 ; i < length ; i ++) {
		for (int k = 0 ; k <= Kmax ; k ++) {
			sum1 += pow(a, k) * cos(PIx2 * pow(b, (double)k) * (x[i] + 0.5));
		}
	}

	double sum2 = 0.0;
	for (int k = 0 ; k <= Kmax ; k ++) {
		sum2 += pow(a, (double)k) * cos(PIx2 * pow(b, (double)k) * (0.5));
	}

	return (sum1 - sum2*((double )(length)));
}

// F8F2
double F8F2(double *x,int length) {

	double sum = 0.0;

	for (int i = 1 ; i < length ; i ++) {
		sum += F8(F2(x[i-1], x[i]));
	}
	sum += F8(F2(x[length-1], x[0]));

	return (sum);
}

// Scaffer's F6 function
double ScafferF6(double x, double y) {
	double temp1 = x*x + y*y;
	double temp2 = sin(sqrt(temp1));
	double temp3 = 1.0 + 0.001 * temp1;
	return (0.5 + ((temp2 * temp2 - 0.5)/(temp3 * temp3)));
}

double EScafferF6(double *x,int length) {

	double sum = 0.0;

	for (int i = 1 ; i < length ; i ++) {
		sum += ScafferF6(x[i-1], x[i]);
	}
	sum += ScafferF6(x[length-1], x[0]);

	return (sum);
}

// Non-Continuous Expanded Scaffer's F6 function
double EScafferF6NonCont(double *x,int length) {

	double sum = 0.0;
	double prevX, currX;

	currX = myXRound(x[0]);
	for (int i = 1 ; i < length ; i ++) {
		prevX = currX;
		currX = myXRound(x[i]);
		sum += ScafferF6(prevX, currX);
	}
	prevX = currX;
	currX = myXRound(x[0]);
	sum += ScafferF6(prevX, currX);

	return (sum);
}

// Elliptic
double elliptic(double *x,int length) {

	double sum = 0.0;
	double a = 1e6;

	for (int i = 0 ; i < length ; i ++) {
		sum += pow(a, (((double )i)/((double )(length-1)))) * x[i] * x[i];
	}

	return (sum);
}

// Hybrid composition
double hybrid_composition(double *x, HCJob &job,int length) {

	int num_func = job.num_func;
	int num_dim = job.num_dim;

	// Get the raw weights
	double wMax = -HUGE_VAL ; // negative inf
	for (int i = 0 ; i < num_func ; i ++) {
		double sumSqr = 0.0;
		shift(job.z[i], x, job.o[i],length);
		for (int j = 0 ; j < num_dim ; j ++) {
			sumSqr += (job.z[i][j] * job.z[i][j]);
		}
		job.w[i] = exp(-1.0 * sumSqr / (2.0 * num_dim * job.sigma[i] * job.sigma[i]));
		if (wMax < job.w[i])
			wMax = job.w[i];
	}

	// Modify the weights
	double wSum = 0.0;
	double w1mMaxPow = 1.0 - pow(wMax, 10.0);
	for (int i = 0 ; i < num_func ; i ++) {
		if (job.w[i] != wMax) {
			job.w[i] *= w1mMaxPow;
		}
		wSum += job.w[i];
	}

	// Normalize the weights
	for (int i = 0 ; i < num_func ; i ++) {
		job.w[i] /= wSum;
	}

	double sumF = 0.0;
	for (int i = 0 ; i < num_func ; i ++) {
		for (int j = 0 ; j < num_dim ; j ++) {
			job.z[i][j] /= job.lambda[i];
		}
		rotate(job.zM[i], job.z[i], job.M[i],length);
		sumF +=
			job.w[i] *
			(
			 job.C * job.basic_func(i, job.zM[i],length) / job.fmax[i] +
			 job.biases[i]
			);
	}
	return (sumF);
}

//
// Elementary operations
//

// Shift
void shift(double *results, double *x, double *o,int length) {
	for (int i = 0 ; i < length ; i ++) {
	    results[i] = x[i] - o[i];
	}
}

// Rotate
void rotate(double *results, double *x, double **matrix,int length) {
	xA(results, x, matrix,length);
}


//
// Matrix & vector operations
//

// (1xD) row vector * (Dx1) column vector = (1) scalar
double xy(double *x, double *y,int length) {
	double result = 0.0;
	for (int i = 0 ; i < length ; i ++) {
		result += (x[i] * y[i]);
	}

	return (result);
}

// (1xD) row vector * (DxD) matrix = (1xD) row vector
void xA(double *result, double *x, double **A,int length) {
	for (int i = 0 ; i < length ; i ++) {
		result[i] = 0.0;
		for (int j = 0 ; j < length ; j ++) {
			result[i] += (x[j] * A[j][i]);
		}
	}
}

// (DxD) matrix * (Dx1) column vector = (Dx1) column vector
void Ax(double *result, double **A, double *x,int length) {
	for (int i = 0 ; i < length ; i ++) {
		result[i] = 0.0;
		for (int j = 0 ; j < length ; j ++) {
			result[i] += (A[i][j] * x[j]);
		}
	}
}

//
// Utility functions for loading data from the given text file
//
void loadTestDataFromFile(string file, int num_test_points, int test_dimension, double **x, double *f) {
	ifstream brSrc(file.c_str());
	loadMatrix(brSrc, num_test_points, test_dimension, x);
	loadColumnVector(brSrc, num_test_points, f);
	brSrc.close();
}

void loadRowVectorFromFile(string file, int columns, double *row) { 
	ifstream brSrc(file.c_str());
	loadRowVector(brSrc, columns, row);
	brSrc.close();
}

void loadRowVector(ifstream &brSrc, int columns, double *row){
	string stToken;
	char *p,*str;
	getline(brSrc , stToken);
	str=new char[stToken.length()+1];
	strcpy(str,stToken.c_str());
	p=strtok(str," ");
	for (int i = 0 ; i < columns ; i++) {
		row[i] = atof(p);		//transfer from char to double
		p=strtok(NULL," ");
	}
	delete [] str;
}

void loadColumnVectorFromFile(string file, int rows, double *column) {
	ifstream brSrc(file.c_str());
	loadColumnVector(brSrc, rows, column);
	brSrc.close();
}

void loadColumnVector(ifstream &brSrc, int rows, double *column) {
	string stToken;
	for (int i = 0 ; i < rows ; i ++) {
		brSrc >> stToken;
		column[i] = atof(stToken.c_str());
	}
}

void loadNMatrixFromFile(string file, int N, int rows, int columns, double ***matrix) {
	ifstream brSrc(file.c_str());
	for (int i = 0 ; i < N ; i ++) {
		loadMatrix(brSrc, rows, columns, matrix[i]);
	}
	brSrc.close();
}

void loadMatrixFromFile(string file, int rows, int columns, double **matrix) {
	ifstream brSrc(file.c_str());
	loadMatrix(brSrc, rows, columns, matrix);
	brSrc.close();
}

void loadMatrix(ifstream &brSrc, int rows, int columns, double **matrix) {
	for (int i = 0 ; i < rows ; i ++) {
		loadRowVector(brSrc, columns, matrix[i]);
	}
}

double interval_f(double x, double y, double z)  ///for building required problems more easily 
{
	if(x>=y && x<=z)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
