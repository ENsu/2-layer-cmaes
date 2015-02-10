#include"F09_shifted_rastrigin.h"
#include"benchmark.h"

// Constructors
F09_shifted_rastrigin::F09_shifted_rastrigin (int dimension, double bias) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME9;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_z = new double[m_dimension];

	// Load the shifted global optimum
	loadRowVectorFromFile(DEFAULT_FILE_DATA9, m_dimension, m_o);
}
F09_shifted_rastrigin::F09_shifted_rastrigin (int dimension, double bias, string file_data) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME9;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_z = new double[m_dimension];

	// Load the shifted global optimum
	loadRowVectorFromFile(file_data, m_dimension, m_o);
}

// Function body
double F09_shifted_rastrigin::f(double *x,int length) {

	double result = 0.0;

	shift(m_z, x, m_o,length);

	result = rastrigin(m_z,length);

	result += m_bias;

	return (result);
}
