#include"F13_shifted_expanded_griewank_rosenbrock.h"
#include"benchmark.h"

// Constructors
F13_shifted_expanded_griewank_rosenbrock::F13_shifted_expanded_griewank_rosenbrock (int dimension, double bias) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME13;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_z = new double[m_dimension];

	// Load the shifted global optimum
	loadRowVectorFromFile(DEFAULT_FILE_DATA13, m_dimension, m_o);

	// z = x - o + 1 = x - (o - 1)
	// Do the "(o - 1)" part first
	for (int i = 0 ; i < m_dimension ; i ++) {
		m_o[i] -= 1.0;
	}
}
F13_shifted_expanded_griewank_rosenbrock::F13_shifted_expanded_griewank_rosenbrock (int dimension, double bias, string file_data) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME13;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_z = new double[m_dimension];

	// Load the shifted global optimum
	loadRowVectorFromFile(file_data, m_dimension, m_o);

	// z = x - o + 1 = x - (o - 1)
	// Do the "(o - 1)" part first
	for (int i = 0 ; i < m_dimension ; i ++) {
		m_o[i] -= 1.0;
	}
}

// Function body
double F13_shifted_expanded_griewank_rosenbrock::f(double *x,int length) {

	double result = 0.0;

	shift(m_z, x, m_o,length);

	result = F8F2(m_z,length);

	result += m_bias;

	return (result);
}
