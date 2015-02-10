#include"F06_shifted_rosenbrock.h"
#include"benchmark.h"
// Constructors
F06_shifted_rosenbrock::F06_shifted_rosenbrock (int dimension, double bias) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME6;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_z = new double[m_dimension];

	// Load the shifted global optimum
	loadRowVectorFromFile(DEFAULT_FILE_DATA6, m_dimension, m_o);

	// z = x - o + 1 = x - (o - 1)
	// Do the "(o - 1)" part first
	for (int i = 0 ; i < m_dimension ; i ++) {
		m_o[i] -= 1.0;
	}

}
F06_shifted_rosenbrock::F06_shifted_rosenbrock (int dimension, double bias, string file_data) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME6;

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
double F06_shifted_rosenbrock::f(double *x,int length) {

	double result = 0.0;

	shift(m_z, x, m_o,length);

	result = rosenbrock(m_z,length);

	result += m_bias;

	return (result);
}
