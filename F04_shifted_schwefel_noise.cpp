#include"F04_shifted_schwefel_noise.h"
#include"benchmark.h"

extern randomG random_;

// Constructors
F04_shifted_schwefel_noise::F04_shifted_schwefel_noise (int dimension, double bias) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME4;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_z = new double[m_dimension];

	// Load the shifted global optimum
	loadRowVectorFromFile(DEFAULT_FILE_DATA4, m_dimension, m_o);
}
F04_shifted_schwefel_noise::F04_shifted_schwefel_noise (int dimension, double bias, string file_data) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME4;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_z = new double[m_dimension];

	// Load the shifted global optimum
	loadRowVectorFromFile(file_data, m_dimension, m_o);
}

// Function body
double F04_shifted_schwefel_noise::f(double *x,int length) {
	double result = 0.0;

	shift(m_z, x, m_o,length);

	result = schwefel_102(m_z,length);

	// NOISE
	// Comment the next line to remove the noise
#ifndef no_noise
	result *= (1.0 + 0.4 * fabs(random_.normal01()));
#endif
	
	result += m_bias;

	return (result);
}
