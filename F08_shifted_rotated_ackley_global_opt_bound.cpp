#include"F08_shifted_rotated_ackley_global_opt_bound.h"
#include"benchmark.h"

// Constructors
F08_shifted_rotated_ackley_global_opt_bound::F08_shifted_rotated_ackley_global_opt_bound (int dimension, double bias) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME8;
	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_matrix = new double*[m_dimension];
	for(int i=0;i<m_dimension;++i)
		m_matrix[i] = new double[m_dimension];

	m_z = new double[m_dimension];
	m_zM = new double[m_dimension];

	// Load the shifted global optimum
	loadRowVectorFromFile(DEFAULT_FILE_DATA8, m_dimension, m_o);
	// Load the matrix
	char tmp[8];
	sprintf(tmp,"%d",dimension);
	loadMatrixFromFile(DEFAULT_FILE_MX_PREFIX8 + tmp + DEFAULT_FILE_MX_SUFFIX8, m_dimension, m_dimension, m_matrix);

	for (int i = 0 ; i < m_dimension ; i += 2) {
		m_o[i] = -32.0;
	}
}
F08_shifted_rotated_ackley_global_opt_bound::F08_shifted_rotated_ackley_global_opt_bound (int dimension, double bias, string file_data, string file_m) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME8;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_matrix = new double*[m_dimension];
	for(int i=0;i<m_dimension;++i)
		m_matrix[i] = new double[m_dimension];

	m_z = new double[m_dimension];
	m_zM = new double[m_dimension];

	// Load the shifted global optimum
	loadRowVectorFromFile(file_data, m_dimension, m_o);
	// Load the matrix
	loadMatrixFromFile(file_m, m_dimension, m_dimension, m_matrix);

	for (int i = 0 ; i < m_dimension ; i += 2) {
		m_o[i] = -32.0;
	}
}

// Function body
double F08_shifted_rotated_ackley_global_opt_bound::f(double *x,int length) {

	double result = 0.0;

	shift(m_z, x, m_o,length);
	rotate(m_zM, m_z, m_matrix,length);

	result = ackley(m_zM,length);

	result += m_bias;

	return (result);
}
