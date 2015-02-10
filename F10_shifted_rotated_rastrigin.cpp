#include"F10_shifted_rotated_rastrigin.h"
#include"benchmark.h"
// Constructors
F10_shifted_rotated_rastrigin::F10_shifted_rotated_rastrigin (int dimension, double bias) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME10;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_matrix = new double*[m_dimension];
	for(int i=0;i<m_dimension;++i)
		m_matrix[i] = new double[m_dimension];

	m_z = new double[m_dimension];
	m_zM = new double[m_dimension];

	// Load the shifted global optimum
	loadRowVectorFromFile(DEFAULT_FILE_DATA10, m_dimension, m_o);
	// Load the matrix
	char tmp[8];
	sprintf(tmp,"%d",dimension);
	loadMatrixFromFile(DEFAULT_FILE_MX_PREFIX10 + tmp + DEFAULT_FILE_MX_SUFFIX10, m_dimension, m_dimension, m_matrix);
}
F10_shifted_rotated_rastrigin::F10_shifted_rotated_rastrigin (int dimension, double bias, string file_data, string file_m) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME10;

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
}

// Function body
double F10_shifted_rotated_rastrigin::f(double *x,int length) {

	double result = 0.0;

	shift(m_z, x, m_o,length);
	rotate(m_zM, m_z, m_matrix,length);

	result = rastrigin(m_zM,length);

	result += m_bias;

	return (result);
}
