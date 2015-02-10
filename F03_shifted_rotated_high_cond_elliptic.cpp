#include"F03_shifted_rotated_high_cond_elliptic.h"
#include"benchmark.h"

// Constructors
F03_shifted_rotated_high_cond_elliptic::F03_shifted_rotated_high_cond_elliptic (int dimension, double bias) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME3;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_matrix = new double*[m_dimension];
	for(int i=0;i<m_dimension;++i)
		m_matrix[i]=new double[m_dimension];

	m_z = new double[m_dimension];
	m_zM = new double[m_dimension];
	// Load the shifted global optimum
	loadRowVectorFromFile(DEFAULT_FILE_DATA3, m_dimension, m_o);
	// Load the matrix
	char tmp[8];
	sprintf(tmp,"%d",dimension);
	loadMatrixFromFile( DEFAULT_FILE_MX_PREFIX3 + tmp + DEFAULT_FILE_MX_SUFFIX3, m_dimension, m_dimension, m_matrix);
	constant = pow(1.0e6, 1.0/(m_dimension-1.0));
}
F03_shifted_rotated_high_cond_elliptic::F03_shifted_rotated_high_cond_elliptic (int dimension, double bias, string file_data, string file_m) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME3;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_matrix = new double*[m_dimension];
	for(int i=0;i<m_dimension;++i)
		m_matrix[i]=new double[m_dimension];

	m_z = new double[m_dimension];
	m_zM = new double[m_dimension];

	// Load the shifted global optimum
	loadRowVectorFromFile(file_data, m_dimension, m_o);
	// Load the matrix
	loadMatrixFromFile(file_m, m_dimension, m_dimension, m_matrix);

	constant = pow(1.0e6, 1.0/(m_dimension-1.0));
}

// Function body
double F03_shifted_rotated_high_cond_elliptic::f(double *x,int length) {
	double result = 0.0;

	shift(m_z, x, m_o,length);
	rotate(m_zM, m_z, m_matrix,length);
	
	double sum = 0.0;

	for (int i = 0 ; i < m_dimension ; i ++) {
		sum += pow(constant, (double)i) * m_zM[i] * m_zM[i];
	}

	result = sum + m_bias;

	return (result);
}
