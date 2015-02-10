#include"F12_schwefel.h"
#include"benchmark.h"
// Constructors
F12_schwefel::F12_schwefel (int dimension, double bias) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME12;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_a = new double*[m_dimension];
	for(int i=0;i<m_dimension;++i)
		m_a[i] = new double[m_dimension];
	m_b = new double*[m_dimension];
	for(int i=0;i<m_dimension;++i)
		m_b[i] = new double[m_dimension];

	m_A = new double[m_dimension];
	m_B = new double[m_dimension];

	// Data:
	//      1. a            100x100
	//      2. b            100x100
	//      3. alpha        1x100
	double **m_data = new double*[100+100+1];
	for(int i=0;i<100+100+1;++i)
		m_data[i] = new double[m_dimension];

	// Load the shifted global optimum
	loadMatrixFromFile(DEFAULT_FILE_DATA12, 100+100+1, m_dimension, m_data);
	for (int i = 0 ; i < m_dimension ; i ++) {
		for (int j = 0 ; j < m_dimension ; j ++) {
			m_a[i][j] = m_data[i][j];
			m_b[i][j] = m_data[100+i][j];
		}
		m_o[i] = m_data[100+100][i];
	}

	for (int i = 0 ; i < m_dimension ; i ++) {
		m_A[i] = 0.0;
		for (int j = 0 ; j < m_dimension ; j ++) {
			m_A[i] += (m_a[i][j] * sin(m_o[j]) + m_b[i][j] * cos(m_o[j]));
		}
	}

}
F12_schwefel::F12_schwefel (int dimension, double bias, string file_data) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME12;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_a = new double*[m_dimension];
	for(int i=0;i<m_dimension;++i)
		m_a[i] = new double[m_dimension];
	m_b = new double*[m_dimension];
	for(int i=0;i<m_dimension;++i)
		m_b[i] = new double[m_dimension];

	m_A = new double[m_dimension];
	m_B = new double[m_dimension];

	// Data:
	//	1. a 		100x100
	//	2. b 		100x100
	//	3. alpha	1x100
	double **m_data = new double*[100+100+1];
	for(int i=0;i<100+100+1;++i)
		m_data[i] = new double[m_dimension];

	// Load the shifted global optimum
	loadMatrixFromFile(file_data, 100+100+1, m_dimension, m_data);
	for (int i = 0 ; i < m_dimension ; i ++) {
		for (int j = 0 ; j < m_dimension ; j ++) {
			m_a[i][j] = m_data[i][j];
			m_b[i][j] = m_data[100+i][j];
		}
		m_o[i] = m_data[100+100][i];
	}

	for (int i = 0 ; i < m_dimension ; i ++) {
		m_A[i] = 0.0;
		for (int j = 0 ; j < m_dimension ; j ++) {
			m_A[i] += (m_a[i][j] * sin(m_o[j]) + m_b[i][j] * cos(m_o[j]));
		}
	}
}

// Function body
double F12_schwefel::f(double *x,int length) {

	double sum = 0.0;

	for (int i = 0 ; i < m_dimension ; i ++) {
		m_B[i] = 0.0;
		for (int j = 0 ; j < m_dimension ; j ++) {
			m_B[i] += (m_a[i][j] * sin(x[j]) + m_b[i][j] * cos(x[j]));
		}

		double temp = m_A[i] - m_B[i];
		sum += (temp * temp);
	}

	return (sum + m_bias);
}
