#include"F05_schwefel_global_opt_bound.h"
#include"benchmark.h"

// Constructors
F05_schwefel_global_opt_bound::F05_schwefel_global_opt_bound (int dimension, double bias) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = DEFAULT_FILE_DATA5;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_A = new double*[m_dimension];
	for(int i=0;i<m_dimension;++i)
		m_A[i]=new double[m_dimension];

	m_B = new double[m_dimension];
	m_z = new double[m_dimension];

	double **m_data = new double*[m_dimension+1];
	for(int i=0;i<m_dimension+1;++i)
		m_data[i]=new double[m_dimension];

	// Load the shifted global optimum
	loadMatrixFromFile(DEFAULT_FILE_DATA5, m_dimension+1, m_dimension, m_data);
	for (int i = 0 ; i < m_dimension ; i ++) {
		if ((i+1) <= ceil(m_dimension / 4.0))
			m_o[i] = -100.0;
		else if ((i+1) >= floor((3.0 * m_dimension) / 4.0))
			m_o[i] = 100.0;
		else
			m_o[i] = m_data[0][i];
	}
	for (int i = 0 ; i < m_dimension ; i ++) {
		for (int j = 0 ; j < m_dimension ; j ++) {
			m_A[i][j] = m_data[i+1][j];
		}
	}
	Ax(m_B, m_A, m_o,m_dimension);
}
F05_schwefel_global_opt_bound::F05_schwefel_global_opt_bound (int dimension, double bias, string file_data) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME5;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	
	m_A = new double*[m_dimension];
	for(int i=0;i<m_dimension;++i)
		m_A[i]=new double[m_dimension];

	m_B = new double[m_dimension];
	m_z = new double[m_dimension];

	double **m_data = new double*[m_dimension+1];
	for(int i=0;i<m_dimension+1;++i)
		m_data[i]=new double[m_dimension];

	// Load the shifted global optimum
	loadMatrixFromFile(file_data, m_dimension+1, m_dimension, m_data);
	for (int i = 0 ; i < m_dimension ; i ++) {
		if ((i+1) <= ceil(m_dimension / 4.0))
			m_o[i] = -100.0;
		else if ((i+1) >= floor((3.0 * m_dimension) / 4.0))
			m_o[i] = 100.0;
		else
			m_o[i] = m_data[0][i];
	}
	for (int i = 0 ; i < m_dimension ; i ++) {
		for (int j = 0 ; j < m_dimension ; j ++) {
			m_A[i][j] = m_data[i+1][j];
		}
	}
	Ax(m_B, m_A, m_o,m_dimension);
}

// Function body
double F05_schwefel_global_opt_bound::f(double *x,int length) {
		/*double temp1 = pow((x[0]-x[2]), 2);
		double temp2 = pow((x[2]-x[3]), 2);
		double temp3 = pow((x[1]-x[4]), 2);
		double temp5 = temp1 + temp2 + temp3;
		return temp5;*/
	double max = -HUGE_VAL; //negative inf

	Ax(m_z, m_A, x,length);

	for (int i = 0 ; i < m_dimension ; i ++) {
		double temp = fabs(m_z[i] - m_B[i]);
		if (max < temp)
			max = temp;
	}

	return (max + m_bias);
}
