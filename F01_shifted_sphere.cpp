#include"F01_shifted_sphere.h"
#include"benchmark.h"

F01_shifted_sphere::F01_shifted_sphere(int dimension, double bias)
{
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME1;
	
	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_z = new double[m_dimension];

	// Load the shifted global optimum
	loadRowVectorFromFile(DEFAULT_FILE_DATA1, m_dimension, m_o);
}


F01_shifted_sphere::F01_shifted_sphere(int dimension, double bias, string file_data)
{
	m_dimension = dimension;
        m_bias = bias;
        m_func_name = FUNCTION_NAME1;

	// Note: dimension starts from 0
	m_o = new double[m_dimension];
	m_z = new double[m_dimension];

	// Load the shifted global optimum
	loadRowVectorFromFile(file_data, m_dimension, m_o);
}

double F01_shifted_sphere::f(double *x,int length)
{
	double result = 0.0;

	shift(m_z, x, m_o,length);

	result = sphere(m_z,length);

	result += m_bias;
	return (result);
}
