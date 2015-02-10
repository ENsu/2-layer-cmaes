#include"F17_rotated_hybrid_composition_1_noise.h"
extern randomG random_;
// Constructors
F17_rotated_hybrid_composition_1_noise::F17_rotated_hybrid_composition_1_noise (int dimension, double bias) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME17;

	// Note: dimension starts from 0
	m_o = new double*[NUM_FUNC];
	for(int i=0;i<NUM_FUNC;++i)
		m_o[i] = new double[m_dimension];
	m_M = new double**[NUM_FUNC];
	for(int i=0;i<NUM_FUNC;++i)
	{
		m_M[i] = new double*[m_dimension];
		for(int j=0;j<m_dimension;++j)
			m_M[i][j] = new double[m_dimension];
	}

	m_testPoint = new double[m_dimension];
	m_testPointM = new double[m_dimension];
	m_fmax = new double[NUM_FUNC];

	m_w = new double[NUM_FUNC];
	m_z = new double*[NUM_FUNC];
	for(int i=0;i<NUM_FUNC;++i)
		m_z[i] = new double[m_dimension];
	m_zM = new double*[NUM_FUNC];
	for(int i=0;i<NUM_FUNC;++i)
		m_zM[i] = new double[m_dimension];

	// Load the shifted global optimum
	loadMatrixFromFile(DEFAULT_FILE_DATA17, NUM_FUNC, m_dimension, m_o);
	// Load the matrix
	char tmp[4];
	sprintf(tmp,"%d",m_dimension);
	loadNMatrixFromFile(DEFAULT_FILE_MX_PREFIX17 + tmp + DEFAULT_FILE_MX_SUFFIX17, NUM_FUNC, m_dimension, m_dimension, m_M);

	// Initialize the hybrid composition job object
	theJob.num_func = NUM_FUNC;
	theJob.num_dim = m_dimension;
	theJob.C = 2000.0;
	theJob.sigma = m_sigma17;
	theJob.biases = m_func_biases17;
	theJob.lambda = m_lambda17;
	theJob.o = m_o;
	theJob.M = m_M;
	theJob.w = m_w;
	theJob.z = m_z;
	theJob.zM = m_zM;
	// Calculate/estimate the fmax for all the functions involved
	for (int i = 0 ; i < NUM_FUNC ; i ++) {
		for (int j = 0 ; j < m_dimension ; j ++) {
			m_testPoint[j] = (5.0 / m_lambda17[i]);
		}
		rotate(m_testPointM, m_testPoint, m_M[i],m_dimension);
		m_fmax[i] = fabs(theJob.basic_func(i, m_testPointM,m_dimension));
	}
	theJob.fmax = m_fmax;
}

F17_rotated_hybrid_composition_1_noise::F17_rotated_hybrid_composition_1_noise (int dimension, double bias, string file_data, string file_m) {
	m_dimension = dimension;
	m_bias = bias;
	m_func_name = FUNCTION_NAME17;

	// Note: dimension starts from 0
	m_o = new double*[NUM_FUNC];
	for(int i=0;i<NUM_FUNC;++i)
		m_o[i] = new double[m_dimension];
	m_M = new double**[NUM_FUNC];
	for(int i=0;i<NUM_FUNC;++i)
	{
		m_M[i] = new double*[m_dimension];
		for(int j=0;j<m_dimension;++j)
			m_M[i][j] = new double[m_dimension];
	}

	m_testPoint = new double[m_dimension];
	m_testPointM = new double[m_dimension];
	m_fmax = new double[NUM_FUNC];

	m_w = new double[NUM_FUNC];
	m_z = new double*[NUM_FUNC];
	for(int i=0;i<NUM_FUNC;++i)
		m_z[i] = new double[m_dimension];
	m_zM = new double*[NUM_FUNC];
	for(int i=0;i<NUM_FUNC;++i)
		m_zM[i] = new double[m_dimension];

	// Load the shifted global optimum
	loadMatrixFromFile(file_data, NUM_FUNC, m_dimension, m_o);
	// Load the matrix
	loadNMatrixFromFile(file_m, NUM_FUNC, m_dimension, m_dimension, m_M);

	// Initialize the hybrid composition job object
	theJob.num_func = NUM_FUNC;
	theJob.num_dim = m_dimension;
	theJob.C = 2000.0;
	theJob.sigma = m_sigma17;
	theJob.biases = m_func_biases17;
	theJob.lambda = m_lambda17;
	theJob.o = m_o;
	theJob.M = m_M;
	theJob.w = m_w;
	theJob.z = m_z;
	theJob.zM = m_zM;
	// Calculate/estimate the fmax for all the functions involved
	for (int i = 0 ; i < NUM_FUNC ; i ++) {
		for (int j = 0 ; j < m_dimension ; j ++) {
			m_testPoint[j] = (5.0 / m_lambda17[i]);
		}
		rotate(m_testPointM, m_testPoint, m_M[i],m_dimension);
		m_fmax[i] = fabs(theJob.basic_func(i, m_testPointM,m_dimension));
	}
	theJob.fmax = m_fmax;
}

double MyHCJob17::basic_func(int func_no, double *x,int length) {
	double result = 0.0;
	switch(func_no) {
		case 0:
		case 1:
			result = rastrigin(x,length);
			break;
		case 2:
		case 3:
			result = weierstrass(x,length);
			break;
		case 4:
		case 5:
			result = griewank(x,length);
			break;
		case 6:
		case 7:
			result = ackley(x,length);
			break;
		case 8:
		case 9:
			result = sphere(x,length);
			break;
		default:
			printf("func_no is out of range.");
			exit(-1);
	}
	return (result);
}

// Function body
double F17_rotated_hybrid_composition_1_noise::f(double *x,int length) {

	double result = 0.0;

	result = hybrid_composition(x, theJob,length);

	result += m_bias;

	// NOISE
	// Comment the next line to remove the noise
#ifndef no_noise	
	result *= (1.0 + 0.2 * fabs(random_.normal01()));
#endif

	return (result);
}