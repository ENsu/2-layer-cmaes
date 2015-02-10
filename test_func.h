#ifndef _test_func
#define _test_func
#include<string>
using namespace std;

class test_func {
protected:
	int	m_dimension;
	double	m_bias;

	string	m_func_name;
public:
/*	// Constructors
	test_func (int dimension, double bias) {
		m_dimension = dimension;
		m_bias = bias;
		m_func_name = "undefined";
	}
	test_func (int dimension, double bias,string func_name) {
		m_dimension = dimension;
		m_bias = bias;
		m_func_name = func_name;
	}
*/
	// Function body to be defined in the child classes
	virtual double f(double *x,int length)=0;

	// Property functions common for all child classes
	int dimension() {
		return (m_dimension);
	}

	double bias() {
		return (m_bias);
	}

	string name() {
		return (m_func_name);
	}
};
#endif
