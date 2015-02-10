#ifndef _HCJob
#define _HCJob
class HCJob {
public:
	// Number of basic functions
	int num_func;
	// Number of dimensions
	int num_dim;

	// Predefined constant
	double C;
	// Coverage range for each basic function
	double *sigma;
	// Biases for each basic function
	double *biases;
	// Stretch / compress each basic function
	double *lambda;
	// Estimated fmax
	double *fmax;
	// Shift global optimum for each basic function
	double **o;
	// Linear transformation matrix for each basic function
	double ***M;

	// Working areas to avoid memory allocation operations
	double *w;
	double **z;
	double **zM;

	HCJob() {
		// Nothing
		// This class is just a place holder.
	}

	virtual double basic_func(int func_no, double *x,int length)=0;
};
#endif
