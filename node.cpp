#include "Eigen/Dense"
#include "benchmark.h"
#include "global.h"
#include "node.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

test_func *testFunc;

Node::Node(int dimension, VectorXd input_allele)
{
	dim = dimension;
	isEvaluated = false;
	allele = input_allele;
	fitness = 0.0;
}

Node::Node(int dimension)
{
	dim = dimension;
	isEvaluated = false;
	allele.setZero(dim); 
	fitness = 0.0;
}

double Node::getFitness()
{
	if(isEvaluated)
		return fitness;
	else{
		double tmp[dim];
		for(int i = 0 ; i < dim ; i++)
			tmp[i] = allele(i);
    	fitness = testFunc->f(tmp , dim);
    	isEvaluated = true;
    	NFE++;
    	return fitness;
	}
}

void Node::print()
{
	cout << "allele: " << endl << this->allele << endl;
	cout << "fitness: " << this->getFitness() << endl;
}
/*Node::Node &operator=(const Node rhs)
{
	if(this == &rhs)
	    return *this;
	allele = rhs.allele;
	fitness = rhs.fitness;
	isEvaluated = rhs.isEvaluated;
	dim = rhs.dim;
	return *this;
}*/

