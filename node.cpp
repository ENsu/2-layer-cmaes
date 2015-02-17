#include "Eigen/Dense"
#include "benchmark.h"
#include "node.h"
#include "global.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

test_func *testFunc;


Node::Node(VectorXd input_allele)
{
	isEvaluated = false;
	allele = input_allele;
	fitness = 0.0;
}

Node::Node()
{
	isEvaluated = false;
	allele.setZero(dimension); 
	fitness = 0.0;
}

double Node::getFitness()
{
	if(isEvaluated)
		return fitness;
	else{
		double tmp[dimension];
		for(int i = 0 ; i < dimension ; i++)
			tmp[i] = allele(i);
    	fitness = testFunc->f(tmp , dimension);
    	isEvaluated = true;
    	NFE++;
    	return fitness;
	}
}

void Node::print()
{
	//cout << "-----------------Node----------------";
	cout << "Node info: " << endl;
	cout << "allele: " << endl << this->allele << endl;
	cout << "fitness: " << this->getFitness() << endl;
	//cout << "---------------Node end--------------";
}

bool Node::outofBound()
{
	for(int i=0 ; i<dimension; i++)
	{
		if(allele(i) > solupbound[funNum-1] || allele(i) < sollowbound[funNum -1])
	    	return true;
	}
    return false;
}

Node& Node::operator=(const Node rhs)
{
	if(this == &rhs)
	    return *this;
	allele = rhs.allele;
	fitness = rhs.fitness;
	isEvaluated = rhs.isEvaluated;
	return *this;
}

