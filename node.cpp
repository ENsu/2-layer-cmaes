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
	isVirtual = false;
	allele = input_allele;
	fitness = 0.0;
}

Node::Node()
{
	isEvaluated = false;
	isVirtual = false;
	allele.setZero(dimension); 
	fitness = 0.0;
}

double Node::getFitness()
{
	if(isEvaluated)
		return fitness;
	else if(!isVirtual)
	{
		double tmp[dimension];
		for(int i = 0 ; i < dimension ; i++)
		{
			tmp[i] = allele(i);
		}
		testFunc=testFunctionFactory(funNum,dimension);
    	fitness = testFunc->f(tmp , dimension);
    	isEvaluated = true;
    	NFE++;
    	return fitness;
	}
	else
	{
		cout << "Should not call a virtual node without fitness assignment" << endl;
		assert(1);
		return NULL;
	}
}

void Node::print()
{
	//cout << "-----------------Node----------------";
	cout << "-----Node info-----" << endl;
	cout << "isEvaluated: " << isEvaluated << endl;
	cout << "isVirtual: " << isVirtual << endl;
	cout << "allele: " << endl << allele << endl;
	cout << "fitness: " << getFitness() << endl;
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

void Node::intoBound()
{
	for(int i=0 ; i<dimension; i++)
	{
		if(allele(i) > solupbound[funNum-1])
			allele(i) = solupbound[funNum-1];
		else if(allele(i) < sollowbound[funNum -1])
			allele(i) = sollowbound[funNum-1];
	}
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

