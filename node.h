#ifndef NODE_H
#define NODE_H
#include "Eigen/Dense"


class Node
{
    public:
	Eigen::VectorXd allele;
	double fitness;
	bool isEvaluated;
	int length;

	Node(int len)
	{
		length = len;
		isEvaluated = false;
		allele.setZero(length);
		fitness = 0.0;
	}
	Node()
	{
		isEvaluated = false;
		fitness = 0.0;
	}

	double getFitness()
	{
		return fitness;
	}

	void setFitness(double val)
	{
		if(!isEvaluated)
		{
			fitness = val;
			isEvaluated = true;
		}
		else
		    fitness = val;
		return ;
	}

	Node &operator=(const Node rhs)
	{
		if(this == &rhs)
		    return *this;
		allele = rhs.allele;
		fitness = rhs.fitness;
		isEvaluated = rhs.isEvaluated;
		length = rhs.length;
		return *this;
	}
};

#endif
