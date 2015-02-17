#ifndef NODE_H
#define NODE_H

#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class Node
{
    public:
		VectorXd allele;
		double fitness;
		bool isEvaluated;

		Node(VectorXd input_allele);
		Node();
		double getFitness();
		void print();
		bool outofBound();
		Node& operator=(const Node rhs);
};

#endif
