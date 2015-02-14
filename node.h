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
		int dim;

		Node(int dimension, VectorXd input_allele);
		Node(int dimension);
		double getFitness();
		void print();
		//Node &operator=(const Node rhs);
};

#endif
