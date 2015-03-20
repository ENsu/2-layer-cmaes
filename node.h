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
		bool isVirtual;
		int group_id; //if its a virtual node, it represent of a group

		Node(VectorXd input_allele);
		Node();
		double getFitness();
		void print();
		bool outofBound();
		void intoBound();
		Node& operator=(const Node rhs);
};

#endif
