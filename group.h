#ifndef GROUP_H
#define GROUP_H

#include <iostream>
#include <list>
#include <utility>
#include "node.h"

using namespace std;

class GROUP
{
	public:
		int dim;
		int gID;
		list<Node> Nodes;
		int AllNodeNum;

		GROUP(int new_gid, int dimension, list<Node> new_Nodes);
		GROUP(int new_gid, int dimension);
		~GROUP();

		double getMin();
		double getMean();
		double getMean(double *weight);
		double getMax();
		double getVariance();
		double getNormVariance(double global_min, double global_max);
		Eigen::MatrixXd getCovar();
		double getUCBVal(int total_played, double global_min, double global_max);
		int getSize();

		Node get_mean_node();		
		Node get_best_node();
		Node get_worst_node();
		
		void add_node(Node a);
		void drop_node(Node a);
		void replace_node(Node a , Node b);
		void sort_node();
		void print();

};


#endif
