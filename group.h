#ifndef GROUP_H
#define GROUP_H

#include <iostream>
#include <list>
#include <utility>
#include "node.h"

using namespace std;

class Group
{
	public:
		list<Node> Nodes;
		int AllNodeNum;

		Group(list<Node> new_Nodes);
		Group();
		~Group();

		double getMin();
		double getMean();
		double getMean(double *weight, int weight_length);
		double getMax();
		double getVariance();
		double getNormVariance(double global_min, double global_max);
		Eigen::MatrixXd getCov();
		Eigen::MatrixXd getCov(double *weight, int weight_length);
		double getUCBVal(int total_played, double global_min, double global_max);
		int getSize();

		Node get_mean_node();		
		Node get_best_node();
		Node get_worst_node();
		
		void add_node(Node a);
		void drop_node(Node a);
		void replace_node(Node a , Node b);

		void truncate_size(int size);
		void sort_node();
		void print();

		Group& operator=(const Group rhs);
};


#endif
