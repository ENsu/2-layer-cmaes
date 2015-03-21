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
		double getMax();
		double getVariance();
		double getNormVariance(double global_min, double global_max);
		double getUCBVal(int total_played, double global_min, double global_max);
		int getSize();

		Node get_mean_node();	
		Node get_mean_node(Eigen::VectorXd weight);	
		Node get_best_node();
		Node get_worst_node();
		Eigen::MatrixXd node_matrix();
		
		void add_node(Node a);
		void drop_node(Node a);
		void replace_node(Node a , Node b);

		void truncate_size(int size);
		void sort_node();
		void sort_node_descend();
		list<Node> random_pick(int pick_num);
		void print();

		Group& operator=(const Group rhs);
};


#endif
