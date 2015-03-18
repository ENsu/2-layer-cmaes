#include <cassert>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include "benchmark.h"
#include "cmaes.h"
#include "group.h"
#include "node.h"
#include "global.h"
#include "util.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

extern test_func *testFunc;

randomG RANDOM;


void update_gen_var(Group *groupArry, int arryLen, int& gen_node_num, double& gen_min, double& gen_max)
{
	gen_node_num = 0;
	gen_max = -1 * DBL_MAX;
	gen_min = DBL_MAX;

	for(int i=0; i<arryLen; i++)
	{
		gen_node_num += groupArry[i].getSize();
		if(gen_max < groupArry[i].getMax())
			gen_max = groupArry[i].getMax();
		if(gen_min > groupArry[i].getMin())
			gen_min = groupArry[i].getMin();
	}
	return;
}

list<Node> random_sample(Group *groupArry, int arryLen, int each)
{
	list<Node> nodeL;
	for(int i=0; i<arryLen; i++)
	{
		list<Node> sampleL = groupArry[i].random_pick(each);
		nodeL.insert(nodeL.end(), sampleL.begin(), sampleL.end());
	}
	return nodeL;
}

list<Node> group_to_virtual_nodes(CMAES *cmaesArry, Group *groupArry, int arryLen)
{
	int gen_node_num = 0; 
	double gen_max = -1 * DBL_MAX;
	double gen_min = DBL_MAX;

	//update value of gen_node_num, gen_min, gen_max according to the groups
	update_gen_var(groupArry, arryLen, gen_node_num, gen_min, gen_max);

	list<Node> virtualNodes;
	for(int i=0; i<arryLen; i++)
	{
		//create a virtual node according to the layer 1 group info
		Node virtual_node = groupArry[i].get_mean_node(cmaesArry[i].weight);
		virtual_node.fitness = groupArry[i].getUCBVal(gen_node_num, gen_min, gen_max);
		virtual_node.isEvaluated = true;
		virtual_node.isVirtual = true;
		virtual_node.group_id = i;  // the arry index of input arrys
		virtualNodes.push_back(virtual_node);
	}


	return virtualNodes;
}

int main( int argc, char *argv[] )
{
	/*if( argc != 4 )   
	{
		cout << "Usage: " << argv[0] << " inputfile outputfile outputfile2" << endl;
		cout << "       Please read the README file." << endl;
		exit(1);
	}*/

// ---------------------initialization------------------------
	//initialize the test function
	RANDOM.randomize(time(NULL));
	initial();

	int generation = 0;

	// global parameters
	dimension = 2;
	funNum = 4;
	double lowerbound = domainlowbound[funNum-1];
	double upperbound = domainupbound[funNum-1];
	int lambda = int(4 + 3 * log(dimension));
	cout << "value of lambda: " << lambda << endl;
	int mu = int(lambda/2);
	cout << "value of mu: " << mu << endl;
	double sigma = (upperbound - lowerbound) * 0.3;

if(true) // pure cmaes
{
	generation = 0;
	list<Node> tmp_node_list;
	tmp_node_list.clear();
	for(int i=0; i<mu; i++)
	{
		Eigen::VectorXd a(dimension);
		for(int j=0; j<dimension; j++)
			a(j) = RANDOM.uniform(lowerbound , upperbound);
		tmp_node_list.push_back(Node(a));
	}

	Group my_group = Group(tmp_node_list);
	CMAES cmaes(mu, lambda, sigma, &my_group);
	Node best_node = my_group.get_best_node();

	for(int i=0; i<1000; i++)
	{
		generation ++;
		cout << "generation: " << generation << endl;
		cmaes.run();
		my_group = *cmaes.group;
		best_node = my_group.get_best_node();
		best_node.print();
	}
}

else if(false) // my cmaes 
{
	// some value used
	int gen_node_num = 0; // used to sum up the node number of current generation
	double gen_max = -1 * DBL_MAX;
	double gen_min = DBL_MAX;
	int sample_each_group = 2;

	// initiate the beginning population.
	Eigen::VectorXd *vectors = new Eigen::VectorXd[mu * mu];
	for(int i=0; i< mu * mu; i++)
	{
		Eigen::VectorXd a(dimension);
		for(int j=0; j<dimension; j++)
			a(j) = RANDOM.uniform(lowerbound , upperbound);
		vectors[i] = a;
	}
	// we store the beginning population in an VectorXD arry
//-------------------------------------------------------------------


//---------------------preparing for iteration----------------------- 
	// run k-means and group the VectorXd arry into 
	// array of VectorXd list, each list represent a vector group
	list<Eigen::VectorXd> *vectorGroup;
	vectorGroup = k_means(vectors, mu*mu, mu);
	//delete vectors;  //error: pointer being freed was not allocated

	// create first layer groups and cmaes and record their UCB value...
	CMAES *Layer1CMAESs = new CMAES[mu];
	Group *Layer1Groups = new Group[mu];

	//we convert the VectorXD list into Node list, assign into a group...
	//binding with a cmaes machine

	list<Node> tmp_node_list; 
	for(int i=0; i<mu; i++)
	{
		list<Eigen::VectorXd>::iterator iter;
		for(iter = vectorGroup[i].begin() ; iter != vectorGroup[i].end() ; ++iter)
			tmp_node_list.push_back(Node(*iter));
		Layer1Groups[i] = Group(tmp_node_list);

		Layer1CMAESs[i] = CMAES(mu, lambda, sigma, &Layer1Groups[i]);
		Layer1CMAESs[i].tune_node_num();
		tmp_node_list.clear();
	}
	// create a second layer group and cmaes...
	list<Node> Layer1VirtualNodes = group_to_virtual_nodes(Layer1CMAESs, Layer1Groups, mu);

	Group Layer2Group(Layer1VirtualNodes);
	CMAES Layer2CMAES(mu, lambda, sigma, &Layer2Group);

//-------------------------------atart iteration-----------------------------
//	while(!termination)

	Group *GroupsPool = new Group[lambda];
	CMAES *CMAESsPool = new CMAES[lambda];

	for(int k=0; k<2; k++)
	{
		//update value of gen_node_num, gen_min, gen_max according to the groups

		// calculate second layer UCB value...
		// first we random sample # nodes in each group
		update_gen_var(Layer1Groups, mu, gen_node_num, gen_min, gen_max);
		list<Node> rSampleNode = random_sample(Layer1Groups, mu, sample_each_group);
		double UCB_global = Group(rSampleNode).getUCBVal(gen_node_num, gen_min, gen_max);

		// find the maximum UCB value among both 1st and 2nd layer...
		double UCBmax = -1 * DBL_MAX;
		int UCBmax_id = 0;

		list<Node>::iterator iter;
		for(iter = Layer2Group.Nodes.begin(); iter != Layer2Group.Nodes.end(); ++iter)
			if(UCBmax < iter->getFitness())
			{
				UCBmax = iter->getFitness();
				UCBmax_id = iter->group_id;
			}
		//cout << "UCBmax = " << UCBmax << endl;
		//cout << "UCBmax_id = " << UCBmax_id << endl;
		//cout << "UCB_global = " << UCB_global << endl;
		if(UCBmax > UCB_global)
		{
			cout << "into 1 layer cmaes process" << endl;
			Layer1CMAESs[UCBmax_id].run();
			Layer1Groups[UCBmax_id] = *Layer1CMAESs[UCBmax_id].group;
			Layer1VirtualNodes = group_to_virtual_nodes(Layer1CMAESs, Layer1Groups, mu);
			Group Layer2Group(Layer1VirtualNodes);
			Layer2CMAES.update_value(Layer2Group);
		}
		else
		{
			cout << "into 2 layer cmaes process" << endl;
			// layer 2 cmaes run...
			// sample (lambda - mu) real nodes according to virtaul nodes
			list<Node> virtualNodesPool = Layer2CMAES.sample_node(lambda - mu);

			for(int i=0; i<mu; i++)
			{
				GroupsPool[i] = Layer1Groups[i];
				CMAESsPool[i] = Layer1CMAESs[i];
			}

			// add new sampled node into group and bind to a new cmaes machine
			int i = mu+1;
			list<Node>::iterator iter;
			for(iter = virtualNodesPool.begin(); iter != virtualNodesPool.end(); ++iter)
			{
				list<Node> new_list;
				new_list.push_back(*iter);
				GroupsPool[i] = Group(new_list);
				CMAESsPool[i] = CMAES(mu, lambda, sigma, &GroupsPool[i]);
				CMAESsPool[i].tune_node_num();
				i++;
			}
			assert(i == lambda);

			// make virtaul node from these newly sampled node
			update_gen_var(GroupsPool, mu, gen_node_num, gen_min, gen_max);
			list<Node> sample_nodes = group_to_virtual_nodes(CMAESsPool, GroupsPool, lambda);
			Group newLayer2Group = Group(sample_nodes);

			// select cirtain virtual node
			newLayer2Group.sort_node();
			newLayer2Group.truncate_size(mu);

			// update relative values
			Layer2CMAES.update_value(newLayer2Group);

			// put the layer 2 result back into layer 1 arry.
			i = 0;
			for(iter = newLayer2Group.Nodes.begin(); iter != newLayer2Group.Nodes.end(); ++iter)
			{
				int id = iter->group_id;
				Layer1Groups[i] = GroupsPool[id];
				Layer1CMAESs[i] = CMAESsPool[id];
				i++;
			}
			assert(i == mu);
		}
	}
}

// --------------------------finish iteration-----------------------------
/*
	Group my_group = Group(tmp_node_list);

	CMAES cmaes(mu, lambda, sigma, my_group);
	Node best_node = my_group.get_best_node();

	for(int i=0; i<500; i++)
	{
		generation ++;
		//cout << "generation: " << generation << endl;
		cmaes.run();
		my_group = *cmaes.group;
		best_node = my_group.get_best_node();
		if(fabs(best_node.getFitness() + 450) < 1E-10 )
			break;
	}
	//my_group.print();
	cout << "best fitness: " << best_node.getFitness() << endl;
	cout << "generation: " << generation << endl;
	cout << "NFE: " << NFE << endl;
	//}
	//my_group.print();	
	//create population

	//group them

	//create bar

	//check condition

	//allocate resource
*/
	return 0;
}
