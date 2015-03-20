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

#define CMAES_MODE "2LCMAES"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

extern test_func *testFunc;

randomG RANDOM;


void update_gen_var(CMAES *cmaesArry, int arryLen, int& gen_node_num, double& gen_min, double& gen_max)
{
	gen_node_num = 0;
	gen_max = -1 * DBL_MAX;
	gen_min = DBL_MAX;

	for(int i=0; i<arryLen; i++)
	{
		gen_node_num += cmaesArry[i].group->getSize();
		if(gen_max < cmaesArry[i].group->getMax())
			gen_max = cmaesArry[i].group->getMax();
		if(gen_min > cmaesArry[i].group->getMin())
			gen_min = cmaesArry[i].group->getMin();
	}
	return;
}

list<Node> random_sample(CMAES *cmaesArry, int arryLen, int each)
{
	list<Node> nodeL;
	for(int i=0; i<arryLen; i++)
	{
		list<Node> sampleL = cmaesArry[i].group->random_pick(each);
		nodeL.insert(nodeL.end(), sampleL.begin(), sampleL.end());
	}
	return nodeL;
}

list<Node> group_to_virtual_nodes(CMAES *cmaesArry, int arryLen)
{
	int gen_node_num = 0; 
	double gen_max = -1 * DBL_MAX;
	double gen_min = DBL_MAX;

	//update value of gen_node_num, gen_min, gen_max according to the groups
	update_gen_var(cmaesArry, arryLen, gen_node_num, gen_min, gen_max);

	list<Node> virtualNodes;
	for(int i=0; i<arryLen; i++)
	{
		//create a virtual node according to the layer 1 group info
		Node virtual_node = cmaesArry[i].group->get_mean_node(cmaesArry[i].weight);
		virtual_node.fitness = cmaesArry[i].group->getUCBVal(gen_node_num, gen_min, gen_max);
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
	funNum = 1;
	double lowerbound = domainlowbound[funNum-1];
	double upperbound = domainupbound[funNum-1];
	int lambda = int(4 + 3 * log(dimension));
	cout << "value of lambda: " << lambda << endl;
	int mu = int(lambda/2);
	cout << "value of mu: " << mu << endl;
	double sigma = (upperbound - lowerbound) * 0.3;

if(CMAES_MODE == "CMAES") // pure cmaes
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
		cout << "fitness: " << best_node.fitness << endl;
		//best_node.print();
		if(cmaes.termination == true)
		{
			cout << "ending sigma value: " << cmaes.sigma << endl;
			break;
		}
	}
}

else if(CMAES_MODE == "2LCMAES") // my cmaes 
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

	//we convert the VectorXD list into Node list, assign into a group...
	//binding with a cmaes machine

	list<Node> tmp_node_list; 
	for(int i=0; i<mu; i++)
	{
		list<Eigen::VectorXd>::iterator iter;
		for(iter = vectorGroup[i].begin() ; iter != vectorGroup[i].end() ; ++iter)
			tmp_node_list.push_back(Node(*iter));
		Group *layer1group = new Group(tmp_node_list);
		Layer1CMAESs[i] = CMAES(mu, lambda, sigma, layer1group);
		Layer1CMAESs[i].tune_node_num();
		tmp_node_list.clear();
	}
	// create a second layer group and cmaes...
	list<Node> Layer1VirtualNodes = group_to_virtual_nodes(Layer1CMAESs, mu);

	Group *layer2group = new Group(Layer1VirtualNodes);
	CMAES Layer2CMAES(mu, lambda, sigma, layer2group);

//-------------------------------atart iteration-----------------------------
//	while(!termination)

	CMAES *CMAESsPool = new CMAES[lambda];

	for(int k=0; k<1000; k++)
	{
		//update value of gen_node_num, gen_min, gen_max according to the groups

		// calculate second layer UCB value...
		// first we random sample # nodes in each group
		update_gen_var(Layer1CMAESs, mu, gen_node_num, gen_min, gen_max);
		list<Node> rSampleNode = random_sample(Layer1CMAESs, mu, sample_each_group);
		double UCB_global = Group(rSampleNode).getUCBVal(gen_node_num, gen_min, gen_max);

		// find the maximum UCB value among both 1st and 2nd layer...
		double UCBmax = -1 * DBL_MAX;
		int UCBmax_id = 0;

		list<Node>::iterator iter;
		for(iter = Layer2CMAES.group->Nodes.begin(); iter != Layer2CMAES.group->Nodes.end(); ++iter)
		{
			assert(iter->isEvaluated == true);
			assert(iter->isVirtual == true);
			if(UCBmax < iter->getFitness())
			{
				UCBmax = iter->getFitness();
				UCBmax_id = iter->group_id;
			}
		}
		//cout << "UCBmax = " << UCBmax << endl;
		//cout << "UCBmax_id = " << UCBmax_id << endl;
		//cout << "UCB_global = " << UCB_global << endl;
		if(false)//(UCBmax > UCB_global)
		{
			// cout << "into 1 layer cmaes process" << endl;
			Layer1CMAESs[UCBmax_id].run();
			if(Layer1CMAESs[UCBmax_id].termination == true)
				break;
			Layer1VirtualNodes = group_to_virtual_nodes(Layer1CMAESs, mu);
			layer2group = new Group(Layer1VirtualNodes);
			Layer2CMAES.update_value(*layer2group); // this step is terrible
			Layer2CMAES.group = layer2group;
		}
		else
		{
			cout << "into 2 layer cmaes process" << endl;
			// layer 2 cmaes run...
			// sample (lambda - mu) real nodes according to virtaul nodes
			list<Node> virtualNodesPool = Layer2CMAES.sample_node(lambda - mu);

			// set pools
			for(int i=0; i<mu; i++)
			{
				CMAESsPool[i] = Layer1CMAESs[i];
			}
			// add new sampled node into group and bind to a new cmaes machine
			int i = mu;
			list<Node>::iterator iter;
			for(iter = virtualNodesPool.begin(); iter != virtualNodesPool.end(); ++iter)
			{
				list<Node> new_groups_nodes;
				new_groups_nodes.push_back(*iter);
				Group *tmp = new Group(new_groups_nodes);
				CMAESsPool[i] = CMAES(mu, lambda, sigma, tmp);
				CMAESsPool[i].tune_node_num();
				i++;
			}
			assert(i == lambda);

			// make virtaul node from these newly created group
			update_gen_var(CMAESsPool, mu, gen_node_num, gen_min, gen_max);
			list<Node> sample_nodes = group_to_virtual_nodes(CMAESsPool, lambda);
			Group *newlayer2group = new Group(sample_nodes);

			// select cirtain virtual node
			newlayer2group->sort_node();
			newlayer2group->truncate_size(mu);
			// update relative values
			Layer2CMAES.group = newlayer2group;
			Layer2CMAES.update_value(*newlayer2group);

			// put the layer 2 result back into layer 1 arry.
			i = 0;
			for(iter = newlayer2group->Nodes.begin(); iter != newlayer2group->Nodes.end(); ++iter)
			{
				int id = iter->group_id;
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
