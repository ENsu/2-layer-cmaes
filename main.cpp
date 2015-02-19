#include <cassert>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include "benchmark.h"
#include "cmaes.h"
#include "group.h"
#include "node.h"
#include "global.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

extern test_func *testFunc;

randomG RANDOM;

int main( int argc, char *argv[] )
{
	/*if( argc != 4 )   
	{
		cout << "Usage: " << argv[0] << " inputfile outputfile outputfile2" << endl;
		cout << "       Please read the README file." << endl;
		exit(1);
	}*/
	// parameters
	dimension = 10;
	funNum = 1;
	double lowerbound = domainlowbound[funNum-1];
	double upperbound = domainupbound[funNum-1];
	// initialize the problem
	RANDOM.randomize(time(NULL));
	initial();
	testFunc=testFunctionFactory(funNum,dimension);

	//Node *nodes = new Node[popSize];
	list<Node> tmp_node_list;
	int lambda = int(4 + 3 * log(dimension));
	int mu = int(lambda/2);
	double sigma = (upperbound - lowerbound) * 0.3;
	cout << "======lambda: " << lambda << "======" << endl;
	cout << "==========mu: " << mu << "======" << endl;
for(int k=0; k<100; k++)
{
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

	for(int i=0; i<500; i++)
	{
		generation ++;
		cout << "generation: " << generation << endl;
		cmaes.run();
		my_group = *cmaes.group;
		best_node = my_group.get_best_node();
		if(fabs(best_node.getFitness() + 450) < 1E-10 )
			break;
	}
	//my_group.print();
	cout << "generation: " << generation << endl;
	cout << "NFE: " << NFE << endl;

}
	//}
	//my_group.print();	
	//create population

	//group them

	//create bar

	//check condition

	//allocate resource

	return 0;
}
