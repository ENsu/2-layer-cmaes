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
	dimension = 2;
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
	cout << "=========" << lambda << "======" << endl;
	cout << "=========" << mu << "======" << endl;

	for(int i=0; i<mu; i++)
	{
		Eigen::VectorXd a(dimension);
		for(int i=0; i<dimension; i++)
			a(i) = RANDOM.uniform(lowerbound , upperbound);
		tmp_node_list.push_back(Node(a));
	}

	Group my_group = Group(tmp_node_list);

	CMAES cmaes(mu, lambda, sigma, &my_group);
	Node best_node = my_group.get_best_node();

	Eigen::MatrixXd covar(2,2);
	//cout << X << endl;
	//Eigen::MatrixXd A = X * X.transpose();
	//cout << A << endl;
	covar <<  693.847, -743.084, -743.084, 795.814;
	//cout << covar << endl;
	Eigen::SelfAdjointEigenSolver<MatrixXd> es(covar / covar.norm());
	//cout << es. << endl;
    cout << es.operatorInverseSqrt() << endl;
	//while(best_node.getFitness() > -449)
	//{

	for(int i=0; i<1000; i++)
	{
		cout << "generation: " << i << endl;
		cmaes.run();
		my_group = *cmaes.group;
		best_node = my_group.get_best_node();
		best_node.print();
		cout <<"yw: " << cmaes.yw << endl;
		cout <<"ps: " << cmaes.ps << endl;
		cout <<"pc: " << cmaes.pc << endl;
		cout <<"covar: " << cmaes.covar << endl;
		cout <<"sigma: " << cmaes.sigma << endl;
	}
	//my_group.print();
	cout << "NFE:" << NFE << endl;
	//}
	//my_group.print();	
	//create population


	//group them

	//create bar

	//check condition

	//allocate resource

	return 0;
}
