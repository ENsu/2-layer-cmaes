#include <cassert>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include "group.h"
#include "benchmark.h"
#include "cmaes.h"

using namespace std;
using Eigen::VectorXd;

extern test_func *testFunc;
extern int NFE;


int main( int argc, char *argv[] )
{
	// process command line
	/*if( argc != 4 )   
	{
		cout << "Usage: " << argv[0] << " inputfile outputfile outputfile2" << endl;
		cout << "       Please read the README file." << endl;
		exit(1);
	}*/
	// read parameters from input file
//	ifstream infile( argv[1] ); 
//	read_parameters( infile );
//	infile.close();
//	// open output file
//	ofstream outfile( argv[2] );  
//	ofstream outfile2( argv[3] );
//	// initilalize random number generator
//	RANDOM.randomize( parameter::seed );
//	RANDOM.randomize();
	initial();		//upload the value of /supportData/fbias_data.txt into m_biases; sqrt(((double )i) + 1.0) into  m_iSqrt
	testFunc=testFunctionFactory(5,2);


	double tmp[2];
	tmp[0] = 0;
	tmp[1] = 0;
	double result = testFunc->f(tmp , 2);
	cout << result << endl;

	VectorXd a(2);
	a << 0, 0;
	Node node(2, a);
	cout << node.getFitness() << endl;
	cout << NFE << endl;

	list<Node> tmp_list;
	for(int i=0; i<10; i++)
	{
		a << i, i;
		Node node(2, a);
		cout << node.getFitness() << endl;
		tmp_list.push_back(node);
	}
	GROUP group(1,2, tmp_list);
	cout << group.getMin() << endl;
	cout << group.getMean() << endl;
	cout << group.getMax() << endl;
	cout <<	group.getVariance() << endl;
	cout <<	group.getUCBVal(10, group.getMin(), group.getMax()) << endl;
	cout <<	group.getSize() << endl;
	group.sort_node();
	group.print();





	// run the ECGA algorithm
//	cout << "ECGA done" << endl;
//	cout << "Evaluated functions:" << parameter::eva_fun << endl;
//	outfile.close();
//	outfile2.close();
	return 0;
}
