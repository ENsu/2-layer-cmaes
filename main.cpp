#include "benchmark.h"
#include <cassert>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include <node.h>

using namespace std;

test_func *testFunc;


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
	testFunc=testFunctionFactory(1,2);
	double tmp[2];
	tmp[0] = -3.9311900e+001;
	tmp[1] = 5.8899900e+001;
	double result = testFunc->f(tmp , 1);
	cout << result << endl;

	Node node
	// run the ECGA algorithm
//	cout << "ECGA done" << endl;
//	cout << "Evaluated functions:" << parameter::eva_fun << endl;
//	outfile.close();
//	outfile2.close();
	return 0;
}
