#include <cassert>
#include <iostream>
#include <cctype>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <string.h>

using namespace std;

vector<double> L;         // a subset

void read_data( ifstream &in )
{
	char  s[100] = "   ";
	double fitness=0;

	while( strcmp( s, "END" ) != 0 ) 
	{
		in>>s>>s>>fitness;
		L.push_back(fitness);
		cout<<"fitness="<<fitness<<endl;
	}
	L.pop_back();
} // end of read_parameters

int compare (const void * a, const void * b)
{
  return ( *(double*)a - *(double*)b );
}



int main( int argc, char *argv[] )
{
	if( argc != 3 )   
	{
		cout << "Usage: " << argv[0] << " inputfile outputfile" << endl;
		cout << "  Statistic the inputfile and output its result in outputfile" << endl;
		exit(1);
	}
	/////////////////////////////////////////
	ifstream infile( argv[1] ); 
	read_data( infile );
	infile.close();
	//////////////////////////////////////////
	for(int i=0; i<L.size(); i++) cout<<L[i]<<"   ";
	cout<<endl;
	double sum=0;
	for(int i=0; i<L.size(); i++){
		sum += L[i];
	}
	double avg=sum/L.size();
	cout<<"avg = "<<avg<<endl;
	double stdev=0;
	for(int i=0; i<L.size(); i++){
		stdev += ((L[i] - avg) * (L[i] - avg));
	}
	stdev = stdev / (L.size());
	stdev = sqrt(stdev);
	cout<<"stdev = "<<stdev<<endl;
	qsort (&L[0], L.size(), sizeof(double), compare);
	double mid = L[(int)(L.size()/2)];
	cout<<"mid num = "<<mid<<endl;

	ofstream outfile( argv[2] );  
	outfile<<"avg = "<<avg;
	outfile<<"	stdev = "<<stdev;
	outfile<<"	mid num = "<<mid<<endl;
	outfile.close();
	return 0;
}
