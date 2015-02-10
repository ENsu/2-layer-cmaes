#include<iostream>
#include<cmath>
#include<cassert>

#include "random.hpp"
//#define SIMPLE_SPRNG
//#include<sprng.h>
using namespace std;

//
// if a seed is not given, ask for one.
// seed must be an integer between 1 and m-1
//
/*void randomG::randomize()
{
	int seed;
	seed = make_sprng_seed();
	init_sprng(DEFAULT_RNG_TYPE,seed,SPRNG_DEFAULT);
}*/
void randomG::randomize( long seed )
{
	gsl_rng_set(r,seed);
}

//
// generate a random number distributed uniformly in 0..1
//
double randomG::uniform01()
{
//	return sprng();
	return gsl_rng_uniform_pos(r);
}

double randomG::uniform( double a, double b )
{
	assert( b >= a );
	return (b-a) * uniform01() + a;
}

int randomG::uniform( int a, int b )
{
	assert( b > a );
//	int i = int( (b-a+1) * uniform01() + a );
//	if( i>b ) i=b;
	int i=gsl_rng_uniform_int(r,b-a+1)+a;
	if( i>b ) i=b;
	return i;
}

//
// generate a random number with an exponential distribution with mean 'mu'
// note: log(x) gives the natural logarithm of x.
//
double randomG::exponential( double mu )
{
	double U = uniform01();
	return -mu*log(U);
}

//
// generate a Bernoulli trial with probability of success 'p'
//
bool randomG::bernoulli( double p )
{
	return ( uniform01() < p );
}

//
// generate a sample from a Geometric random variable with parameter 'p'.
// returns how many successes not including the first failure. 
// 'p' is the probability of failure.
//
long randomG::geometric( double p )
{
	double U = uniform01();
	return int( log(U)/log(1-p) );
}

//
// generate a sample from the standard normal distribution (mean=0, var=1).
// see book "Simulation Modeling and Analysis" by Law & Kelton, pag. 490-491.
//
double randomG::normal01()
{
	double U1, U2, V1, V2, W, Y, X1, X2;

	do {
		U1 = uniform01(); 
		U2 = uniform01();
		V1 = 2*U1 - 1;
		V2 = 2*U2 - 1;
		W  = V1*V1 + V2*V2;
	} while( W > 1 );

	Y = sqrt( (-2*log(W)) / W );
	X1 = V1 * Y;
	X2 = V2 * Y;

	return X1;
}

//
// generate a sample from a normal distribution with 
// mean 'mean' and variance 'variance'.
//
double randomG::normal( double mean, double variance )
{
	double X = normal01();
	double stddev = sqrt( variance );

	return mean + stddev * X;
}
