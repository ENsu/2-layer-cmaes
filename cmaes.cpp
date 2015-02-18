#include "cmaes.h"
#include <algorithm>
#include <math.h> 
#include "Eigen/Dense"
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/normal_distribution.hpp"
#include "global.h" 

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

CMAES::CMAES(int mu_ref, int lambda_ref, double sigma_ref, Group *group_ref)
{
	group = group_ref;
	mu = mu_ref;
	sigma = sigma_ref;
	lambda = lambda_ref;

	//-------init dynamic parameters---------
    pc.setZero(dimension);
    ps.setZero(dimension);
    yw.setZero(dimension);
    covar.setIdentity(dimension , dimension);
    init_static_parameters();
	// initialize covar with group's node's distribution
    //covar = group->getCov(weight, mu);
    //covar = covar / covar.norm();
}

void CMAES::init_static_parameters()  // these values won't be changed during cma-es
{
	double weight_sum_2 = 0;
    weight = new double[mu];

    double weight_sum = 0.0;
    for(int i = 0 ; i < mu ; i++)
    {
		weight[i] = log(mu+0.5) - log(i+1);
		weight_sum += weight[i];
    }
    for(int i = 0 ; i < mu ; i++)
    {
		weight[i] /= weight_sum;
		weight_sum_2 += ( weight[i] * weight[i]) ;
    }
    mu_w = 1 / weight_sum_2;

    int n = dimension;
	cs = (mu_w + 2) / (n + mu_w + 5);
	double ds_tmp = sqrt((mu_w - 1)/(n + 1)) - 1;
	ds = 1 + 2 * max(0.0, ds_tmp) + cs;
    cc = (4.0 + mu_w/n) / (n + 4 + 2*mu_w/n);
	c1 = 2 / ((n+1.3) * (n+1.3) + mu_w);
	int apha_mu = 2;
	double cmu_tmp = apha_mu*(mu_w - 2 + 1/mu_w) / ((n + 2)*(n + 2) + apha_mu*mu_w/2);
	cmu = min(1 - c1, cmu_tmp);
	e_n01 = sqrt(dimension) * ( 1.0 - 1.0 / (4*dimension) + 1.0 / (21 * dimension * dimension) );


	cout << "----------init static parameters----------" << endl;
	cout << "mu_w :" << mu_w << endl;
	cout << "cc: " << cc << endl;
	cout << "cs: " << cs << endl;
	cout << "c1: " << c1 << endl;
	cout << "cmu: " << cmu << endl;
	cout << "ds:" << ds << endl;
	cout << "----------------------------------------------" << endl;
}

CMAES::~CMAES()
{
	delete[] weight;
}


void CMAES::run() // run an cma-es iteration
{
	//-----------step 1: sample new lambda - mu points---------
	int new_sample_num = lambda - mu;
	assert(new_sample_num > 0);
	Eigen::MatrixXd old_mean = group->get_mean_node().allele;
	Node *y = new Node[new_sample_num];

	for(int i=0; i<new_sample_num; i++) //prepare the new (lambda - mu) sample points
	{
		Eigen::MatrixXd sample = MVNsample(covar).col(0);
		sample = (sample * sigma) + old_mean;
		Node tmp(sample);
		while(tmp.outofBound())
		{
			sample = MVNsample(covar).col(0);
			sample = (sample * sigma) + old_mean;
			tmp = Node(sample);
		}
		y[i] = tmp;
		cout << "--------new node--------";
		tmp.print();

	}
	group->print();
	//-----------step 2: add these node in group and do selection---------
	for(int i=0; i<new_sample_num; i++)
	{
		group->add_node(y[i]);
	}
	group->sort_node();
	group->truncate_size(mu);

	//-----------step 3: update relative variables, pc, ps, C, sigma---------
	//mean is updated in Group class
	Eigen::MatrixXd new_mean = group->get_mean_node().allele;
	value_update(new_mean, old_mean); // update the relativa value, 
	return ;
}

// update pc, ps, covar, sigma
void CMAES::value_update(Eigen::VectorXd new_mean, Eigen::VectorXd old_mean)
{
	update_yw(new_mean, old_mean);
	update_pc();
	update_ps();
	update_covar();
	update_sigma();
}

int CMAES::h_sig()
{
	double first_term = ps.norm() / sqrt(1 - pow((1-cs), 2*(generation+1)));
	double second_term = (1.4 + 2/(dimension+1)) * e_n01;
    if(first_term < second_term)
		return 1;
	else
		return 0;
}

void CMAES::update_yw(Eigen::VectorXd new_mean, Eigen::VectorXd old_mean)
{
	yw = (new_mean - old_mean) / sigma;
	assert(yw == yw);
	return ;
}

void CMAES::update_pc()
{
	int hsig = h_sig();
    pc = (1.0-cc)*pc + hsig * sqrt(cc*(2-cc)*mu_w) * yw;
    assert(pc == pc);
    return ;
}

void CMAES::update_ps()
{
	cout << covar << endl;
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(covar);
    Eigen::MatrixXd tmp = es.operatorInverseSqrt();
    ps = (1.0-cs)*ps + sqrt(cs*(2-cs)*mu_w) * tmp * yw;
    assert(ps == ps);
    return ;

}

void CMAES::update_covar()
{
	Eigen::MatrixXd groupCov = group->getCov(weight, mu);
	int hsig = h_sig();
	double delta_hsig = (1 - hsig) * cc * (2 - cc);
    covar = (1-c1-cmu) * covar + c1*(pc*pc.transpose() + delta_hsig*covar) + cmu * groupCov;
    //covar = covar / covar.norm();
    assert(covar == covar);
    return ;
}

void CMAES::update_sigma()
{
    double val = exp ( cs / ds	* ((ps.norm() / e_n01 )- 1 ) );
    sigma = sigma * val;
    assert(sigma == sigma);
    return ;
}

namespace Eigen {  
namespace internal {
template<typename Scalar>

struct scalar_normal_dist_op 
{
static boost::mt19937 rng;    // The uniform pseudo-random algorithm
mutable boost::normal_distribution<Scalar> norm;  // The gaussian combinator

EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

    template<typename Index>
    inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
};

template<typename Scalar> boost::mt19937 scalar_normal_dist_op<Scalar>::rng;
template<typename Scalar>
struct functor_traits<scalar_normal_dist_op<Scalar> >
{ enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false }; };

} // end namespace internal
} // end namespace Eigen
//ref from http://stackoverflow.com/questions/6142576/sample-from-multivariate-normal-gaussian-distribution-in-c


// multivariate noraml sampling with certain covar
Eigen::MatrixXd CMAES::MVNsample(Eigen::MatrixXd covar)
{
    int size = dimension; // Dimensionality (rows)
    int nn = 1;

    Eigen::internal::scalar_normal_dist_op<double> randN; // Gaussian functor
    //Eigen::internal::scalar_normal_dist_op<double>::rng.seed(ranseed); // Seed the rng
    Eigen::MatrixXd normTransform(size,size);
	Eigen::VectorXd zero;
	zero.setZero(dimension);

    Eigen::LLT<Eigen::MatrixXd> cholSolver(covar);

    if (cholSolver.info()==Eigen::Success) 
    {
		normTransform = cholSolver.matrixL();
    }
    else
    {
		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
		normTransform = eigenSolver.eigenvectors() 
		    * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }

    Eigen::MatrixXd samples = (normTransform 
	    * Eigen::MatrixXd::NullaryExpr(size,nn,randN)).colwise()
    + zero;

    return samples;
} // ref from http://stackoverflow.com/questions/6142576/sample-from-multivariate-normal-gaussian-distribution-in-c
