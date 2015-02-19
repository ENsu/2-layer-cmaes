#include "cmaes.h"
#include <algorithm>
#include <math.h> 
#include "Eigen/Dense"
#include "global.h"
#include "random.hpp" 

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

extern randomG RANDOM;

CMAES::CMAES(int mu_ref, int lambda_ref, double sigma_ref, Group *group_ref)
{
	group = group_ref;
	mu = mu_ref;
	sigma = sigma_ref;
	lambda = lambda_ref;

	//-------init dynamic parameters---------
    pc.setZero(dimension);
    ps.setZero(dimension);
    covar.setIdentity(dimension , dimension);
    B.setIdentity(dimension , dimension);
    D.setIdentity(dimension , dimension);
    init_static_parameters();
}

void CMAES::init_static_parameters()  // these values won't be changed during cma-es
{
	double weight_sum_2 = 0.0;
    double weight_sum = 0.0;
    weight.setZero(mu);

    for(int i = 0 ; i < mu ; i++)
    {
		weight(i) = log(mu+0.5) - log(i+1);
		weight_sum += weight(i);
    }
    for(int i = 0 ; i < mu ; i++)
    {
		weight(i) /= weight_sum;
		weight_sum_2 += (weight(i) * weight(i)) ;
    }
    mu_w = 1 / weight_sum_2;

    weight_sum = 0.0;
    for(int i=0; i<mu; i++)
    {
    	weight_sum += weight(i);
    }
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

CMAES::~CMAES(){}


void CMAES::run() // run an cma-es iteration
{
	//-----------step 1: sample new lambda - mu points---------
	int new_sample_num = lambda - mu;
	assert(new_sample_num > 0);
	Eigen::MatrixXd old_xi = group->node_matrix();
	Node *y = new Node[new_sample_num];
	Node mean = group->get_mean_node(weight);
	// assert(old_xi * weight == group->get_mean_node(weight).allele);
	for(int i=0; i<new_sample_num; i++) //prepare the new (lambda - mu) sample points
	{
		Eigen::VectorXd sample = MVNsample();

		sample = (sample * sigma) + mean.allele;
		Node tmp(sample);
		while(tmp.outofBound()) // just resample when node out of boundary
		{
			sample = MVNsample();
			sample = (sample * sigma) + mean.allele;
			tmp = Node(sample);
		}
		y[i] = tmp;
	}

	//-----------step 2: add these node in group and do selection---------
	for(int i=0; i<new_sample_num; i++)
	{
		group->add_node(y[i]);
	}
	group->sort_node();
	group->truncate_size(mu);

	//-----------step 3: update relative variables, pc, ps, C, sigma---------
	//mean is updated in Group class
	Eigen::MatrixXd new_xi = group->node_matrix();
	Eigen::MatrixXd zi = D.inverse() * B.transpose() * (new_xi - old_xi) / sigma;
	Eigen::VectorXd z_mean = zi * weight;
	value_update(z_mean, zi); // update the relativa value, 
	return ;
}

// update pc, ps, covar, sigma
void CMAES::value_update(Eigen::VectorXd z_mean, Eigen::MatrixXd zi)
{
	update_ps(z_mean);
	update_pc(z_mean);
	update_covar(zi);
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

void CMAES::update_ps(Eigen::VectorXd z_mean)
{
    ps = (1.0-cs)*ps + sqrt(cs*(2-cs)*mu_w) * B * z_mean;
    assert(ps == ps);
    return ;
}

void CMAES::update_pc(Eigen::VectorXd z_mean)
{
	int hsig = h_sig();
    pc = (1.0-cc)*pc + hsig * sqrt(cc*(2-cc)*mu_w) * B * D * z_mean;
    assert(pc == pc);
    return ;
}

void CMAES::update_covar(Eigen::MatrixXd zi)
{
	int hsig = h_sig();
	double delta_hsig = (1 - hsig) * cc * (2 - cc);
	Eigen::MatrixXd yi = B * D * zi;
	MatrixXd diag_weight = MatrixXd(weight.asDiagonal());
    covar = (1-c1-cmu) * covar + c1*(pc*pc.transpose() + delta_hsig*covar) + cmu * yi * diag_weight * yi.transpose();
    assert(covar == covar);

    // force covar to be symmetric and positive!!
    covar = covar.cwiseAbs(); // this line is not used in the demo matlab code!!!
   	Eigen::MatrixXd Covs = Eigen::MatrixXd(covar.triangularView<Eigen::StrictlyUpper>());
	Eigen::MatrixXd Covd = Eigen::MatrixXd(covar.triangularView<Eigen::Upper>());
	covar = Covd + Covs.transpose();

	//update B & D
	cout << "==========covar========" << endl << covar << endl;
	Eigen::EigenSolver<MatrixXd> es(covar);
	assert(es.eigenvalues().imag().isZero());
	D = es.eigenvalues().real().cwiseSqrt().asDiagonal();
	assert(es.eigenvectors().imag().isZero());
	B = es.eigenvectors().real();
	cout << "==========B: ==========" << endl << B << endl;
	for(int i=0; i<dimension; i++)
	{
		cout << "i = " << i << endl;
		for(int j=i+1; j<dimension; j++)
		{
			cout << "j = " << j << endl;
			cout << "dot: " << B.col(i).dot(B.col(j)) << endl; 
			assert(B.col(i).dot(B.col(j)) < 1E-10); 
		}
		cout << "norm: " << B.col(i).norm() << endl;
		assert(fabs(B.col(i).norm() - 1) < 1E-10);
	}
	cout << "B * BT = " << (B * B.transpose()).transpose() << endl; 
	cout << "BT * B = " << B.transpose() * B << endl;
	Eigen::MatrixXd tmp;
	tmp.setIdentity(dimension, dimension);
	assert(tmp.isApprox(B * B.transpose()));
	assert(tmp.isApprox(B.transpose() * B));
	assert(covar.isApprox(B*D*D*B.transpose()));
//	cout << "check0: " << endl << B.inverse() << endl << "-----" << endl << B.transpose() << endl;
//	cout << "check1: " << endl << B * D * D * B.transpose() << endl;
    return ;
}

void CMAES::update_sigma()
{
    double val = exp ( cs / ds	* ((ps.norm() / e_n01 )- 1 ) );
    sigma = sigma * val;
    assert(sigma == sigma);
    return ;
}

// multivariate noraml sampling with certain covar
Eigen::VectorXd CMAES::MVNsample()
{
	Eigen::VectorXd sample(dimension);
	for(int i=0; i<dimension; i++)
		sample(i) = RANDOM.normal01();
	return B * D * sample;
} // ref from http://stackoverflow.com/questions/6142576/sample-from-multivariate-normal-gaussian-distribution-in-c
