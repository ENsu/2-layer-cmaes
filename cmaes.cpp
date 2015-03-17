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

CMAES::CMAES(){}

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
}

CMAES::~CMAES(){}

// if its Group's Node number not equal to mu
// cmaes changes its number of Node to be equal to mu
void CMAES::tune_node_num()
{
	// if the number of nodes in group is different from mu
	// we add/substract nodes to meet mu
	if(group->Nodes.size() > mu)
	{
		group->sort_node();
		group->truncate_size(mu);
	}
	else if(group->Nodes.size() < mu)
	{
		int node_added = mu - group->Nodes.size();
		Eigen::VectorXd mean = group->get_mean_node().allele;
		for(int i=0; i < node_added; i++)
		{
			Node tmp(MVNsample(mean));
			while(tmp.outofBound()) // just resample when node out of boundary
			{
				tmp = Node(MVNsample(mean));
			}
			group->add_node(tmp);
		}
	}
	assert(group->Nodes.size() == mu);
}

void CMAES::run() // run an cma-es iteration
{
	//check group's node number == mu value
	//if not, you shoule run function tune_node_num first
	assert(group->Nodes.size() == mu);
	//-----------step 1: sample new lambda - mu points and add to group---------
	int new_sample_num = lambda - mu;
	list<Node> new_sample_node = sample_node(new_sample_num);

	//--------step 2: create new group with selected node---------
	list<Node> all_node_list = group->Nodes;
	all_node_list.splice(all_node_list.end(), new_sample_node);
	Group *new_group = new Group(all_node_list);

	new_group->sort_node();
	new_group->truncate_size(mu);
	//-----------step 3: update relative variables, pc, ps, C, sigma---------
	update_value(*new_group);
	group = new_group;
	return ;
}

list<Node> CMAES::sample_node(int node_num)
{
	list<Node> new_sample_node;
	Eigen::VectorXd mean = group->get_mean_node(weight).allele;
	for(int i=0; i<node_num; i++) //prepare the new (lambda - mu) sample points
	{	
		Node tmp(MVNsample(mean));
		while(tmp.outofBound()) // just resample when node out of boundary
		{
			tmp = Node(MVNsample(mean));
		}
		new_sample_node.push_back(tmp);
	}
	return new_sample_node;
}

// cmaes can passively receive a new group as sampling result
// and update its related parameters, B, D, covar, pc, ps and sigma
void CMAES::update_value(Group new_group)
{
	Eigen::MatrixXd old_xi = group->node_matrix();
	assert(old_xi * weight == group->get_mean_node(weight).allele);
	Eigen::MatrixXd new_xi = new_group.node_matrix();
	assert(new_xi * weight == new_group.get_mean_node(weight).allele);
	group = &new_group;
	Eigen::MatrixXd zi = D.inverse() * B.transpose() * (new_xi - old_xi) / sigma;
	Eigen::VectorXd z_mean = zi * weight;

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
	Eigen::EigenSolver<MatrixXd> es(covar);
	assert(es.eigenvalues().imag().isZero());
	Eigen::VectorXd eigenvalues = es.eigenvalues().real();
	D = es.eigenvalues().real().cwiseSqrt().asDiagonal();
	assert(es.eigenvectors().imag().isZero());
	B = es.eigenvectors().real();
	for(int i=0; i<dimension; i++)
	{
		assert((covar * B.col(i)).isApprox(eigenvalues(i) * B.col(i)));
		assert(fabs(B.col(i).norm() - 1) < 1E-10);
	}
	assert(covar.isApprox(B*D*D*B.inverse()));
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
Eigen::VectorXd CMAES::MVNsample(VectorXd mean)
{
	Eigen::VectorXd sample(dimension);
	for(int i=0; i<dimension; i++)
		sample(i) = RANDOM.normal01();
	return mean + sigma * B * D * sample;
} // ref from http://stackoverflow.com/questions/6142576/sample-from-multivariate-normal-gaussian-distribution-in-c

CMAES& CMAES::operator=(const CMAES rhs)
{
    if(this == &rhs)
        return *this;

    cc = rhs.cc;
    cs = rhs.cs;
    c1 = rhs.c1;
    cmu = rhs.cmu;
    ds = rhs.ds;
    e_n01 = rhs.e_n01;
    mu_w = rhs.mu_w;

    group = rhs.group; // they share the same group, right?
    mu = rhs.mu;
    lambda = rhs.lambda;
    sigma = rhs.sigma;
    weight = rhs.weight;
    B = rhs.B;
    D = rhs.D;
    covar = rhs.covar;
    pc = rhs.pc;
    ps = rhs.ps;

    return *this;
}