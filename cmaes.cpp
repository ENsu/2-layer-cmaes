#include "cmaes.h"
#include <algorithm>
#include <math.h> 
#include "Eigen/Dense"
#include "global.h"
#include "random.hpp" 
#include <cfloat>

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;

extern randomG RANDOM;

#define DEBUG false

CMAES::CMAES(){}

CMAES::CMAES(int mu_ref, int lambda_ref, double sigma_ref, Group *group_ref)
{
	group = group_ref;
	mu = mu_ref;
	sigma = sigma_ref;
	lambda = lambda_ref;

	termination = false;
	termination_count = 0;
	best_fitness_val = group->getMin();

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
    	//cout << "weight["<<i<<"]: " << weight(i) << endl; 
    }
    assert(weight_sum-1 < 1E-10);
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

list<Node> CMAES::sample_node(int node_num)
{
	list<Node> new_sample_node;
	Eigen::VectorXd mean = group->get_mean_node(weight).allele;
	if(DEBUG)
	{
		cout << "-----sampling condition:-----" << endl;
		cout << "mean: " << mean << endl;
		cout << "outofBound?? " << (mean(0) - 3.14159265359) << endl;
		cout << "D: " << D << endl;
		cout << "sigma: " << sigma << endl;
	}
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

	//----------step 4: check best fitness improvement------------
	if(best_fitness_val > group->getMin())
	{
		best_fitness_val = group->getMin();
		termination_count = 0;
	}
	else
	{
		termination_count++;
		if(termination_count >= (10 + 30 * dimension/lambda))
			termination=true;
	}

	return ;
}

// cmaes can passively receive a new group as sampling result
// and update its related parameters, B, D, covar, pc, ps and sigma
void CMAES::update_value(Group new_group)
{
	Eigen::MatrixXd old_xi = group->node_matrix();
	assert((old_xi * weight).isApprox(group->get_mean_node(weight).allele));
	Eigen::MatrixXd new_xi = new_group.node_matrix();
	assert((new_xi * weight).isApprox(new_group.get_mean_node(weight).allele));
	Eigen::MatrixXd zi = D.inverse() * B.transpose() * (new_xi - old_xi) / sigma;
	Eigen::VectorXd z_mean = zi * weight;
	
	double diver_check_old = sigma * D.maxCoeff();
	if(DEBUG)
	{
		cout << "D.inverse " << endl << D.inverse() << endl;
		cout << "B.transpose " << endl << B.transpose() << endl;
		cout << "(new_xi - old_xi) " << endl << (new_xi - old_xi) << endl;
		cout << "sigma: " << sigma << endl;
		cout << "weight: " << endl << weight << endl;
		cout << "z_mean: " << endl << z_mean << endl;
	}

	update_ps(z_mean);
	update_pc(z_mean);
	update_covar(zi);
	update_sigma();

	//-----------------termination criterion-------------
	double diver_check_new = sigma * D.maxCoeff();
	if((diver_check_new/diver_check_old) > 1E4)
		termination = true;
	//----------------------------------------------------
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
	if(DEBUG)
		cout << "pre ps: " << ps << endl;
    ps = (1.0-cs)*ps + sqrt(cs*(2-cs)*mu_w) * B * z_mean;
    if(DEBUG)
    	cout << "post ps: " << ps << endl;
    assert(ps == ps);
    return ;
}

void CMAES::update_pc(Eigen::VectorXd z_mean)
{
	int hsig = h_sig();
	if(DEBUG)
		cout << "pre pc: " << pc << endl;
    pc = (1.0-cc)*pc + hsig * sqrt(cc*(2-cc)*mu_w) * B * D * z_mean;
	if(DEBUG)
		cout << "post pc: " << pc << endl;
    assert(pc == pc);
    return ;
}

void CMAES::update_covar(Eigen::MatrixXd zi)
{
	if(DEBUG)
	{
		cout << "pre B: " << endl << B << endl;
		cout << "pre D: " << endl << D << endl;
		cout << "pre covar: " << endl << covar << endl; 
	}
	int hsig = h_sig();
	double delta_hsig = (1 - hsig) * cc * (2 - cc);
	Eigen::MatrixXd yi = B * D * zi;
	MatrixXd diag_weight = MatrixXd(weight.asDiagonal());
    covar = (1-c1-cmu) * covar + c1*(pc*pc.transpose() + delta_hsig*covar) + cmu * yi * diag_weight * yi.transpose();
    assert(covar == covar);

    // force covar to be symmetric and positive!!
    // covar = covar.cwiseAbs(); // this line is not used in the demo matlab code!!!
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
	if(DEBUG)
	{
		cout << "post B: " << endl << B << endl;
		cout << "post D: " << endl << D << endl;
		cout << "post covar: " << endl << covar << endl; 
	}
	for(int i=0; i<dimension; i++)
	{
		//assert((covar * B.col(i) - eigenvalues(i) * B.col(i)).isMuchSmallerThan(eigenvalues(i) * B.col(i)));
		assert(fabs(B.col(i).norm() - 1) < 1E-10);
	}
	assert(covar.isApprox(B*D*D*B.inverse()));
    return ;
}

void CMAES::update_sigma()
{
	if(DEBUG)
    	cout << "pre sigma: " << sigma << endl;
    double val = exp ( cs / ds	* ((ps.norm() / e_n01 )- 1 ) );
    sigma = sigma * val;
    if(DEBUG)
    	cout << "post sigma: " << sigma << endl;
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