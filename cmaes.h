#include "Eigen/Dense"
#include "group.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;

class CMAES
{
	private:
		double cc , cs , c1 , cmu , ds;
		double e_n01;
		double mu_w;

    public:
    	Group *group;

		int mu;
		int lambda;
		double sigma;
		bool termination;
		CMAES();
		CMAES(int mu_ref, int lambda_ref, double sigma_ref, Group *group_ref);
		void init_static_parameters();
		~CMAES();

		Eigen::VectorXd weight;
		Eigen::MatrixXd B, D;
		Eigen::MatrixXd covar;
		Eigen::VectorXd pc , ps;
		// this can be obtained by (mean_new - mean_old) / sigma
		// actually don't need to be remembered during iteration

		void tune_node_num();
		void run();
		list<Node> sample_node(int node_num);
		void update_value(Group new_group);
		int h_sig();
		void update_ps(Eigen::VectorXd z_mean);
		void update_pc(Eigen::VectorXd z_mean);
		void update_covar(Eigen::MatrixXd zi);
		void update_sigma();

		Eigen::VectorXd MVNsample(VectorXd mean);

		CMAES& operator=(const CMAES rhs);
};

