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
		Eigen::VectorXd weight;

    public:
    	Group *group;

		int mu;
		int lambda;
		double sigma;
		CMAES(int mu_ref, int lambda_ref, double sigma_ref, Group *group_ref);
		void init_static_parameters();
		~CMAES();

		Eigen::MatrixXd B, D;
		Eigen::MatrixXd covar;
		Eigen::VectorXd pc , ps;
		// this can be obtained by (mean_new - mean_old) / sigma
		// actually don't need to be remembered during iteration

		void value_update(Eigen::VectorXd z_mean, Eigen::MatrixXd zi);
		int h_sig();
		void update_ps(Eigen::VectorXd z_mean);
		void update_pc(Eigen::VectorXd z_mean);
		void update_covar(Eigen::MatrixXd zi);
		void update_sigma();
		void run();
		Eigen::VectorXd MVNsample();
};

