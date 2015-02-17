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
		double *weight;

    public:
    	Group *group;

		int mu;
		int lambda;
		double sigma;
		CMAES(int mu_ref, int lambda_ref, double sigma_ref, Group *group_ref);
		void init_unchangable_values();
		~CMAES();

		Eigen::MatrixXd covar;
		Eigen::VectorXd pc , ps;
		Eigen::VectorXd yw; 
		// this can be obtained by (mean_new - mean_old) / sigma
		// actually don't need to be remembered during iteration

		void value_update(Eigen::VectorXd new_mean, Eigen::VectorXd old_mean);
		int h_sig();
		void update_yw(Eigen::VectorXd new_mean, Eigen::VectorXd old_mean);
		void update_pc();
		void update_ps();
		void update_covar();
		void update_sigma();
		void run();
		Eigen::MatrixXd MVNsample(Eigen::MatrixXd covar);
};

