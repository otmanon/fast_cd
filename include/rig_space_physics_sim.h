#pragma once
#include "fast_complementary_dynamics_sim.h"

class RigPhysicsSim : public FastCDSim
{

public: 
	RigPhysicsSim(Eigen::MatrixXd& X, Eigen::MatrixXi& T, Eigen::SparseMatrix<double>& J_user, Eigen::SparseMatrix<double>& J_physics, double ym, double pr, double dt, int num_clusters, std::string clusters_dir, bool do_clustering);

	Eigen::VectorXd step(Eigen::VectorXd& p_user_next, Eigen::VectorXd& p_user_curr, Eigen::VectorXd& p_user_prev, Eigen::VectorXd& p_sim_curr, Eigen::VectorXd& p_sim_prev);


	Eigen::VectorXd parameters_to_displacements(Eigen::VectorXd& p_user, Eigen::VectorXd& p_sim);
};