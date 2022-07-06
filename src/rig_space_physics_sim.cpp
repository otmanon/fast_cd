#include "rig_space_physics_sim.h"
#include <igl/boundary_facets.h>
#include <igl/unique.h>
#include <igl/colon.h>
RigPhysicsSim::RigPhysicsSim(Eigen::MatrixXd& X, Eigen::MatrixXi& T, Eigen::SparseMatrix<double>& J_user, Eigen::SparseMatrix<double>& J_physics, double ym, double pr, double dt, int num_clusters, std::string clusters_dir, bool do_clustering)
{
	B_full = J_physics;
	B = J_physics; //maybe need a transpose

	L_full.resize(B_full.cols());
	L_full.setOnes();
	L = L_full;
	J = J_user;
	
	this->X = X;
	this->T = T;
	this->J = J;
	this->num_modes = num_modes;
	this->num_clusters = num_clusters;
	this->do_clustering = do_clustering;
	this->num_modal_features = J_physics.cols();
	this->do_reduction = true; //this is always true for this simulation
	this->dt = dt;

	this->clusters_file_dir = clusters_dir;

	stiffness = ym / (2.0 * (1.0 + pr));
	incompressibility = ym * pr / ((1.0 + pr) * (1.0 - 2.0 * pr));

	constraints.use_constraints = false;
	Eigen::VectorXi _n; //throwaway var
	igl::boundary_facets(T, F, FiT, _n);
	igl::unique(F, ext_ind);
	this->num_modal_features = num_modal_features;


	reduced_newton_solver = new DenseQuasiNewtonSolver();
	full_newton_solver = new QuasiNewtonSolver();

	constraints.reduced_newton_solver = new DenseQuasiNewtonSolver();
	constraints.full_newton_solver = new QuasiNewtonSolver();

	if (do_clustering)
		init_clusters(num_clusters, num_modal_features);
	else
		igl::colon(0, T.rows() - 1, labels);
	this->num_clusters = labels.maxCoeff() + 1;
	init_system_matrices();
	precompute_solvers();

}




Eigen::VectorXd RigPhysicsSim::step(Eigen::VectorXd& p_user_next, Eigen::VectorXd& p_user_curr, Eigen::VectorXd& p_user_prev, Eigen::VectorXd& p_sim_curr, Eigen::VectorXd& p_sim_prev)
{
	return reduced_step(p_user_next, p_user_curr, p_user_prev, p_sim_curr, p_sim_prev, false, 20);  //same as reduced step but B space (allowable motion) is rig jacobian
}