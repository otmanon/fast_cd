#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>

struct fast_cd_arap_sim 
{

public:

	//Geometry
	Eigen::MatrixXd V;
	//Tets
	Eigen::MatrixXi T;

	//Subspace
	Eigen::MatrixXd B;

	//labels
	Eigen::VectorXi labels;

	fast_cd_sim_params params;

	local_global_solver solver;

	fast_cd_static_precomp sp;



	void init(fast_cd_sim_params params, local_global_solver solver, fast_cd_static_precomp sp)
	{
		this->params = params;
		this->solver = solver;
		this->sp = sp;
	}


	void step(Eigen::VectorXd& z, Eigen::VectorXd& p, fast_cd_sim_state state, Eigen::VectorXd& f_ext, Eigen::VectorXd& bc)
	{
		/*fast_cd_arap_dynamic_precomp dp;
		dp.precomp(z, p, state, f_ext, bc);

		z = solver.solve(z, dpre, spre);*/
	}

};

struct fast_cd_sim_state 
{
	Eigen::VectorXd z_curr, z_prev, p_curr, p_prev;

	void init(Eigen::VectorXd& z_curr, Eigen::VectorXd& z_prev, Eigen::VectorXd& p_curr, Eigen::VectorXd& p_prev)
	{
		this->z_curr = z_curr;
		this->z_prev = z_prev;
		this->p_curr = p_curr;
		this->p_prev = p_prev;
	}

	void update(Eigen::VectorXd& z_next, Eigen::VectorXd& p_next)
	{
		this->z_prev = z_curr;
		this->z_curr = z_next;
		this->p_prev = p_curr;
		this->p_curr = p_next;
	}
};

struct fast_cd_sim_params
{
	//timestep
	double h;
	//youngs modulus (one per each tet)
	Eigen::VectorXd ym;
	//poisson's ratio
	Eigen::VectorXd pr;


	double invh2;

	int num_modes;
	int num_clusters;
	

};

struct local_global_solver_params
{

};

struct local_global_solver
{
	local_global_solver_params params;
};


struct fast_cd_static_precomp
{

};


struct fast_cd_arap_dynamic_precomp
{

};