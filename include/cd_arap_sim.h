#pragma once
#include "cd_sim_params.h"
#include "cd_sim_state.h"
#include "cd_arap_local_global_solver.h"
#include "cd_arap_precomp.h"



#include <Eigen/Core>
#include <Eigen/Sparse>

using namespace std; using namespace Eigen;



struct cd_arap_sim
{
public:

	cd_sim_params* params;

	cd_arap_local_global_solver* sol;

	cd_arap_static_precomp* sp;

	cd_arap_dynamic_precomp* dp;

	cd_arap_sim() {};
	cd_arap_sim(cd_sim_params* params, cd_arap_local_global_solver* solver, cd_arap_dynamic_precomp* dp, cd_arap_static_precomp* sp)
	{
		this->params = params;
		this->sol =  solver;
		this->sp = sp;
		this->dp = dp;
	}

	Eigen::VectorXd step(VectorXd& z, VectorXd& p, cd_sim_state state, VectorXd& f_ext, VectorXd& bc)
	{
		//needs to be updated every timestep
		assert(params->Aeq.rows() == bc.rows() && "Need rhs of linear constraint to match lhs");
		assert(f_ext.rows() == z.rows() && "Force needs to be of same dimensionality as D.O.F.'s");

		dp->precomp(z, p, state, f_ext, bc, *sp);
		z = sol->solve(z, *params, *dp, *sp);
		return z;
	}
};