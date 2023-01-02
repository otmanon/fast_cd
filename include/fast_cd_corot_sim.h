#pragma once
#include "fast_cd_corot_sim_params.h"
#include "fast_cd_corot_local_global_solver.h"
#include "fast_cd_corot_precomp.h"
struct fast_cd_corot_sim 
{

public:

	fast_cd_sim_params* params;

	fast_cd_corot_local_global_solver* sol;

	fast_cd_corot_static_precomp* sp;

	fast_cd_corot_dynamic_precomp* dp;

	fast_cd_corot_sim() {};

	fast_cd_corot_sim(fast_cd_sim_params& params, local_global_solver_params& solver_params) {
		this->params = &params;
		this->sp = new fast_cd_corot_static_precomp(*this->params);
		this->dp = new fast_cd_corot_dynamic_precomp();

		this->sol = new fast_cd_corot_local_global_solver(sp->BAB, this->params->AeqB, solver_params);

	};
	fast_cd_corot_sim(fast_cd_sim_params& params, local_global_solver_params& solver_params, fast_cd_corot_dynamic_precomp& dp, fast_cd_corot_static_precomp& sp)
	{
		this->params = &params;
		this->sp = &sp;
		this->dp = &dp;
		this->sol = new fast_cd_corot_local_global_solver(this->sp->BAB, this->params->AeqB, solver_params);
	}

	virtual Eigen::VectorXd step(const VectorXd& z, const VectorXd& p, const cd_sim_state& state, const  VectorXd& f_ext, const  VectorXd& bc)
	{
		//needs to be updated every timestep
		assert(params->Aeq.rows() == bc.rows() && "Need rhs of linear constraint to match lhs");
		assert(f_ext.rows() == z.rows() && "Force needs to be of same dimensionality as D.O.F.'s");

		dp->precomp(z, p, state, f_ext, bc, *sp);
		VectorXd z_next = sol->solve(z, p, *params, *dp, *sp);
		return z_next;
	}


	/*
	Updates the equality constraint matrix to be C, where C slices out known rows of the full space system
	*/
	virtual void set_equality_constraint(SparseMatrix<double>& C)
	{
		params->Aeq = C;

		sp = new fast_cd_corot_static_precomp(*params);
		local_global_solver_params old_params = sol->sol_p;

		sol = new fast_cd_corot_local_global_solver(sp->BAB, params->AeqB, old_params);
	}

	/*
	Projects z to full space
	*/
	virtual VectorXd full_space(VectorXd& z)
	{
		return z; //for a full space sim, z is already full space
	}


	fast_cd_sim_params* parameters()
	{
		return this->params;
	}

};