#pragma once
#include "sim_params.h"
#include "sim_state.h"
#include "corot_local_global_solver.h"
#include "corot_precomp.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

using namespace std; using namespace Eigen;



struct corot_sim
{
public:

	corot_sim_params* params;

	corot_local_global_solver* sol;

	corot_static_precomp* sp;

	corot_dynamic_precomp* dp;

	corot_sim() {};
	corot_sim(corot_sim_params& params, local_global_solver_params& solver_params) {
		this->params = &params;
		this->sp = new corot_static_precomp(*this->params);
		this->dp = new corot_dynamic_precomp();

		this->sol = new corot_local_global_solver(sp->A, this->params->Aeq, solver_params);

	};
	corot_sim(corot_sim_params& params, local_global_solver_params& solver_params, corot_dynamic_precomp& dp, corot_static_precomp& sp)
	{
		this->params = &params;
		this->sp = &sp;
		this->dp = &dp;
		this->sol = new corot_local_global_solver(sp.A, params.Aeq, solver_params);
	}

	virtual Eigen::VectorXd step(const VectorXd& z, const sim_state& state, const  VectorXd& f_ext, const  VectorXd& bc)
	{
		//needs to be updated every timestep
		assert(params->Aeq.rows() == bc.rows() && "Need rhs of linear constraint to match lhs");
		assert(f_ext.rows() == z.rows() && "Force needs to be of same dimensionality as D.O.F.'s");

		dp->precomp(z,  state, f_ext, bc, *sp);
		VectorXd z_next = sol->solve(z, *params, *dp, *sp);
		return z_next;
	}


	/*
	Updates the equality constraint matrix to be C, where C slices out known rows of the full space system
	*/
	virtual void set_equality_constraint(SparseMatrix<double>& C)
	{
		params->Aeq = C;

		sp = new corot_static_precomp(*params);
		local_global_solver_params old_params = sol->p;

		sol = new corot_local_global_solver(sp->A, params->Aeq, old_params);
	}

	/*
	Projects z to full space
	*/
	virtual VectorXd full_space(VectorXd& z)
	{
		return z; //for a full space sim, z is already full space
	}


	corot_sim_params* parameters()
	{
		return this->params;
	}

};