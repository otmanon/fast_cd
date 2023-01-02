#pragma once
#include "cd_sim_params.h"
#include "cd_sim_state.h"
#include "momentum_leaking_matrix.h"

#include <igl/massmatrix.h>
using namespace std;
/*
A helper class that aids in organizing different external force types. Construct this before your simulation,
then each simulation step query the external force with the .get() function
*/
struct cd_external_force
{


	// which external force are we using, so far, just "none" and "momentum_leak" is supported
	string external_force_type;

	double external_force_magnitude;

	// precomputation for "momentum_leak" force if active
	SparseMatrix<double> invh2MDJ;

	cd_external_force(cd_sim_params& sim_params, string external_force_type = "none", double external_force_magnitude = 0)
	{
		this->external_force_type = external_force_type;
		this->external_force_magnitude = external_force_magnitude;

		if (external_force_type == "momentum_leak")
			init_momentum_leak(sim_params);
	}

	virtual void init_momentum_leak(cd_sim_params& sim_params)
	{
		//calculate diagonal momentum diffusion matrix with default params
		SparseMatrix<double> D, M;
		momentum_leaking_matrix(sim_params.X, sim_params.T, fast_cd::MOMENTUM_LEAK_DIFFUSION, D);
		igl::massmatrix(sim_params.X, sim_params.T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
		invh2MDJ = sim_params.invh2 * (sim_params.J.transpose() *
			igl::repdiag(M, 3) * igl::repdiag(D, 3)).transpose();
	}

	/*
	Returns the external force being supplied to the fast complementary dynamics system.
	Inputs:
	step - which timestep of the simulation are we in. This is useful for forces that have a time-varying component
	sol_p -   12|B|x1 flattened rig parameters at next timestep
	state - fast_cd_state struct that contains info on z_curr, z_prev, p_curr and p_prev. Useful for inertial-like external forces
	*/
	virtual VectorXd get(int step, VectorXd& p, cd_sim_state& state)
	{
		if (external_force_type == "none")
		{
			return VectorXd::Zero(state.z_curr.rows());
		}
		else if (external_force_type == "momentum_leak")
		{
			assert(state.z_curr.rows() == invh2MDJ.rows() && "external force does not match state dimensinoality!");
			return external_force_magnitude * invh2MDJ * (2.0 * state.p_curr - state.p_prev - p);
		}
	}
};