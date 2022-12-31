#pragma once
using namespace Eigen;
struct local_global_solver_params
{
	int max_iters;
	bool to_convergence;
	double threshold;

	local_global_solver_params() {};

	local_global_solver_params(bool to_convergence, int max_iters, double threshold)
	{
		this->to_convergence = to_convergence;
		this->max_iters = max_iters;
		this->threshold = threshold;
	}
};
