#pragma once
#include <Eigen/Core>
using namespace Eigen;
struct cd_sim_state
{
	VectorXd z_curr, z_prev, p_curr, p_prev;

	void init(VectorXd& z_curr, VectorXd& z_prev, VectorXd& p_curr, VectorXd& p_prev)
	{
		this->z_curr = z_curr;
		this->z_prev = z_prev;
		this->p_curr = p_curr;
		this->p_prev = p_prev;
	}

	void update(VectorXd& z_next, VectorXd& p_next)
	{
		this->z_prev = z_curr;
		this->z_curr = z_next;
		this->p_prev = p_curr;
		this->p_curr = p_next;
	}
};
