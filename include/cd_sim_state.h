#pragma once
#include <Eigen/Core>
using namespace Eigen;
struct cd_sim_state
{
	VectorXd z_curr, z_prev, p_curr, p_prev;

	cd_sim_state() {};
	cd_sim_state(const VectorXd& z_curr, const VectorXd& z_prev,const  VectorXd& p_curr, const VectorXd& p_prev)
	{
		this->z_curr = z_curr;
		this->z_prev = z_prev;
		this->p_curr = p_curr;
		this->p_prev = p_prev;
	}

	
void init(const VectorXd& z_curr, const VectorXd& z_prev, const VectorXd& p_curr, const VectorXd& p_prev)
	{
		this->z_curr = z_curr;
		this->z_prev = z_prev;
		this->p_curr = p_curr;
		this->p_prev = p_prev;
	}

	void update(const VectorXd& z_next, const VectorXd& p_next)
	{
		this->z_prev = z_curr;
		this->z_curr = z_next;
		this->p_prev = p_curr;
		this->p_curr = p_next;
	}
};
