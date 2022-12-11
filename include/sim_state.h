#pragma once
#include "cd_sim_state.h"
#include <Eigen/Core>

using namespace Eigen;
struct sim_state :cd_sim_state
{

	sim_state(){};
	
	sim_state(const VectorXd& z_curr, const VectorXd& z_prev)
	{
		this->z_curr = z_curr;
		this->z_prev = z_prev;
		this->p_curr = VectorXd(0);
		this->p_prev = VectorXd(0);
	};

	void init(const VectorXd& z_curr, const VectorXd& z_prev)
	{
		this->z_curr = z_curr;
		this->z_prev = z_prev;
		this->p_curr = VectorXd(0);
		this->p_prev = VectorXd(0);
	}


	void update(const VectorXd& z_next)
	{
		this->z_prev = z_curr;
		this->z_curr = z_next;
	}
};