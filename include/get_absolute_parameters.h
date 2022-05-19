#pragma once
#include <Eigen/Core>
#include "matrix4f_from_parameters.h"
#include "update_parameters_at_handle.h"

/*
Given rest pose params and relative params, computes absolute params
*/
void get_absolute_parameters(Eigen::VectorXd& p0, Eigen::VectorXd& p_rel, Eigen::VectorXd& p)
{
	int num_b = p0.rows() / 12;
	Eigen::Matrix4f A0, A_rel, A;
	p.resizeLike(p0);
	for (int i = 0; i < num_b; i++)
	{
		A0 = matrix4f_from_parameters(p0, i);
		A_rel = matrix4f_from_parameters(p_rel, i);
		A = A_rel * A0;
		update_parameters_at_handle(p, A, i);
	}

}
