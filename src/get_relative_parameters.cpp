#include "get_relative_parameters.h"
#include "matrix4f_from_parameters.h"
#include "update_parameters_at_handle.h"
/*
From rest pose parameters p0, and deformed pose parameters p, get p_rel
*/
void get_relative_parameters(Eigen::VectorXd& p0, Eigen::VectorXd& p, Eigen::VectorXd& p_rel)
{
	int num_b = p0.rows() / 12;
	Eigen::Matrix4f A0, A_rel, A;
	p_rel.resizeLike(p0);
	for (int i = 0; i < num_b; i++)
	{
		A0 = matrix4f_from_parameters(p0, i);
		A = matrix4f_from_parameters(p, i);
		A_rel = A * A0.inverse();
		
		update_parameters_at_handle(p_rel, A_rel, i);
	}


}