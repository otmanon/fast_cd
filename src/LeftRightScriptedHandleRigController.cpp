#include "LeftRightScriptedHandleRigController.h"
#include "update_parameters_at_handle.h"

void LeftRightScriptedHandleRigController::set_scripted_motion(int step)
{
	int num_handles = p_rel.rows() / 12;
	float  x = 0.5*sin(0.1 * step);

	Eigen::Matrix4f A_rel, A_rest, A;
	// moves every handle left-and right according to a sinusoid, starting at 0 motion and after a period of 2pi
	for (int i = 0; i < num_handles; i++)
	{
		A_rest = matrix4f_from_parameters(p_rest_inv, i);
		A_rel = matrix4f_from_parameters(p_rel, i);
		A_rel(0, 3) = x;
		update_parameters_at_handle(p_rel, A_rel, i);
		A = A_rel * A_rest;
		V.row(i) = A.block(0, 3, 3, 1).transpose().cast<double>();
	}
}