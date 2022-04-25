#include "NullRigControllerScriptedConstrainedBall.h"
#include <igl/slice.h>
/*
Rig controller class for a null rig that simply sets and identifies motion of vertices centered at the C.O.M of X inside a ball of radius r
*/

NullRigControllerScriptedConstrainedBall::NullRigControllerScriptedConstrainedBall(Eigen::MatrixXd& X, double r) :  r(r)
{
	Eigen::RowVector3d center = X.colwise().mean();

	for (int i = 0; i < X.rows(); i++)
	{
		const double d = (X.row(i) - center).squaredNorm();
		if (d < r * r)
		{
			bI.conservativeResize(bI.rows() + 1);
			bI(bI.rows() - 1) = i;
		}
	}

	igl::slice(X, bI, 1, bc0);

}

void NullRigControllerScriptedConstrainedBall::scripted_motion(int step, Eigen::MatrixXd& u)
{
	u.setZero();
	u.resizeLike(bc0);

	u.col(0) = 0.1* sin(step * 0.1) * Eigen::VectorXd::Ones(bc0.rows());

}


