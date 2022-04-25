#include "inverse_kinematics.h"
#include "end_effectors_objective_and_gradient.h"
#include "projected_gradient_descent.h"
#include "levenberg_marquardt.h"
void ik::inverse_kinematics(Skeleton& skeleton, Eigen::VectorXi& b, Eigen::VectorXd& xb0, int max_iter, double tol, IK_SOLVER solver , double mu)
{
	//expose this		
	// Gather initial angles
	Eigen::VectorXd A(skeleton.size() * 3);
	for (int si = 0; si < skeleton.size(); si++)
	{
		A.block(si * 3, 0, 3, 1) = skeleton[si].xzx;
	}
	std::function<double(const Eigen::VectorXd&)> f;
	std::function<Eigen::VectorXd(const Eigen::VectorXd&)> grad_f;
	std::function<Eigen::MatrixXd(const Eigen::VectorXd&)> dxda;
	std::function<void(Eigen::VectorXd&)> proj_z;
	end_effectors_objective_and_gradient(skeleton, b, xb0, f, grad_f, proj_z);

	// Optimize angles
	if (solver == IK_SOLVER::LM)
	{
		end_effectors_kinematic_jacobian(skeleton, b, dxda);
		levenberg_marquardt(f, grad_f, dxda, proj_z, max_iter, tol, mu, A);
	}
	else
		projected_gradient_descent(f, grad_f, proj_z, max_iter, tol, A);
	//projected_gradient_descent(f, grad_f, proj_z, max_iters, tol, A);
	// Distribute optimized angles
	for (int si = 0; si < skeleton.size(); si++)
	{
		skeleton[si].xzx = A.block(si * 3, 0, 3, 1);
	}
};