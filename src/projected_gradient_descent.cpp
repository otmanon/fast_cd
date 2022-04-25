#include "projected_gradient_descent.h"
#include "line_search.h"

void projected_gradient_descent(
	const std::function<double(const Eigen::VectorXd&)>& f,
	const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f,
	const std::function<void(Eigen::VectorXd&)>& proj_z,
	const int max_iters,
	const double tol,
	Eigen::VectorXd& z)
{

	Eigen::VectorXd dz;
	for (int i = 0; i < max_iters; i++)
	{
		Eigen::VectorXd grad = grad_f(z);
		if (grad.squaredNorm() < tol) return;
		dz = -grad;
		double alpha = line_search_pgd(f, proj_z, z, dz, 2000);
		z += alpha * dz;
		proj_z(z);
	}
	/////////////////////////////////////////////////////////////////////////////
	// Add your code here
	/////////////////////////////////////////////////////////////////////////////
}
