#include "levenberg_marquardt.h"
#include "line_search.h"
#include <Eigen/Dense>
void levenberg_marquardt(
	const std::function<double(const Eigen::VectorXd&)>& f,
	const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& dfda,
	const std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& dxda,		
	const std::function<void(Eigen::VectorXd&)>& proj_z,
	const int max_iters,
	const double tol,
	const double mu,					//parameter that blends between gradient descent and gauss-newton
	Eigen::VectorXd& z)
{

	Eigen::VectorXd dz;

	for (int i = 0; i < max_iters; i++)
	{
		Eigen::VectorXd grad = dfda(z);
		if (grad.squaredNorm() < tol) return;
		
		Eigen::MatrixXd J = dxda(z);
		Eigen::MatrixXd I = Eigen::MatrixXd::Identity(z.rows(), z.rows());
		Eigen::MatrixXd A = J.transpose() * J + mu * I;
		dz = A.llt().solve(-grad);
	//	dz = A.ldlt().solve_with_constraints(grad); // negative grad? check to ensure factor of 1/2 isn't msising.
		//dz = -grad;

		double alpha = line_search_pgd(f, proj_z, z, dz, 2000);
		z += alpha * dz;
		proj_z(z);
	}
}
