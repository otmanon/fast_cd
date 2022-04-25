#include "complementary_equality_constraint.h"
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <momentum_leaking_matrix.h>
#include <iostream>
void complementary_equality_constraint(const Eigen::MatrixXd& X, const Eigen::MatrixXi& T,
	const Eigen::SparseMatrix<double>& J, Eigen::SparseMatrix<double>& Aeq)
{
	Eigen::SparseMatrix<double> M, D;
	igl::massmatrix(X, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
	M = igl::repdiag(M, 3);
	

	momentum_leaking_matrix(X, T, fast_cd::MOMENTUM_LEAK_DIFFUSION, D);
	D = igl::repdiag(D, 3);
	
	if (J.rows() > 0)
	{
		//augment hessian with complementary constraint jacobians.
		Aeq = J.transpose() * (M * D);
	}
	else {
		std::cout << "Rig Jacobian has no columns! assuming null rig (no complementary constraint)" << std::endl;
		Aeq.resize(0, M.rows());
		Aeq.setZero();
	}


}
