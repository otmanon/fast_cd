#include "complementary_equality_constraint.h"
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <momentum_leaking_matrix.h>
#include <iostream>

#include "arap_hessian.h"
#include "remove_zero_cols.h"
void complementary_equality_constraint(const Eigen::MatrixXd& X, const Eigen::MatrixXi& T,
	const Eigen::SparseMatrix<double>& J, Eigen::SparseMatrix<double>& Aeq, double dt, std::string metric_type, double alpha)
{
	Eigen::SparseMatrix<double> M, D, K;
	
	if (metric_type == "momentum")
	{
		igl::massmatrix(X, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
		M = igl::repdiag(M, 3);
		momentum_leaking_matrix(X, T, fast_cd::MOMENTUM_LEAK_DIFFUSION, D, dt);
		D = igl::repdiag(D, 3);
		K = M * D;
	}
	else if (metric_type == "elastic")
	{
		Eigen::SparseMatrix<double> H;
		arap_hessian(X, T, H);
		K = H;
	}
	else if (metric_type == "mix")
	{

		igl::massmatrix(X, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
		M = igl::repdiag(M, 3);
	//	momentum_leaking_matrix(X, T, fast_cd::MOMENTUM_LEAK_DIFFUSION, D, dt);
	//	D = igl::repdiag(D, 3);

		Eigen::SparseMatrix<double> H;
		arap_hessian(X, T, H);

		K = (M +  alpha * H);
	}

	if (J.rows() > 0)
	{

	
		//augment hessian with complementary constraint jacobians.
		Aeq = J.transpose() * K;
		
		Eigen::VectorXi I, CI;
		Eigen::SparseMatrix<double> A1, A2;
		//A1 = Aeq.transpose();
		remove_zero_cols(Aeq, A2, I, CI);
		Aeq = A2;
	}
	else {
		std::cout << "Rig Jacobian has no columns! assuming null rig (no complementary constraint)" << std::endl;
		Aeq.resize(0, K.rows());
		Aeq.setZero();
	}


}
