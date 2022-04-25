#include "augment_with_linear_constraints.h"
#include "igl/cat.h"


void augment_with_linear_constraints(const SparseMatrix<double>& A, VectorXd& b, const SparseMatrix<double>& C, const VectorXd& beq, SparseMatrix<double>& Q, VectorXd& f)
{
	if (beq.rows() > 0 && C.rows() > 0)
	{
		f = igl::cat(1, b, beq);
		augment_with_linear_constraints(A, C, Q);
	}
	else
	{
		f = b;
		Q = A;
	}
}


void augment_with_linear_constraints(const SparseMatrix<double>& A, const SparseMatrix<double>& C, SparseMatrix<double>& Q)
{
	if (C.rows() > 0)
	{
		SparseMatrix<double> top_row, bot_row, z, CT;
		z.resize(C.rows(), C.rows());
		z.setZero();
		CT = C.transpose();
		top_row = igl::cat(2, A,CT);
		bot_row = igl::cat(2, C, z);
		Q = igl::cat(1, top_row, bot_row);
	}
	else
	{
		Q = A;
	}
}


void augment_with_linear_constraints(const Eigen::MatrixXd& A, const Eigen::MatrixXd& C, Eigen::MatrixXd& Q)
{
	if (C.rows() > 0)
	{
		Eigen::MatrixXd top_row, bot_row, z, CT;
		z.resize(C.rows(), C.rows());
		z.setZero();
		CT = C.transpose();
		top_row = igl::cat(2, A, CT);
		bot_row = igl::cat(2, C, z);
		Q = igl::cat(1, top_row, bot_row);
	}
	else
	{
		Q = A;
	}
}