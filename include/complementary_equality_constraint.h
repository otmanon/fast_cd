#pragma once
#include <Eigen\Core>
#include <Eigen\Sparse>
void complementary_equality_constraint(const Eigen::MatrixXd& X, const Eigen::MatrixXi& T,
	const Eigen::SparseMatrix<double>& J, Eigen::SparseMatrix<double>& Aeq);
