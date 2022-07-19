#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
void complementary_equality_constraint(const Eigen::MatrixXd& X, const Eigen::MatrixXi& T,
	const Eigen::SparseMatrix<double>& J, Eigen::SparseMatrix<double>& Aeq, double dt = 1e-6, std::string metric_type="momentum", double alpha=1.0);
