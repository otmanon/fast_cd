#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
void fast_complementary_dynamics_constrained_geivp(
	const Eigen::MatrixXd& X, const Eigen::MatrixXi& T, const Eigen::SparseMatrix<double>& J,
	Eigen::SparseMatrix<double>& H, Eigen::SparseMatrix<double>& M, double dt=1e-6, std::string metric_type="momentum", double alpha=1.0);