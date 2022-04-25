#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <functional>
#include <igl/min_quad_with_fixed.h>
class CDLocalGlobalSolver
{
public:
	CDLocalGlobalSolver() {};

	CDLocalGlobalSolver(int max_iters, double tol, const Eigen::SparseMatrix<double>& Q, const Eigen::SparseMatrix<double>& Qeq);

	CDLocalGlobalSolver(int max_iters, double tol);

	void precompute(const Eigen::SparseMatrix<double>& Q, const Eigen::SparseMatrix<double>& Qeq);

	Eigen::VectorXd solve(const Eigen::VectorXd& z,
		std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& local_step,
		std::function<Eigen::VectorXd(const Eigen::MatrixXd&, const igl::min_quad_with_fixed_data<double>&)>& global_step);

public:

	int max_iters;
	double tol;

	igl::min_quad_with_fixed_data<double> precomp;          //prefactorization for linear equality constraint solvers
};