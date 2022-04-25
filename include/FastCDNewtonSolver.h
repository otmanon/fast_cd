#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <functional>
#include <igl/min_quad_with_fixed.h>
//This q
class FastCDNewtonSolver
{
public:
	FastCDNewtonSolver(const int max_iter, const double tolerance, const Eigen::MatrixXd& Q, double alpha = 1, int max_iter_line_search = 1);

	FastCDNewtonSolver(const int max_iter, const double tolerance, double alpha = 1, int max_iter_line_search = 1);

	void precompute(const Eigen::MatrixXd& Q);

	void precompute_with_constraints(const Eigen::MatrixXd& Q, const Eigen::MatrixXd& Qeq);

	Eigen::VectorXd solve(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
		std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f); //this hess_f should just evaluate the laplacian, in the rest state... remains constant for now and doesn't change with z
	
	Eigen::VectorXd solve_with_constraints(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
		std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f, const Eigen::VectorXd& uc_bc);
public:
	int max_iters;
	int max_iter_line_search;

	Eigen::MatrixXd A, Qeq, Q;
	double alpha;

	double tolerance;

	Eigen::LLT<Eigen::MatrixXd> precomp;          //prefactorization for linear equality constraint solvers
	Eigen::LDLT<Eigen::MatrixXd> ldlt_precomp;
};