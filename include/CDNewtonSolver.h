#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <functional>
#include <igl/min_quad_with_fixed.h>
//This q
class CDNewtonSolver
{

public:
	CDNewtonSolver(const int max_iter, const double tolerance,const Eigen::SparseMatrix<double>& Q,  const Eigen::SparseMatrix<double>& Qeq, double alpha = 1, int max_iter_line_search = 1);
	
	CDNewtonSolver(const int max_iter, const double tolerance, double alpha = 1, int max_iter_line_search = 1);

	void precompute(const Eigen::SparseMatrix<double>& Q, const Eigen::SparseMatrix<double>& Qeq);
	void precompute_with_constraints(const Eigen::SparseMatrix<double>& Q, const Eigen::SparseMatrix<double>& Aeq, const Eigen::SparseMatrix<double>& Seq);
	
	Eigen::VectorXd solve(const Eigen::VectorXd& z,  std::function<double(const Eigen::VectorXd&)>& f,
		std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f); //this hess_f should just evaluate the laplacian, in the rest state... remains constant for now and doesn't change with z
	
	//solves for z using a newton iterative scheme with constraints bc.
	Eigen::VectorXd solve_with_constraints(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
		std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f, const Eigen::VectorXd& bc0);
public:
	int max_iters;
	int max_iter_line_search;

	double alpha;

	double tolerance;
	Eigen::SparseMatrix<double> Qeq, CompEq;
	igl::min_quad_with_fixed_data<double> precomp;          //prefactorization for linear equality constraint solvers
	igl::min_quad_with_fixed_data<double> precomp_constraints;


	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt_precomp;
};