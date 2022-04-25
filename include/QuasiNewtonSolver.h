#pragma once
#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <functional>
#include <igl\min_quad_with_fixed.h>
//This q
class QuasiNewtonSolver
{
public:

	QuasiNewtonSolver(const int max_iter, const double tolerance, double alpha = 1, int max_iter_line_search = 1);

	void precompute(const Eigen::SparseMatrix<double>& Q);
	
	void precompute_with_equality_constraints(const Eigen::SparseMatrix<double>& Q, const Eigen::SparseMatrix<double>& S);

	Eigen::VectorXd solve(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
		std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f); //this hess_f should just evaluate the laplacian, in the rest state... remains constant for now and doesn't change with z

	//solves for z using a newton iterative scheme with constraints bc.
	Eigen::VectorXd solve_with_equality_constraints(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
		std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f, const Eigen::VectorXd& bc0);
public:
	int max_iters;
	int max_iter_line_search;

	double alpha;

	double tolerance;
	Eigen::SparseMatrix<double>  S; //equality constraint matrix 3c x 3n
	
	Eigen::SparseLU<Eigen::SparseMatrix<double>> lu_precomp;
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt_precomp;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt_precomp;

	//lhs matrix of constrained system to project z0 onto constraint space
	Eigen::MatrixXd CCT;
	Eigen::LLT<Eigen::MatrixXd> llt_proj;
};