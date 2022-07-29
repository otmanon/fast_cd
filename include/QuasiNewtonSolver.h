#pragma once
#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <functional>
#include <igl/min_quad_with_fixed.h>
//This q
class QuasiNewtonSolver
{
public:

	QuasiNewtonSolver();

	void precompute(const Eigen::SparseMatrix<double>& Q);
	
	void precompute_with_equality_constraints(const Eigen::SparseMatrix<double>& Q, const Eigen::SparseMatrix<double>& S);

	Eigen::VectorXd solve(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
		std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f, bool do_line_search = true, bool to_convergence = false, double max_iters = 10, double  convergence_threshold = 1e-6); //this hess_f should just evaluate the laplacian, in the rest state... remains constant for now and doesn't change with z

	//solves for z using a newton iterative scheme with constraints bc.
	Eigen::VectorXd solve_with_equality_constraints(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
		std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f, const Eigen::VectorXd& bc0, bool do_line_search = true, bool to_convergence = false, double max_iters = 10, double convergence_threshold = 1e-6);
public:


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