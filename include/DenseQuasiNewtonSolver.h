#pragma once
#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>
#include <functional>
//This q
class DenseQuasiNewtonSolver
{
public:

	DenseQuasiNewtonSolver();

	void precompute(const Eigen::MatrixXd& Q);

	void precompute_with_equality_constraints(const Eigen::MatrixXd& Q, const Eigen::MatrixXd& S);

	Eigen::VectorXd solve(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
		std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f, bool do_line_search = true, bool to_convergence = false, double max_iters=20); //this hess_f should just evaluate the laplacian, in the rest state... remains constant for now and doesn't change with z

	//solves for z using a newton iterative scheme with constraints bc.
	Eigen::VectorXd solve_with_equality_constraints(const Eigen::VectorXd& z, std::function<double(const Eigen::VectorXd&)>& f,
		std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f, const Eigen::VectorXd& bc0, bool do_line_search = true, bool to_convergence = false,double max_iters = 20);
public:
	int max_iters;
	int max_iter_line_search;

	double alpha;

	double tolerance;
	Eigen::MatrixXd  S; //equality constraint matrix 3c x 3n

	
	Eigen::LLT<Eigen::MatrixXd> llt_precomp;
	Eigen::LDLT<Eigen::MatrixXd> ldlt_precomp;
	Eigen::FullPivLU<Eigen::MatrixXd> lu_precomp;
	Eigen::LLT<Eigen::MatrixXd> llt_proj;
};