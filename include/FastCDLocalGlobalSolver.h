#pragma once
#include <Eigen\Core>
#include <Eigen\Sparse>
#include <Eigen\Dense>
#include <functional>
class FastCDLocalGlobalSolver
{
public:
	FastCDLocalGlobalSolver() {};

	FastCDLocalGlobalSolver( int max_iters, double tol, const Eigen::MatrixXd& Q);

	FastCDLocalGlobalSolver(int max_iters, double tol);

	void precompute(const Eigen::MatrixXd& Q);
	
	Eigen::VectorXd solve(const Eigen::VectorXd& z,
		std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& local_step,
		std::function<Eigen::VectorXd(const Eigen::MatrixXd&, const Eigen::LLT<Eigen::MatrixXd>&)>& global_step);

public: 

	int max_iters;
	double tol;

	Eigen::LLT<Eigen::MatrixXd> precomp;
};