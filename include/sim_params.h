#pragma once
#include "ympr_to_lame.h"

#include <igl/massmatrix.h>
#include <igl/repdiag.h>

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <string>
using namespace std;
using namespace Eigen;
struct sim_params
{
public:
	//geometry
	MatrixXd X;

	//tet indices
	MatrixXi T;

	//timestep
	double h;
	double invh2;

	// first lame parameter per tet
	VectorXd mu;

	//second lame parameter per tet
	VectorXd lambda;

	//Equality constraints given as input from user
	SparseMatrix<double> Aeq;

	//whether to activate inertia or not. if false, will add a bit of regularization to quadratic term
	bool do_inertia;

	//Extra user defined quadratic term (like for mass springs)
	SparseMatrix<double> Q;

	string sim_constraint_type;

	sim_params() {
	};

	sim_params(const MatrixXd& X, const MatrixXi& T,
		double ym, double pr, double h,
		bool do_inertia)
	{
		this->X = X;
		this->T = T;

		double mu, lambda;
		ympr_to_lame(ym, pr, mu,lambda);
		
		this->mu = mu * VectorXd::Ones(T.rows());
		this->lambda = lambda * VectorXd::Ones(T.rows());
		this->h = h;
		this->invh2 = 1.0 / (h * h);
		this->do_inertia = do_inertia;
		this->sim_constraint_type = sim_constraint_type;
		Q.resize(X.rows() * 3, X.rows() * 3);
		Q.setZero();
		this->Aeq.resize(0, X.rows() * X.cols());
	}


	void set_equality_constraint(SparseMatrix<double>& Aeq)
	{
		assert(Aeq.cols() == X.rows() * X.cols());
		this->Aeq = Aeq;
	}
};