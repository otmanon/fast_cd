#pragma once
#include <igl/massmatrix.h>
#include <igl/repdiag.h>

#include "momentum_leaking_matrix.h"

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <string>
using namespace std;
using namespace Eigen;
struct cd_sim_params
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

	// linear rig jacobian! #V*dim x #bones*dim*dim+1
	SparseMatrix<double> J;

	//whether to activate inertia or not. if false, will add a bit of regularization to quadratic term
	bool do_inertia;

	//Extra user defined quadratic term (like for mass springs)
	SparseMatrix<double> Q;

	string sim_constraint_type;

	cd_sim_params() {

		mu = VectorXd::Zero(0.0);
		lambda = VectorXd::Zero(0.0);
		h = 0.0;
		do_inertia = false;
	};
	cd_sim_params(const MatrixXd& X, const MatrixXi& T,
		const SparseMatrix<double>& J, double mu, double lambda, double h,
		bool do_inertia, string sim_constraint_type = "none") :
		cd_sim_params(X, T, mu, lambda, h, do_inertia, sim_constraint_type)
	{

		this->J = J;
	
		if (sim_constraint_type == "cd")
		{
			SparseMatrix<double>  M;
			igl::massmatrix(X, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
			this->Aeq = (J.transpose() * igl::repdiag(M, 3));
		}
		else if (sim_constraint_type == "cd_momentum_leak")
		{
			SparseMatrix<double> D, M;
			momentum_leaking_matrix(X, T, fast_cd::MOMENTUM_LEAK_DIFFUSION, D);
			igl::massmatrix(X, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
			this->Aeq = (J.transpose() * igl::repdiag(M, 3) * igl::repdiag(D, 3));

		}
		else
		{
			this->Aeq.resize(0, X.rows() * X.cols());
		}
	}

	cd_sim_params(const MatrixXd& X, const MatrixXi& T, 
		double mu, double lambda, double h,
		bool do_inertia, string sim_constraint_type = "none")
	{
		this->X = X;
		this->T = T;
		this->mu = mu * VectorXd::Ones(T.rows());
		this->lambda = VectorXd::Zero(T.rows());
		this->h = h;
		this->invh2 = 1.0 / (h * h);
		this->do_inertia = do_inertia;
		this->sim_constraint_type = sim_constraint_type;
		Q.resize(X.rows() * 3, X.rows() * 3);
		Q.setZero();
	}

};