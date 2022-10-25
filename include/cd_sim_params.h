#pragma once
#include <igl/massmatrix.h>
#include <igl/repdiag.h>

#include <Eigen/Sparse>
#include <Eigen/Core>

using namespace Eigen;
struct cd_sim_params
{
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

	cd_sim_params() {};
	cd_sim_params(MatrixXd& X, MatrixXi& T,
		SparseMatrix<double>& J, double mu, double lambda, double h,
		bool do_inertia)
	{
		this->X = X;
		this->T = T;
		this->mu = mu * VectorXd::Ones(T.rows());
		this->lambda = VectorXd::Zero(T.rows());
		this->h = h;
		this->invh2 = 1.0 / (h * h);
		this->do_inertia = do_inertia;

		Q.resize(X.rows() * 3, X.rows() * 3);
		Q.setZero();

		this->J = J;

		Eigen::SparseMatrix<double> M;
		massmatrix(X, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
		M = igl::repdiag(M, 3);
		this->Aeq = J.transpose().eval() * M; // no D-matrix here!
	}

};