#pragma once

#include "fast_sim_params.h"
#include "read_clusters_from_cache.h"
#include "read_modes_from_cache.h"
#include "lbs_jacobian.h"
#include "momentum_leaking_matrix.h"

#include <igl/massmatrix.h>

using namespace Eigen;
using namespace std;


struct fast_cd_sim_params : fast_sim_params
{
	SparseMatrix<double> J;

	fast_cd_sim_params() {};


	fast_cd_sim_params(const MatrixXd& X, const MatrixXi& T,
		const MatrixXd& B, const VectorXi& l,
		const SparseMatrix<double>& J, double ym,
		double pr, double h,
		bool do_inertia) : fast_sim_params(X, T, B, l, ym, pr, h, do_inertia)
	{
		this->J = J;	
	}

	/*
	Contains all the parameters required to build a
	fast Complementary Dynamics simulator.

	Inputs:
		X - n x 3 vertex geometry
		T - T x 4 tet indices
		B - 3n x m subspace matrix
		l - T x 1 clustering labels for each tet
		J - 3n x sol_p rig jacobian. If sim_constraint_type = none, this also
								doubles as the equality constraint matrix
		mu - (double) first lame parameter
		lambda - (double) second lame parameter
		do_inertia - (bool) whether or not sim should have inertia
				(if no, then adds Tik. regularizer to laplacian
		sim_constraint_type - (string) "none" or "cd" or "cd_momentum_leak"
	*/
	fast_cd_sim_params(const MatrixXd& X, const MatrixXi& T,
		const MatrixXd& B, const VectorXi& l,
		const SparseMatrix<double>& J, double mu,
		double lambda, double h,
		bool do_inertia, string sim_constraint_type)
	{
		this->X = X;
		this->T = T;
		this->J = J;
		this->mu = mu * VectorXd::Ones(T.rows());
		this->lambda = VectorXd::Zero(T.rows());
		this->h = h;
		this->invh2 = 1.0 / (h * h);
		this->do_inertia = do_inertia;
		this->sim_constraint_type = sim_constraint_type;

		this->B = B;
		this->labels = l;
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
		if (sim_constraint_type == "custom")
		{
			printf("Please provide a custom  \
					#c x dimN equality constraint  \
					matrix to fast_cd_sim_params like so \n \
					fast_cd_sim_params(X, T, B, l, J, \
					 Aeq, mu, lambda, h, do_inertia, \
					sim_constraint_type) \
					");
		}
		else
		{
			this->Aeq.resize(0, X.rows() * X.cols());
		}
	}


	/*
	Contains all the parameters required to build a
	fast Complementary Dynamics simulator.

	Inputs:
		X - n x 3 vertex geometry
		T - T x 4 tet indices
		B - 3n x m subspace matrix
		l - T x 1 clustering labels for each tet
		J - 3n x sol_p rig jacobian. If sim_constraint_type = none, this also
								doubles as the equality constraint matrix
		Aeq -c x 3n linear equality constraint matrix
		mu - (double) first lame parameter
		lambda - (double) second lame parameter
		do_inertia - (bool) whether or not sim should have inertia
				(if no, then adds Tik. regularizer to laplacian

	*/
	fast_cd_sim_params(const MatrixXd& X, const MatrixXi& T,
		const MatrixXd& B, const VectorXi& l,
		const SparseMatrix<double>& J, const SparseMatrix<double>& Aeq,
		double mu, double lambda, double h,
		bool do_inertia, string sim_constraint_type = "none")
	{
		this->X = X;
		this->T = T;
		this->J = J;
		this->mu = mu * VectorXd::Ones(T.rows());
		this->lambda = VectorXd::Zero(T.rows());
		this->h = h;
		this->invh2 = 1.0 / (h * h);
		this->do_inertia = do_inertia;
		this->sim_constraint_type = "custom";
		this->B = B;
		this->labels = l;

		this->Aeq = Aeq;
		assert(Aeq.cols() == X.rows() * X.cols());

	}



	void set_equality_constraint(SparseMatrix<double>& Aeq)
	{
		this->Aeq = Aeq;
	}


};