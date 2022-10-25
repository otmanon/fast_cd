#pragma once
#include "cd_sim_params.h"

struct fast_cd_sim_params : cd_sim_params
{
	MatrixXd B;
	VectorXi labels;
	fast_cd_sim_params() {};
	fast_cd_sim_params(MatrixXd& X, MatrixXi& T, MatrixXd& B, VectorXi& l,
		SparseMatrix<double>& J, double mu, double lambda, double h,
		bool do_inertia) : cd_sim_params(X, T, J, mu, lambda, h, do_inertia)
	{
		this->B = B;
		this->labels = l;

		//usually no equality constraint.
		this->Aeq.resize(0, X.rows() * 3);
	}
};