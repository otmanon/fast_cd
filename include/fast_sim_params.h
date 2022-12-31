#pragma once
#include "sim_params.h"

#include <igl/massmatrix.h>

using namespace Eigen;
struct fast_sim_params : sim_params
{
	MatrixXd B;
	VectorXi labels;

	MatrixXd AeqB;
	fast_sim_params() {};

	fast_sim_params(const MatrixXd& X, const MatrixXi& T,
		const MatrixXd& B, const VectorXi& l,
		double mu, double lambda, double h,
		bool do_inertia) : sim_params(X, T,  mu, lambda, h, do_inertia)
	{
		this->B = B;
		this->labels = l;
		this->AeqB.resize(0, B.cols());
	}
	void set_equality_constraint(MatrixXd & AeqB)
	{
		assert(AeqB.cols() == this->B.cols() && "constraints must have same dimensionality as subspace!\n");
		this->AeqB = AeqB;
	}

};