#pragma once
#include "corot_sim_params.h"

#include <igl/massmatrix.h>

using namespace Eigen;
struct fast_corot_sim_params : corot_sim_params
{
	MatrixXd B;
	VectorXi labels;

	MatrixXd AeqB;
	fast_corot_sim_params() {};

	fast_corot_sim_params(const MatrixXd& X, const MatrixXi& T,
		const MatrixXd& B, const VectorXi& l,
		double ym, double pr, double h,
		bool do_inertia) : corot_sim_params(X, T,  ym, pr, h, do_inertia)
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