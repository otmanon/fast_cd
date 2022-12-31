#pragma once
#include "vector_gradient_operator.h"
#include "repeat_for_each_entry.h"
#include "vectorized_transpose.h"
#include "vectorized_trace.h"

#include <igl/volume.h>
#include <igl/diag.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

void linear_elasticity_hessian(MatrixXd& X, MatrixXi& T,
	double mu, double lambda, SparseMatrix<double>& H)
{
	VectorXd muv = mu * VectorXd::Ones(T.rows());
	muv = repeat_for_each_entry(muv, 3, 3);
	

	SparseMatrix<double> Mu, Lam;
	igl::diag(muv, Mu);
	
	SparseMatrix<double> K;
	vector_gradient_operator(X, T, K);

	VectorXd tet_volume, tet_volume_exp; 
	igl::volume(X, T, tet_volume);
	tet_volume_exp = repeat_for_each_entry(tet_volume, 3, 3);
	SparseMatrix<double> V;
	igl::diag(tet_volume_exp, V);


	SparseMatrix<double> Transposer;
	vectorized_transpose(3, 3, T.rows(), Transposer);


	SparseMatrix<double> shear = K.transpose() * V * Mu * K + 
		K.transpose() * V * Mu * Transposer * K;


	VectorXd lamv = lambda * VectorXd::Ones(T.rows());
	igl::diag(lamv, Lam);
	igl::diag(tet_volume, V);

	SparseMatrix<double> Trace;
	vectorized_trace(3, T.rows(), Trace);
	SparseMatrix<double> vol = K.transpose() * Trace.transpose() * Lam * V * Trace * K;
	//now to add the volume preserving parts

	H = vol + shear;
}