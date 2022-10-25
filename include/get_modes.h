#pragma once
#include "arap_hessian.h"
#include "displacement_modes.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
using namespace Eigen;
/*
V->n x 3 geometry
T->f x 4 tets
W->n x b weights
J -> 3n x 12b lbs jacobian, may act as constraint on modes
mode_type->string specifiying mode type("displacement" or "skinning")
num_modes -> int, number of modes to search for/compute. If cache have less modes than this, will recompute

Output
B -> either 3n x num_modes for displacement modes, or n x num_modes/12 for skinning modes
L -> either num_modes vector or num_modes/12 vector of eigenvalues
*/
void get_modes(MatrixXd& V, MatrixXi& T, MatrixXd& W, SparseMatrix<double>& J, std::string mode_type,
 int num_modes, MatrixXd& B, VectorXd& L)
{

	if (mode_type == "displacement")
	{ //these are arap hessian modes
		SparseMatrix<double>  H, M;
		massmatrix(V, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
		M = igl::repdiag(M, 3);
		arap_hessian(V, T, H);
		displacement_modes(H, M, J.transpose() * M, num_modes, B, L);
	}
	else{
		std::cout << "mode type " << mode_type << " is not yet supported" << std::endl;
		exit(0);
	}
}