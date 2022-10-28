#pragma once
#include "arap_hessian.h"
#include "displacement_modes.h"
#include "selection_matrix.h"
#include "skinning_modes.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <igl/massmatrix.h>
#include <igl/repdiag.h>
#include <igl/colon.h>
#include <igl/cat.h>
using namespace Eigen;
using namespace igl;
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
Ws -> If we picked the skinning subspace, returns the weights here. otherwise leaves it empty
*/
void get_modes(MatrixXd& V, MatrixXi& T, MatrixXd& W, SparseMatrix<double>& J, std::string mode_type,
 int num_modes, MatrixXd& B, VectorXd& L, MatrixXd& Ws)
{
	SparseMatrix<double> M, M2;
	massmatrix(V, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
	M2 = igl::repdiag(M, 3);
	SparseMatrix<double>  H;
	arap_hessian(V, T, H);

	if (mode_type == "displacement")
	{ //these are arap hessian modes
		displacement_modes(H, M, J.transpose() * M2, num_modes, B, L);
	}
	else if (mode_type == "skinning")
	{
		skinning_modes( V, H, M, J.transpose() * M2, num_modes, B, Ws,  L);	
	}
	else{
		std::cout << "mode type " << mode_type << " is not yet supported" << std::endl;
		exit(0);
	}
}