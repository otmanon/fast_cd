//#pragma once
//#include "arap_hessian.h"
//#include "displacement_modes.h"
//#include "selection_matrix.h"
//#include "skinning_modes.h"
//#include "linear_elasticity_hessian.h"
//
//#include <Eigen/Core>
//#include <Eigen/Sparse>
//#include <iostream>
//#include <igl/massmatrix.h>
//#include <igl/repdiag.h>
//#include <igl/colon.h>
//#include <igl/cat.h>
//
//
//
//using namespace Eigen;
//using namespace igl;
///*
//V->n x 3 geometry
//T->f x 4 tets
//W->n x b weights
//J ->  12bx3n lbs jacobian, may act as constraint on modes
//mode_type->string specifiying mode type("displacement" or "skinning")
//num_modes -> int, number of modes to search for/compute. If cache have less modes than this, will recompute
//(Optional)
//debug -> if true, will write some debugging files to .DMAT in directory output_dir
//output_dir -> where to write debug files
//Output
//B -> either 3n x num_modes for displacement modes, or n x num_modes/12 for skinning modes
//L -> either num_modes vector or num_modes/12 vector of eigenvalues
//Ws -> If we picked the skinning subspace, returns the weights here. otherwise leaves it empty
//*/
//using namespace std;
//void get_modes(const MatrixXd& V, const MatrixXi& T, const SparseMatrix<double>& J, std::string mode_type,
// int num_modes, MatrixXd& B, VectorXd& L, MatrixXd& Ws, bool debug=false, string output_dir="")
//{
//	assert(J.cols() == V.rows() * V.cols() && "Rig Jacobian must have as many cols as there are full space degrees of freedom 3#V");
//	SparseMatrix<double> M;
//	massmatrix(V, T, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
//
//	SparseMatrix<double>  H;
//
//
//	if (mode_type == "displacement")
//	{ //these are arap hessian modes
//		arap_hessian(V, T, H);
//		displacement_modes(H, M, J, num_modes, B, L, debug, output_dir);
//	}
//	else if (mode_type == "skinning")
//	{
//		arap_hessian(V, T, H);
//		skinning_modes( V, H, M, J , num_modes, B, Ws,  L, debug, output_dir);
//	}
//	else{
//		std::cout << "mode type " << mode_type << " is not yet supported" << std::endl;
//		exit(0);
//	}
//
//}