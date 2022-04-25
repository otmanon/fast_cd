#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
/*
Build prefactorized matrices that let you easily calculated a flattened covariance matrix from a given flattened displacement field. 

Once these matrices are calculates we should be able to get F from u^c like so. 

f = H + K u
F = mat(f) 
Note that F is an integrated deformation gradient over all tets. To get the constant version just divide by volume of the tet.
Assumes here that u^c is column order flattened, and same with F
*/
void deformation_gradient_from_u_prefactorized_matrices(const Eigen::MatrixXd& X, const Eigen::MatrixXi& T, Eigen::VectorXd& H, Eigen::SparseMatrix<double>& K, Eigen::SparseMatrix<double>& M);