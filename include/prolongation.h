#pragma once
#include <Eigen\Dense>
#include <Eigen\Sparse>
/*
Builds an interpolation matrix that maps from a coarse matrix to a fine version of the same matrix. 

returns W such that
X_fine = W X_coarse 

If a fine element is not inside the coarse mesh, we project it to the closest point. 
*/
void prolongation(Eigen::MatrixXd& X_fine, Eigen::MatrixXd& X_coarse, Eigen::MatrixXi& T_coarse, Eigen::SparseMatrix<double>& W);

/*
Builds an interpolation matrix that maps from a coarse matrix to a fine version of the same matrix.

returns W such that
X_fine = W X_coarse

If a fine element is not inside the coarse mesh, we project it to the closest point.
*/
void prolongation(Eigen::MatrixXd& X_fine, Eigen::MatrixXd& X_coarse, Eigen::MatrixXi& T_coarse, Eigen::SparseMatrix<double>& W, Eigen::SparseMatrix<double>& WI);