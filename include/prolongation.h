#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
/*
Builds an interpolation matrix that maps from a coarse mesh to a fine version of the same mesh. 

returns W such that
X_fine = W X_coarse 

If a fine element is not inside the coarse mesh, we project it to the closest point. 
*/
void prolongation(Eigen::MatrixXd& X_fine, Eigen::MatrixXd& X_coarse, Eigen::MatrixXi& T_coarse, Eigen::SparseMatrix<double>& W);

/*
Builds an interpolation matrix that maps from a coarse mesh to a fine version of the same mesh.

returns W such that
X_fine = W X_coarse

If a fine element is not inside the coarse mesh, we project it to the closest point.

Input: 
X_fine -> Geometry of fine tet mesh  nf x 3
X_coarse -> Geometry of coarse tet mesh nc x 3
T_coarse ->  tet indices of coarse tet mesh F x 4. We only map the surface of this coarse mesh to the fine mesh.
Output:
W -> the coarse to fine interpolation matrix s.t. X_fine = W X_coarse
WI -> the coarse to fine interpolation matrix with 1's instead of interpolation weights
*/
void prolongation(Eigen::MatrixXd& X_fine, Eigen::MatrixXd& X_coarse, Eigen::MatrixXi& T_coarse, Eigen::SparseMatrix<double>& W, Eigen::SparseMatrix<double>& WI);