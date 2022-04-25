#pragma once
#include <Eigen/Core>


/*
for a mesh with rest X with tets T, and a displacement field at each vertex U, calculate a list of stacked deformation gradients F for each tet.

Following Sifakes notes.
*/
void deformation_gradient(Eigen::MatrixXd& X, Eigen::MatrixXi& T, Eigen::MatrixXd& U, Eigen::MatrixXd& F);