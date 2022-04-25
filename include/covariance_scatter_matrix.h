#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
using namespace Eigen;

/*
Builds the covariance scatter matrix K, where V^T K gives you a stacked list of covariance matrices.
*/
void covariance_scatter_matrix(MatrixXd& V, MatrixXi& T, SparseMatrix<double>& K);