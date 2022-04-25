#pragma once
#include <eigen/core>
#include <eigen/sparse>
using namespace Eigen;

/*
Builds the covariance scatter matrix K, where V^T K gives you a stacked list of covariance matrices.
*/
void covariance_scatter_matrix(MatrixXd& V, MatrixXi& T, SparseMatrix<double>& K);