
#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>


/*
Computes the r smallest eigenmodes of the generalized eigenvalue problem:
A U = L B U

U is a matrix of stacked eigenvectors. L are the r eigenvalues. still need to expose/return L, but that'as barely relevant for this project
*/
void compute_modes_spectra(Eigen::SparseMatrix<double>&A, Eigen::SparseMatrix<double>&B, int r, Eigen::MatrixXd & U, Eigen::VectorXd & S) {
	printf("Computing Modes with Spectra Is Not Implemented in This Branch!");
};