#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>

using namespace Eigen;

/*
	Matrix that maps from row-order flattening to column-order flattening and vice versa.
	Imagine we have V, a rows x cols matrix that we wish to flatten.

	Then :
	col_flattening(V) = S * row_flattening(V)
	row_flattening(V) = S.transpose() * col_flattening(V)

	*/
void interweaving_matrix(int rows, int cols, SparseMatrix<double>& S);