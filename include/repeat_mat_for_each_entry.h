#pragma once
#include "interweaving_matrix.h"

#include <igl/repdiag.h>
#include <Eigen/Sparse>
using namespace Eigen;
/*
Repeat mat for each entry
given a sparse matrix like this
G = [v1 0 0]
	[0 v2 0]
Repeats the values so that they get multiplied against corresponding entries in a flattened matrix


A - SparseMatrix whose entries are to be repeated
n - rows of matrix whose entries need to be multiplied by each value of A
m - number of cols of matrix whose entreis need to be multiplied by each value of A */
SparseMatrix<double> repeat_mat_for_each_entry(SparseMatrix<double>& A, int n, int m)
{
	SparseMatrix<double>A_tmp = igl::repdiag(A, n);
	SparseMatrix<double> S_cols, S_rows;
	interweaving_matrix(A.cols(), n, S_cols);
	interweaving_matrix(A.rows(), n, S_rows);
	A_tmp = S_rows.transpose() * A_tmp * S_cols;
	SparseMatrix<double> Ar= igl::repdiag(A_tmp, m);
	return Ar;
}