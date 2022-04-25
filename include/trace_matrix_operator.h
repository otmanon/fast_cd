#pragma once
#include <Eigen\Sparse>
/*
	Assume we have k  dim x dim matrices stacked up in a (k*dim x dim) matrix list. Flatten that matrix list to a vector (k*dim*dim x 1)
	using row order flattening. This function returns a sparse matrix that mimics the trace operator when multiplied against that flattened vector
*/
void trace_matrix_operator(const int k, const int dim, Eigen::SparseMatrix<double>& Tr);