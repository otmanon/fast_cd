#pragma once
#include <Eigen/Sparse>

/*
Computes the vectorized trace, a matrix that when multiplied against a flattened
stack of matrices, returns for each of those stacked matrices its stacked trace

Input:
n - number of rows in the original submatrices
d - number of submatrices stacked
Output
Tr - d x nnd vectorized trace operator
*/
using namespace Eigen;
void vectorized_trace(int n, int d, SparseMatrix<double> Tr)
{

	std::vector<Triplet<int>> tripletList;
	tripletList.reserve(n * n * d);
	for (int i = 0; i < d; i++)
	{
		//three entries for each row i.
		tripletList.emplace_back(i, d * n + i * n + 0, 1);
		tripletList.emplace_back(i, d * n + i * n + 1, 1);
		tripletList.emplace_back(i, d * n + i * n + 2, 1);
	}
	Tr.resize(d, n * n * d);
	Tr.setFromTriplets(tripletList.begin(), tripletList.end());
}