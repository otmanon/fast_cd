#pragma once

#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Core>

using namespace Eigen;
void selection_matrix(VectorXi& bI, int n, SparseMatrix<double>& S)
{
	S.resize(bI.rows(), n);
	S.setZero();
	std::vector<Triplet<double>> tripletList;
	tripletList.reserve(bI.rows());
	for (int i = 0; i < bI.rows(); i++)
	{
		tripletList.emplace_back(i, bI(i), 1);
	}
	S.setFromTriplets(tripletList.begin(), tripletList.end());
}