#pragma once

#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Core>

using namespace Eigen;
/*
Creats a selection matrix S,
that selects out rows bI of an nx1 vector

Inputs : 
   bI - bI x 1 indices to select
   n - number of full vertices 
		(used to determine the number of columns of S)
Outputs :
	S - c x n selection matrix
*/
void selection_matrix(const VectorXi& bI, int n, SparseMatrix<double>& S)
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