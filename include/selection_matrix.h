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


/*
Creats a selection matrix S,
that selects out rows corresponding to vertices of index
bI of an n*dim x 1 vector

Inputs :
   bI - bI x 1 indices to select
   n - number of full vertices
		(used to determine the number of columns of S)
   dim - if bI indexes vertices, but each vertex has dim components,
		then bI is first expanded into bIx, bIy and bIz and the selection
		matrix indexes each of the components individually. Only dim=2 and 3 are supported
Outputs :
	S - c*dim x n*dim selection matrix
*/
void selection_matrix(const VectorXi& bI, int n, int dim, SparseMatrix<double>& S)
{
	
//	bIy = bIx.array() + n;
//	VectorXi I;
//	if (dim == 2)
//		I = igl::cat(1, bIx, bIy);
//	else if (dim == 3)
//	{
//VectorXi bIx = bI, bIy = bI, bIz = bI;
//		bIz = bIy.array() + n;
//		I = igl::cat(1, bIx, igl::cat(1, bIy, bIz));
//	}
//	else {
//		printf("dimension of selection matrix not supported! \n");
//		return;
//	}
//	int nd = n * dim;
	selection_matrix(bI, n, S);
	SparseMatrix<double>S2;
	igl::repdiag(S, dim, S2);
	S = S2;
}