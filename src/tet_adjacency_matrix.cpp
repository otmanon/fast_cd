#include "tet_adjacency_matrix.h"
#include <igl/tet_tet_adjacency.h>

void tet_adjacency_matrix_face_share(const Eigen::MatrixXi& T, Eigen::SparseMatrix<bool, 0, int>& A)
{
	assert(T.cols() == 4 && "tets have four vertices");
	Eigen::MatrixXi tet_tet_adj, _n;
	igl::tet_tet_adjacency(T, tet_tet_adj, _n);

	std::vector<Eigen::Triplet<int>> tripletList;

	for (int i = 0; i < T.rows(); i++)
	{
		for (int j = 0; j < T.cols(); j++)
		{
			if (tet_tet_adj(i, j) >= 0) //if not adjacent, will be a -1
			{
				tripletList.push_back(Eigen::Triplet<int>(i, tet_tet_adj(i, j), 1));
			}
		}
	}

	A.resize(T.rows(), T.rows());
	A.setFromTriplets(tripletList.begin(), tripletList.end());
}



void tet_adjacency_matrix_vert_share(const Eigen::MatrixXi& T, Eigen::SparseMatrix<bool, 0, int>& A)
{
	assert(T.cols() == 4 && "tets have four vertices");

	std::vector<Eigen::Triplet<double>> V2TIJV;
	for (int i = 0; i < T.rows(); i++)
	{
		for (int j = 0; j < T.cols(); j++)
		{
			V2TIJV.emplace_back(T(i, j), i, 1);
		}
	}
	Eigen::SparseMatrix<double> V2T(T.maxCoeff() + 1, T.rows());
	V2T.setFromTriplets(V2TIJV.begin(), V2TIJV.end());
	Eigen::SparseMatrix<double> Adj = V2T.transpose() * V2T; //

	for (int i = 0; i < T.rows(); i++)
	{
		Adj.coeffRef(i, i) = 0; //set htis to 0
	}

	A = Adj.cast<bool>();
}

// careful, A will contain diagonal entries and non-zeros might be >1 