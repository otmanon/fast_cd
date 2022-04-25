#include <sparse_diag.h>
SparseMatrix<double> sparse_diag(VectorXd& diag)
{
	std::vector<Triplet<double>> tripletList;
	tripletList.reserve(diag.rows());
	for (int i = 0; i < diag.rows(); i++)
	{
		Triplet<double> t = Triplet<double>(i, i, diag[i]);
		tripletList.emplace_back(t);
	}

	SparseMatrix<double> S = SparseMatrix<double>(diag.rows(), diag.rows());
	S.setFromTriplets(tripletList.begin(), tripletList.end());
	return S;
}