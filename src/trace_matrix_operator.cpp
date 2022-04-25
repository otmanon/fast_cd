#include "trace_matrix_operator.h"

void trace_matrix_operator(const int k, const int dim, Eigen::SparseMatrix<double>& Tr)
{
	Tr.resize(k, k * dim * dim);
	
	std::vector<Eigen::Triplet<double >> tripletList;
	tripletList.reserve(dim * k); // 3 entries per row
	for (int i = 0; i < k; i++)
	{
		for (int d = 0; d < dim; d++)
		{
			tripletList.emplace_back(i, dim*i + dim*k*d + d, 1);
		}
	}
	Tr.setFromTriplets(tripletList.begin(), tripletList.end());
}