#include "interweaving_matrix.h"
#include "igl/repmat.h"
#include "igl/colon.h"

void interweaving_matrix(int rows, int cols, SparseMatrix<double>& S)
{
	int flattened_dim = rows * cols;

	std::vector<Triplet<double>> tripletList;
	tripletList.reserve(flattened_dim);
	for (int i_c = 0; i_c < flattened_dim; i_c++)
	{
		//column flattened index to matrix indices
		int j = i_c / rows;
		int i = i_c - j * rows;

		//matrix index to row flattened index
		int i_r = i * cols + j;

		tripletList.emplace_back(Triplet<double>(i_c, i_r, 1));
	}
	
	S.resize(flattened_dim, flattened_dim);
	S.setFromTriplets(tripletList.begin(), tripletList.end());

}