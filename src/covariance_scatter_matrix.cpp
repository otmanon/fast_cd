#include "covariance_scatter_matrix.h"
#include "igl/cotmatrix_entries.h"
#include "igl/arap_linear_block.h"
#include "igl/cat.h"
#include "interweaving_matrix.h"
void covariance_scatter_matrix(MatrixXd& V, MatrixXi& T, SparseMatrix<double>& K)
{
	SparseMatrix<double> KX, KY, KZ;
	igl::arap_linear_block(V, T, 0, igl::ARAP_ENERGY_TYPE_ELEMENTS, KX);
	igl::arap_linear_block(V, T, 1, igl::ARAP_ENERGY_TYPE_ELEMENTS, KY);
	igl::arap_linear_block(V, T, 2, igl::ARAP_ENERGY_TYPE_ELEMENTS, KZ);

	K = igl::cat(2, KX, igl::cat(2, KY, KZ));

	//this assumes a col-wise flattening of rotation. we will usually be given as input row-wise
	//convert row-wise to col-wise using interweaving matrix
	SparseMatrix<double> S;
	interweaving_matrix(KX.cols(), V.cols(), S);

	K = K * S;

}