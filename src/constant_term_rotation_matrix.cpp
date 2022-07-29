#include "constant_term_rotation_matrix.h"
#include "interweaving_matrix.h"
void constant_term_rotation_matrix(const Eigen::MatrixXd& W,const const int u, const int v, const int b, Eigen::SparseMatrix<double>& E_uv)
{

	Eigen::SparseMatrix<double> S;

	interweaving_matrix(W.rows(), 3, S);
	std::vector<Eigen::Triplet<double>> tripletList;
	//fill out _uv in row order  for each entry of the rotation matrix
	for (int vi = 0; vi < W.rows(); vi++)
	{
		tripletList.emplace_back(3 * vi + u, 3 * vi + v, W(vi, b));
	}

	E_uv.resize(3 * W.rows(), 3 * W.rows());
	E_uv.setFromTriplets(tripletList.begin(), tripletList.end());
	E_uv = S * E_uv * S.transpose(); //everything is in column order flattening, but build E as row order flattening out of convenience...

}