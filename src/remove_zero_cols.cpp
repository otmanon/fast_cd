#include "remove_zero_cols.h"
#include <igl/slice.h>
void remove_zero_cols(Eigen::SparseMatrix<double>& A, Eigen::SparseMatrix<double>& B, Eigen::VectorXi& I, Eigen::VectorXi& CI)
{
	//sum through each column
	Eigen::VectorXd sum = A * Eigen::VectorXd::Ones(A.cols());
	double thresh = 1e-12;
	int n = 0;

	CI = -Eigen::VectorXi::Ones(A.cols());
	for (int row = 0; row < sum.rows(); row++)
	{
		if (sum(row) > thresh) // great we make the cut
		{
			I.conservativeResize(I.rows() + 1);
			I(n) = row;
			CI(row) = n;
			n += 1;
		}
	}

	igl::slice(A, I, 1, B);
	B = A;
}