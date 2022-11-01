#include "kron.h"

void kron(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, Eigen::MatrixXd& C)
{
	int na = B.rows();
	int ma = B.cols();
	C.resize(na * A.rows(), ma * A.cols());
	C.setZero();
	for (int i = 0; i < A.rows(); i++)
	{
		for (int j = 0; j < A.cols(); j++)
		{
			int ci = i * na; int cj = j * ma;
			C.block(ci, cj, na, ma) = B * A(i, j);
		}
	}
}