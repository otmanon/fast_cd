#include "kron.h"

void kron(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, Eigen::MatrixXd& C)
{
	int na = A.rows();
	int ma = A.cols();
	C.resize(na * B.rows(), ma * B.cols());
	C.setZero();
	for (int i = 0; i < B.rows(); i++)
	{
		for (int j = 0; j < B.cols(); j++)
		{
			int ci = i * na; int cj = j * ma;
			C.block(ci, cj, na, ma) = A * B(i, j);
		}
	}
}