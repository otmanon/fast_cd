#include "deformation_gradient.h"
#include "igl/slice.h"
#include <Eigen\Dense>
void deformation_gradient(Eigen::MatrixXd& X, Eigen::MatrixXi& T, Eigen::MatrixXd& U, Eigen::MatrixXd& F)
{
	Eigen::MatrixXd tet_X, tet_P;

    Eigen::Matrix3d PB, XB, XBi;

	Eigen::MatrixXd P = X + U;
	Eigen::MatrixXd B = Eigen::MatrixXd(4, 3);
	B <<	-1,-1,-1,
			1, 0, 0,
			0, 1, 0,
			0, 0, 1;

	F.resize(3 * T.rows(), 3);
	for (int i = 0; i < T.rows(); i++)
	{
		const Eigen::VectorXi ind = T.row(i);
		igl::slice(X, ind, 1, tet_X);
		igl::slice(P, ind, 1, tet_P);
		tet_X = tet_X.transpose().eval();
		tet_P = tet_P.transpose().eval();
		
		PB = tet_P * B;
		XB = tet_X * B;
		XBi = XB.inverse();
		
		F.block(3 * i, 0, 3, 3) = (PB * XBi).transpose();
	}
}