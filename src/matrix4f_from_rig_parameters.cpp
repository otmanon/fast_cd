#include "matrix4f_from_parameters.h"

Eigen::MatrixX4f matrix4f_from_parameters(Eigen::VectorXd& p, int ind)
{
	int c = 12;
	int num_p = p.rows() / c;

	Eigen::VectorXf a = p.middleRows(12 * ind, 12).cast<float>();
	Eigen::Matrix4f A = Eigen::Matrix4f::Identity();



	Eigen::MatrixXf T = Eigen::Map<Eigen::MatrixXf>(a.data(), 4, 3);
	A.topRows(3) = T.transpose();
	//for (int i = 0; i < 3; i++)
	//{
	//	Eigen::VectorXd row = p.middleRows(ind * 4 + i * 4 * (num_p),4);
	//	T.row(i) = row.transpose().cast<float>();
//
	//}

	return A;
}
