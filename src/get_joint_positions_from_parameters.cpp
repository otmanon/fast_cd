#include "get_joint_positions_from_parameters.h"
#include "matrix4f_from_parameters.h"
void get_joint_positions_from_parameters(Eigen::VectorXd& p0, Eigen::MatrixXd& joints)
{
	int num_b = p0.rows() / 12;
	joints.resize(num_b, 3);

	Eigen::Matrix4d A;
	Eigen::RowVector4d l;
	for (int i = 0; i < num_b; i++)
	{
		A = matrix4f_from_parameters(p0, i).cast<double>();
		joints.row(i) = A.block(0, 3, 3, 1).transpose();
	}
}
