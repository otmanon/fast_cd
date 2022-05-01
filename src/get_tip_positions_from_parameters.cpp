#include "get_tip_positions_from_parameters.h"
void get_tip_positions_from_parameters(Eigen::VectorXd& p0, Eigen::VectorXd& bl, Eigen::MatrixXd& tips)
{
	int num_b = p0.rows() / 12;
	tips.resize(num_b, 3);

	Eigen::Matrix4d A;
	Eigen::RowVector4d l;
	for (int i = 0; i < num_b; i++)
	{
		A = matrix4f_from_parameters(p0, i).cast<double>();
		l = Eigen::RowVector4d(bl(i), 0, 0, 1);
		tips.row(i) = (A * l.transpose()).topRows(3);
	}
}
