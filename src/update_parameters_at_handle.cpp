#include "update_parameters_at_handle.h"
void update_parameters_at_handle(Eigen::VectorXd& p, Eigen::Matrix4f T, int ind)
{
	int c = 12;
	int num_p = p.rows() / c;
	Eigen::MatrixXd A = T.topRows(3).transpose().cast<double>();
	p.middleRows(12 * ind, 12) = Eigen::Map<Eigen::VectorXd>(A.data(), 12);
}