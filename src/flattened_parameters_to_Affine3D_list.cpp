#include "flattened_parameters_to_Affine3D_list.h"
#include "matrix4f_from_parameters.h"

void flattened_parameters_to_Affine3d_list(Eigen::VectorXd& p, std::vector<Eigen::Affine3d>& A)
{
	int num_b = p.rows() / 12;
	A.reserve(num_b);
	Eigen::MatrixXd T; Eigen::Affine3d TG;
	for (int i = 0; i < num_b; i++)
	{
		T = matrix4f_from_parameters(p, i).cast<double>();
		TG.matrix() = T;
		A.emplace_back(TG);
	}
}