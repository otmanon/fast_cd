#include "compute_handle_positions_from_rig_parameters.h"

//these are row order flattened
void compute_handle_positions_from_parameters(Eigen::VectorXd& p, Eigen::MatrixXd& V)
{
	int c = 12;
	int num_p = p.rows() / c;
	V.resize(num_p, 3);
	for (int i = 0; i < num_p; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			V(i, j) = p(12 * i + j * 4 + 3);
		}
	}
}
