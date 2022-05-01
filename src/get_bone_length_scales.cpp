#include "get_bone_length_scales.h"
#include "get_tip_positions_from_parameters.h"
#include "get_joint_positions_from_parameters.h"
double get_bone_length_scales(Eigen::VectorXd& p, Eigen::VectorXi& pI)
{
	Eigen::MatrixXd C;
	get_joint_positions_from_parameters(p, C);

	int num_b = p.rows() / 12;

	Eigen::VectorXd lengths;
	//look through skeleton hierarchy
	for (int i = 0; i < num_b; i++)
	{
		if (pI(i) >= 0)
		{
			const Eigen::RowVector3d disp = C.row(i) - C.row(pI(i));
			const double dist = disp.squaredNorm();
			if (dist)
			{
				lengths.conservativeResize(lengths.rows() + 1);
				lengths(lengths.rows() - 1) = dist;
			}
		}
	}

	//get_tip_positions_from_parameters(p_tmp, pl, tips); //damn Im done... we actually can't do this because pl is just the bone length... need to loop through each bone, go to its parent, calc the diff.
	double ls = lengths.sum() / lengths.rows();
	return ls;
}