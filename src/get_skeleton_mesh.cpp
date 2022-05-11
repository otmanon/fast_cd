#include "get_skeleton_mesh.h"

#include "matrix4f_from_parameters.h"


/*
From absolute world rig parameters p, build your skeleton mesh.
*/
void get_skeleton_mesh(float thickness, Eigen::VectorXd& p, Eigen::VectorXd& bl, Eigen::MatrixXd& renderV, Eigen::MatrixXi& renderF, Eigen::MatrixXd& renderC)
{


	Eigen::MatrixXd BV(5, 3);
	BV <<
		0, -1, -1,
		0, 1, -1,
		0, 1, 1,
		0, -1, 1,
		1, 0, 0;
	BV.rightCols(2) *= thickness;
	Eigen::MatrixXi BF(6, 3);
	BF <<
		0, 2, 1,
		0, 3, 2,
		0, 1, 4,
		1, 2, 4,
		2, 3, 4,
		3, 0, 4;
	Eigen::MatrixXd BC(6, 3);
	Eigen::RowVector3d red(1, 0, 0), green(0, 1, 0), blue(0, 0, 1);
	BC <<
		1 - red.array(),
		1 - red.array(),
		1 - blue.array(),
		green,
		blue,
		1 - green.array();
	int num_bones = p.rows() / 12;
	//for (int b = 0; b < num_b; b++)
	//{
	//	if (bone_list[b].parent_index >= 0) { num_bones++; }
	//}
	renderV.resize(BV.rows() * num_bones, 3);
	renderF.resize(BF.rows() * num_bones, 3);
	renderC.resize(BF.rows() * num_bones, 3);
	{
		int k = 0;
		for (int b = 0; b < num_bones; b++)
		{

			//Eigen::Matrix4f parent_A_rest = matrix4f_from_parameters(p_rest, );

			const double len = bl[b];
			Eigen::MatrixXd BVk(BV.rows(), 3);
			Eigen::Matrix4d A = matrix4f_from_parameters(p, b).cast<double>();
			for (int v = 0; v < BV.rows(); v++)
			{
				const Eigen::Vector4d p =
					(A * Eigen::Vector4d(len * BV(v, 0), BV(v, 1), BV(v, 2), 1));
				BVk.row(v) = p.topRows(3).transpose();
			}
			renderV.block(k * BV.rows(), 0, BV.rows(), 3) = BVk;
			renderF.block(k * BF.rows(), 0, BF.rows(), 3) = (BF.array() + k * BV.rows()).matrix();
			renderC.block(k * BC.rows(), 0, BC.rows(), 3) = BC;
			k++;
		}
	}
}