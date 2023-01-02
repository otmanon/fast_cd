#include "compute_bbw_weights.h"
#include "rig_parameters.h"
#include "get_joint_positions_from_parameters.h"

#include <igl/point_mesh_squared_distance.h>
#include <igl/bbw.h>
#include <igl/colon.h>
#include <igl/unique_rows.h>

/*
For a rig with rest pose parameters sol_p, attached to a volumetric mesh V, T, compute bounded biharmonic weights
by finding the closest vertex point to each rig joint in sol_p, and calling bbw.h
*/
Eigen::MatrixXd compute_bbw_weights(Eigen::VectorXd& p, Eigen::MatrixXd& V, Eigen::MatrixXi& T)
{
	Eigen::MatrixXd joint_positions;
	get_translations_from_rig_parameters(p, joint_positions);
	
	Eigen::MatrixXi VI = igl::colon<int>(0, V.rows() - 1);

	Eigen::VectorXd sqrD;
	Eigen::VectorXi closestI;
	Eigen::MatrixXd closest_point;
	igl::point_mesh_squared_distance(joint_positions, V, VI, sqrD, closestI, closest_point);
	
	Eigen::VectorXi _n1_, _n2_;
	Eigen::VectorXi uniqueI;
	igl::unique_rows(closestI, uniqueI, _n1_, _n2_);
	
	assert(uniqueI.rows() == closestI.rows() && "Rig joints are too close together and got mapped to the same vertex!");
	Eigen::MatrixXd bc = Eigen::MatrixXd::Identity(uniqueI.rows(), uniqueI.rows());
	Eigen::MatrixXd W;
	igl::BBWData data;
	
	igl::bbw(V, T, closestI, bc, data, W);
	// normalize at the end.

	for (int i = 0; i < W.rows(); i++)
	{
		W.row(i) = W.row(i) / W.row(i).sum();
	}
	return W;
}