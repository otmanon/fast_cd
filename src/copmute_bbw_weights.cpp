#include "compute_bbw_weights.h"
#include <igl/bbw.h>
#include <igl/point_mesh_squared_distance.h>
#include "get_joint_positions_from_parameters.h"
#include <igl/colon.h>
#include <igl/unique_rows.h>
Eigen::MatrixXd compute_bbw_weights(Eigen::VectorXd& p, Eigen::MatrixXd& V, Eigen::MatrixXi& T)
{
	Eigen::MatrixXd joint_positions;
	get_joint_positions_from_parameters(p, joint_positions);
	
	Eigen::MatrixXi VI = igl::colon<int>(0, V.rows() - 1);

	Eigen::VectorXd sqrD;
	Eigen::VectorXi closestI;
	Eigen::MatrixXd closest_point;
	igl::point_mesh_squared_distance(joint_positions, V, VI, sqrD, closestI, closest_point);
	
	Eigen::VectorXi _n1_, _n2_;
	Eigen::VectorXi uniqueI;
	igl::unique_rows(closestI, uniqueI, _n1_, _n2_);
	Eigen::MatrixXd bc = Eigen::MatrixXd::Identity(uniqueI.rows(), uniqueI.rows());
	Eigen::MatrixXd W;
	igl::BBWData data;
	
	igl::bbw(V, T, uniqueI, bc, data, W);

	return W;
}