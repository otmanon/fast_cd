#pragma once
#include <Eigen/Core>

//Eigen::MatrixXd scale_and_center_geometry(Eigen::MatrixXd& V, const double h = 1, const Eigen::RowVector3d& c = Eigen::RowVector3d(0, 0, 0));

/// <summary>
///  Scales and Centers geometry V so that it has height h and has it's mean position be c.
/// </summary>
/// <param name="V"></param> Geometry, list of vertex positions (3D)
/// <param name="so"></param>   returned so, how much V needs to be scaled by to have height h
/// <param name="to"></param>  returned to, how much V needs to be translated to have height c
/// <param name="h"></param>   desired height, default=1
/// <param name="c"></param>  desired center of mass, default=[0 0 0]
/// <returns></returns>
Eigen::MatrixXd scale_and_center_geometry(Eigen::MatrixXd& V, const  double h, const Eigen::RowVector3d c , double& so, Eigen::RowVector3d& to);
