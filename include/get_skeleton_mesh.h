#pragma once
#include <Eigen/Core>
/*
From absolute world rig parameters p, build your skeleton mesh.
*/
void get_skeleton_mesh(float thickness, Eigen::VectorXd& p, Eigen::VectorXd& bl, Eigen::MatrixXd& renderV, Eigen::MatrixXi& renderF, Eigen::MatrixXd& renderC);