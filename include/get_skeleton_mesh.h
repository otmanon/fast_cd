#pragma once
#include <Eigen/Core>
/*
From absolute world rig parameters sol_p, build your skeleton mesh.
*/
void get_skeleton_mesh( const float thickness,  const Eigen::VectorXd& p,  const Eigen::VectorXd& bl, Eigen::MatrixXd& renderV, Eigen::MatrixXi& renderF, Eigen::MatrixXd& renderC);

void get_skeleton_mesh( const float thickness,  const Eigen::VectorXd& p, const  Eigen::VectorXd& p0, const Eigen::VectorXd& bl, Eigen::MatrixXd& renderV, Eigen::MatrixXi& renderF, Eigen::MatrixXd& renderC);