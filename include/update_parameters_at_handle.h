#pragma once
#include "Eigen/Core"

/*
Given a list of row-order flattened 12x1 affine handle matrices p, update the entries corresponding to the ind's handle with the affine transformation matrix T
*/
void update_parameters_at_handle(Eigen::VectorXd& p, Eigen::Matrix4f T, int ind);