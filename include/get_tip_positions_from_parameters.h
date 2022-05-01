#pragma once
#include <Eigen/Core>

#include "matrix4f_from_parameters.h"
void get_tip_positions_from_parameters(Eigen::VectorXd& p0, Eigen::VectorXd& bl, Eigen::MatrixXd& tips);