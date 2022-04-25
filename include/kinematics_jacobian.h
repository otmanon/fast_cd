#pragma once
#include "Skeleton.h"
void kinematics_jacobian(
    const Skeleton& bone_list,
    const Eigen::VectorXi& b,
    Eigen::MatrixXd& J);