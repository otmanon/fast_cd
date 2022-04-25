#pragma once
#include "Skeleton.h"
Eigen::VectorXd transformed_tips(
    const Skeleton& bone_list,
    const Eigen::VectorXi& b);
