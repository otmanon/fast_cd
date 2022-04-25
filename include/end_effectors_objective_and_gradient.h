#pragma once
#include "Skeleton.h"
void end_effectors_objective_and_gradient(
    const Skeleton& bone_list,
    const Eigen::VectorXi& b,
    const Eigen::VectorXd& xb0,
    std::function<double(const Eigen::VectorXd&)>& f,
    std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& grad_f,
    std::function<void(Eigen::VectorXd&)>& proj_z);


void end_effectors_kinematic_jacobian(
    const Skeleton& bone_list,
    const Eigen::VectorXi& b,
    std::function<Eigen::MatrixXd(const Eigen::VectorXd&)>& dxda);