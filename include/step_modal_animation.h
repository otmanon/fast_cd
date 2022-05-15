#pragma once
#include <Eigen/Core>
#include <igl/PI.h>
#include <iostream>
Eigen::MatrixXd step_modal_animation(Eigen::MatrixXd& B, int step, int modei, double period, double scale)
{
    Eigen::VectorXd B_flat = B.col(modei);
    double mod = scale * (1.0 - cos((2 * igl::PI / (period)) * step));
    std::cout << scale << std::endl;
    Eigen::VectorXd u = mod * B_flat;
    Eigen::MatrixXd U = Eigen::Map<Eigen::MatrixXd>(u.data(), u.rows() / 3, 3);
    return U;
}
