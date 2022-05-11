#pragma once
#include <Eigen/Core>

void get_rig_pinned_motion(Eigen::VectorXd& p_rest, Eigen::VectorXd& p_rel, Eigen::MatrixXd& X, std::vector<Eigen::VectorXi> bI, Eigen::VectorXd& bc);
