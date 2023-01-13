#pragma once
#include <Eigen/Core>

Eigen::MatrixXd compute_bbw_weights(const Eigen::VectorXd& p, const Eigen::MatrixXd& V, const Eigen::MatrixXi& T); // , Eigen::MatrixXd);