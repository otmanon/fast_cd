#pragma once
#include <Eigen/Core>

Eigen::VectorXd compute_bbw_weights(Eigen::VectorXd& p, Eigen::MatrixXd& V, Eigen::MatrixXi& T); // , Eigen::MatrixXd);