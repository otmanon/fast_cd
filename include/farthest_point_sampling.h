#pragma once


#include <Eigen/Core>

void farthest_point_sampling(const Eigen::MatrixXd& V, const Eigen::MatrixXi& T, const int num_samples, Eigen::VectorXi& I, Eigen::MatrixXd& P);