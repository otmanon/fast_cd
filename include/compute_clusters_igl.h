#pragma once
#include <Eigen/Core>

void compute_clusters_igl(const Eigen::MatrixXi& T, const Eigen::MatrixXd& B, const Eigen::VectorXd& L, const int num_clusters, const int num_feature_modes, Eigen::VectorXi& labels, Eigen::MatrixXd& C);