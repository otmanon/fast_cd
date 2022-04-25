#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
namespace igl
{
  void old_kmeans(
    const Eigen::MatrixXd & P,
    const Eigen::VectorXd & W,
    const int k,
    Eigen::MatrixXd & C,
    Eigen::VectorXi & I);
  void old_kmeans(
    const Eigen::MatrixXd & P,
    const int k,
    Eigen::MatrixXd & C,
    Eigen::VectorXi & I);
}


