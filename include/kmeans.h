#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>
namespace igl
{
  // Inputs:
  //   P  #P by dim list of input points
  //   W  #P list of positive weights
  //   k  number of clusters
  //   nearest_neighbor  function so that nearest_neighbor(C,P,I,sqrD) computes
  //     the nearest neighbor in C for each point in P, storing indices in I and
  //     squared distances in sqrD.
  // Outputs
  //   C  k by dim list of cluster centroids
  //   I  #P list of indices into rows of C
   void kmeans(
    const Eigen::MatrixXd & P,
    const Eigen::VectorXd & W,
    const int k,
    const std::function<void(const Eigen::MatrixXd&,const Eigen::MatrixXd&,Eigen::VectorXi&,Eigen::VectorXd&)> & nearest_neighbor,
    Eigen::MatrixXd & C,
    Eigen::VectorXi & I);
  // Inputs:
  //   P  #P by dim list of input points
  //   W  #P list of positive weights
  //   k  number of clusters
  // Outputs
  //   C  k by dim list of cluster centroids
  //   I  #P list of indices into rows of C
  void kmeans(
    const Eigen::MatrixXd & P,
    const Eigen::VectorXd & W,
    const int k,
    Eigen::MatrixXd & C,
    Eigen::VectorXi & I);
  // Inputs:
  //   P  #P by dim list of input points
  //   k  number of clusters
  // Outputs
  //   C  k by dim list of cluster centroids
  //   I  #P list of indices into rows of C
   void kmeans(
    const Eigen::MatrixXd & P,
    const int k,
    Eigen::MatrixXd & C,
    Eigen::VectorXi & I);
}
