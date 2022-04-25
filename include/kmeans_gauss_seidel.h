#pragma once
#include <Eigen/Core>
namespace igl
{
  // Gauss-Seidel-style updates. These are potentially efficient if
  // nearest_neighbor is O((k+n) log k).
  // 
  // Inputs:
  //   P  #P by dim list of **unique** input points
  //   W  #P list of positive weights
  //   k  number of clusters
  //   nearest_neighbor  function so that nearest_neighbor(C,P,I,sqrD) computes
  //     the nearest neighbor in C for each point in P, storing indices in I and
  //     squared distances in sqrD.
  //   C  k by dim list of cluster centroids initial guesses
  // Outputs
  //   C  updated
  //   I  #P list of indices into rows of C
  //   sqrD  #P list of squared distances to corresponding cluster center
  //   Wcount  k list of weighted counts
  //   obj  current scalar objective function 
  // Returns true if converged.
     bool kmeans_gauss_seidel(
    const Eigen::MatrixXd & P,
    const Eigen::VectorXd & W,
    const int k,
    const std::function<void(const Eigen::MatrixXd&,const Eigen::MatrixXd&,Eigen::VectorXi&,Eigen::VectorXd&)> & nearest_neighbor,
    Eigen::MatrixXd & C,
    Eigen::VectorXi & I,
    Eigen::VectorXd & sqrD,
    Eigen::VectorXd & Wcount,
    double & obj);
}


