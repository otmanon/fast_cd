#pragma once
#include <Eigen/Core>
namespace igl
{
  // Jacobi-style online updates. These are worst-case costly because a single
  // iteration through all points could be O(n k log k).
  // 
  // Inputs:
  //   P  #P by dim list of input points
  //   W  #P list of positive weights
  //   k  number of clusters
  //   C  k by dim list of cluster centroids initial guesses
  //   I  #P list of indices into rows of C
  //   sqrD  #P list of squared distances to corresponding cluster center
  //   Wcount  k list of weighted counts
  //   obj  current scalar objective function 
  // Outputs:
  //   C,I,sqr,Wcount,obj  updated to converged values
  void kmeans_jacobi(
    const Eigen::MatrixXd & P,
    const Eigen::VectorXd & W,
    const int k,
    Eigen::MatrixXd & C,
    Eigen::VectorXi & I,
    Eigen::VectorXd & sqrD,
    Eigen::VectorXd & Wcount,
    double & obj);
}



