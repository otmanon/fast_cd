#ifndef IGL_SPLIT_COMPONENTS_H
#define IGL_SPLIT_COMPONENTS_H

#include <Eigen/Sparse>
#include <Eigen/Core>
#include <vector>

namespace igl
{
  // Given a binary labeling (I) of nodes in a graph with adjacency (A), label
  // the connected components of the true regions
  //
  // Inputs:
  //   I  #I list of binary values 
  //   A  #I by #I adjacency matrix
  // Outputs:
  //   S  #S list of lists of labeled nodes so that S[i][j] = k indicates that
  //     node k is the jth node in component i
  void split_components(
    const Eigen::Array<bool,Eigen::Dynamic,1> & I,
    const Eigen::SparseMatrix<bool> & A,
    std::vector<std::vector<int> > & S);
  // Given an integer labeling (I) of nodes in a graph with adjacency (A),
  // output a new labeling where each label is a maximal connected sub-component of an
  // original label.
  //
  // Inputs:
  //   I  #I list of labels
  //   A  #I by #I adjacency matrix
  // Outputs:
  //   J  #I new labels
  void split_components(
    const Eigen::VectorXi & I,
    const Eigen::SparseMatrix<bool> & A,
    Eigen::VectorXi & J);

  // Inputs:
  //   I  #I list of labels
  //   F  #F by 3 list of triangle indices into rows of I
  // Outputs:
  //   J  #I new labels
  void split_components(
    const Eigen::VectorXi & I,
    const Eigen::MatrixXi & F,
    Eigen::VectorXi & J);
  // Given a list of scalar functions (columns of B) over a mesh with faces (F),
  // output a new list of functions where each is a maximal connected (in-terms
  // of non-zero scalar values) sub-component of an original function.
  //
  // Inputs:
  //   B  #V by #B list of input functions
  //   F  #F by 3 list of triangle indices into rows of V
  // Outputs:
  //   C  #V by #C≥#B list of output functions
  void split_components(
    const Eigen::MatrixXd & B,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd & C);
}

#ifndef IGL_STATIC_LIBRARY
#  include "split_components.h"
#endif

#endif 
