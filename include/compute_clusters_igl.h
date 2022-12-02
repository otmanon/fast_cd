#pragma once
#include <Eigen/Core>

void compute_clusters_igl(const Eigen::MatrixXi& T, const Eigen::MatrixXd& B, const Eigen::VectorXd& L, const int num_clusters, const int num_feature_modes, Eigen::VectorXi& labels, Eigen::MatrixXd& C);

/*
Computes clusters over our tet mesh, with tet indices T, using displacement faetures.

Inputs:
    T -> |T|x4 tet indices
    B -> dim|V|xnum_modes subspace basis
    L -> num_modes x num_modes eigenvalues sorted in ascending order with the columns of B
    num_clusters -> (int) number of clusters to compute
    num_feature_modes -> (int) number of columns of our subspace basis to use. Can truncate this to small number for efficiency
    split_components -> (bool) whether to split clusters into separate components. Some clusters may group tets that are far away together, this will amend that
*/
void compute_clusters_displacement_features(const Eigen::MatrixXi& T, const Eigen::MatrixXd& B, const Eigen::VectorXd& L, const int num_clusters, const int num_feature_modes, Eigen::VectorXi& labels, Eigen::MatrixXd& C, bool split_components = false);

/*
Computes clusters over our tet mesh, with tet indices T, using weight faetures.

Inputs:
    T -> |T|x4 tet indices
    W -> |V|xnum_modes skinning subspace weights
    L -> num_modes x num_modes eigenvalues sorted in ascending order with the columns of B
    num_clusters -> (int) number of clusters to compute
    num_feature_modes -> (int) number of columns of our subspace basis to use. Can truncate this to small number for efficiency
    split_components -> (bool) whether to split clusters into separate components. Some clusters may group tets that are far away together, this will amend that
*/
void compute_clusters_weight_features(const Eigen::MatrixXi& T, const Eigen::MatrixXd& W, const Eigen::VectorXd& L, const int num_clusters, const int num_feature_modes, Eigen::VectorXi& labels, Eigen::MatrixXd& C, bool split_components = false);