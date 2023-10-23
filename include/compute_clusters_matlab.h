//#pragma once
//#include <Eigen/Core>
//
//void compute_clusters_matlab(const Eigen::MatrixXi& T, const Eigen::MatrixXd& B, const Eigen::VectorXd& L, const int num_clusters, const int num_feature_modes, Eigen::VectorXi& labels, Eigen::MatrixXd& C);
//
//
///*
//Computes the clusters above with an initial guess of labels0, Where labels 0 may have less rows than labels.
//
//T - |T|x4 tet simplex
//B - modes
//L - eigenvalue for each mode
//num_clusters - number of clusters to compute
//num_feature modes - number of feature modes to compute with (diminishing returns on this, default should be 10)
//labels0 - |T|x1 list of cluster indeces
//
//*/
//void compute_clusters_matlab_initguess(const Eigen::MatrixXi& T, const Eigen::MatrixXd& B, const Eigen::MatrixXd& L, const int num_clusters, const int num_feature_modes, const Eigen::MatrixXd& C0, Eigen::VectorXi& labels, Eigen::MatrixXd& C);