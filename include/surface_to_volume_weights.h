#pragma once
#include <Eigen/Core>

Eigen::MatrixXd surface_to_volume_weights(  Eigen::MatrixXd& surfaceW,   Eigen::VectorXi& bI,   Eigen::MatrixXd& X,   Eigen::MatrixXi& T);


Eigen::MatrixXd surface_to_volume_weights( const Eigen::MatrixXd& surfaceW, const  Eigen::MatrixXd& surfaceV,  const Eigen::MatrixXd& X, const  Eigen::MatrixXi& T);