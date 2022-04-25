#pragma once
#include "Eigen/Dense"
#include "Eigen/SparseCore"
void arap_hessian(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    Eigen::SparseMatrix<double>& H);