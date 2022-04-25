#pragma once
#include <Eigen/Core>
#include <Eigen/Sparse>

void grouping_matrix_from_clusters(Eigen::VectorXi& I, Eigen::SparseMatrix<double>& G);